"""UVM-RAD DICOM patient viewer (read-only).

Usage:
  python uvm_dataset_reader.py --patient 63
  python uvm_dataset_reader.py --patient UVM-063 --root 'D:\\...\\UVM_RAD_DATASET_FINAL'

Shows CT, available MR series and RTSTRUCT contours as Napari labels.
Images use their physical DICOM transforms (LPS) for display; no MRI-to-CT
registration is performed. RT masks are rasterized on the CT grid.
"""
from __future__ import annotations

import argparse
import inspect
import logging
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pydicom
import SimpleITK as sitk
from skimage.draw import polygon

LOG = logging.getLogger('uvm-reader')
DEFAULT_ROOT = Path(r'D:\FNO_image_segmentation\TCIA uveal - anonym\UVM_RAD_DATASET_FINAL')
MR_PREFIXES = ('MR1_', 'MR2_', 'MR3_', 'MR4_')
COLORS = {
    'GTV': '#ff3b30', 'PTV': '#ffb000', 'Eye_R': '#00bcd4',
    'Eye_L': '#1976d2', 'Lens_R': '#ac72ff', 'Lens_L': '#7943dd',
    'OpticNrv_R': '#44ca75', 'OpticNrv_L': '#007d45',
    'CiliaryBody_R': '#ff78ac', 'CiliaryBody_L': '#d81b60',
}


def normalize_patient(value: str) -> str:
    text = str(value).strip().upper()
    match = re.fullmatch(r'(?:UVM[-_])?(\d{1,3})', text)
    if not match:
        raise ValueError('Patient must be a numeric ID (e.g. 63) or UVM-063')
    return f'UVM-{int(match.group(1)):03d}'


def read_dicom(path: Path):
    return pydicom.dcmread(str(path), stop_before_pixels=True)


def discover_series(folder: Path) -> list[list[str]]:
    """Return all image series in a directory; no assumptions on filenames."""
    if not folder.is_dir():
        return []
    ids = sitk.ImageSeriesReader.GetGDCMSeriesIDs(str(folder)) or []
    series = []
    for series_id in ids:
        names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(str(folder), series_id)
        if names:
            series.append(list(names))
    return series


def load_volume(names: list[str]):
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(names)
    image = reader.Execute()
    if image.GetDimension() != 3:
        raise ValueError(f'Expected 3D image, got dimension {image.GetDimension()}')
    arr = sitk.GetArrayFromImage(image)
    if arr.ndim != 3:
        raise ValueError(f'Expected ZYX image, got shape {arr.shape}')
    ds = pydicom.dcmread(names[0], stop_before_pixels=True)
    return image, arr, ds


def affine_zyx(image: sitk.Image) -> np.ndarray:
    """Napari voxel ZYX -> patient LPS physical coordinates ordered ZYX (mm)."""
    direction = np.asarray(image.GetDirection(), dtype=np.float64).reshape(3, 3)
    spacing = np.asarray(image.GetSpacing(), dtype=np.float64)
    reverse = np.eye(3)[::-1, :]
    a = np.eye(4, dtype=np.float64)
    a[:3, :3] = reverse @ direction @ np.diag(spacing) @ reverse
    a[:3, 3] = reverse @ np.asarray(image.GetOrigin(), dtype=np.float64)
    return a


def find_rtstruct(patient_dir: Path) -> Path | None:
    folder = patient_dir / 'RTSTRUCT'
    if not folder.is_dir():
        return None
    results = []
    for candidate in sorted(folder.rglob('*')):
        if not candidate.is_file():
            continue
        try:
            ds = read_dicom(candidate)
            if str(getattr(ds, 'Modality', '')) == 'RTSTRUCT':
                results.append(candidate)
        except Exception as exc:
            LOG.warning('Skipping unrecognized RTSTRUCT candidate %s: %s', candidate, exc)
    if not results:
        return None
    if len(results) > 1:
        LOG.warning('Multiple RTSTRUCT files found. Using %s only.', results[0])
    return results[0]


def rasterize_rtstruct(rt_path: Path, ct_image: sitk.Image, ct_sop_uids: set[str]):
    """Rasterize closed planar RT contours on the CT grid in patient coordinates.

    Nonplanar/unsupported contours are skipped with warnings. A rasterized mask
    approximates physical vector contours, does not alter their source geometry.
    """
    rt = pydicom.dcmread(str(rt_path), stop_before_pixels=True)
    if str(getattr(rt, 'Modality', '')) != 'RTSTRUCT':
        raise ValueError(f'Not an RTSTRUCT: {rt_path}')
    size_x, size_y, size_z = ct_image.GetSize()
    by_number = {int(r.ROINumber): str(r.ROIName) for r in getattr(rt, 'StructureSetROISequence', [])}
    masks = {}
    stats = defaultdict(lambda: {'contours': 0, 'skipped': 0, 'unmatched_refs': 0})
    for roi in getattr(rt, 'ROIContourSequence', []):
        roi_num = int(roi.ReferencedROINumber)
        roi_name = by_number.get(roi_num, f'ROI_{roi_num}')
        mask = np.zeros((size_z, size_y, size_x), dtype=bool)
        for contour in getattr(roi, 'ContourSequence', []):
            stats[roi_name]['contours'] += 1
            geo_type = str(getattr(contour, 'ContourGeometricType', '')).upper()
            if geo_type not in ('CLOSED_PLANAR', 'CLOSEDPLANAR_XOR'):
                stats[roi_name]['skipped'] += 1
                LOG.warning('%s: unsupported contour geometric type %s', roi_name, geo_type)
                continue
            refs = [str(im.ReferencedSOPInstanceUID) for im in getattr(contour, 'ContourImageSequence', [])]
            if refs and not set(refs).intersection(ct_sop_uids):
                stats[roi_name]['unmatched_refs'] += 1
                LOG.warning('%s: contour references image outside CT series; not rasterizing', roi_name)
                continue
            points = np.asarray(getattr(contour, 'ContourData', []), dtype=np.float64)
            if points.size < 9 or points.size % 3 != 0:
                stats[roi_name]['skipped'] += 1
                continue
            points = points.reshape(-1, 3)
            coords = np.asarray([
                ct_image.TransformPhysicalPointToContinuousIndex(tuple(map(float, p)))
                for p in points
            ])  # x,y,z continuous voxel indices
            k = int(np.rint(np.median(coords[:, 2])))
            max_plane_error = float(np.max(np.abs(coords[:, 2] - k)))
            if k < 0 or k >= size_z or max_plane_error > 0.35:
                stats[roi_name]['skipped'] += 1
                LOG.warning('%s: contour outside CT plane (slice=%d error=%.3f voxels)', roi_name, k, max_plane_error)
                continue
            rows, cols = polygon(coords[:, 1], coords[:, 0], shape=(size_y, size_x))
            if geo_type == 'CLOSEDPLANAR_XOR':
                mask[k, rows, cols] ^= True
            else:
                mask[k, rows, cols] = True
        if roi_name in masks:
            LOG.warning('Duplicate ROI name %s in RTSTRUCT: merging masks', roi_name)
            masks[roi_name] |= mask
        else:
            masks[roi_name] = mask
    return masks, stats


def volume_frame_uids(files: list[str]) -> tuple[set[str], set[str]]:
    sop_uids, frame_uids = set(), set()
    for file in files:
        ds = read_dicom(Path(file))
        sop_uids.add(str(getattr(ds, 'SOPInstanceUID', '')))
        frame_uids.add(str(getattr(ds, 'FrameOfReferenceUID', '')))
    return sop_uids, frame_uids


def launch(root: Path, patient: str, no_rtstruct: bool = False) -> None:
    patient_dir = root / patient
    if not patient_dir.is_dir():
        raise FileNotFoundError(f'Patient folder not found: {patient_dir}')
    import napari

    ct_series = discover_series(patient_dir / 'CT')
    if not ct_series:
        raise FileNotFoundError(f'No CT DICOM series in {patient_dir / "CT"}')
    ct_image, ct_data, ct_ds = load_volume(ct_series[0])
    ct_uids, ct_for = volume_frame_uids(ct_series[0])
    if len(ct_series) > 1:
        LOG.warning('More than one CT series; using the first one.')
    viewer = napari.Viewer(title=f'UVM-RAD | {patient} | CT and MR + RTSTRUCT', ndisplay=2)
    ct_layer = viewer.add_image(ct_data, name='CT', colormap='gray',
                                affine=affine_zyx(ct_image), visible=True)
    # Sensible starting window for head CT without altering source data.
    ct_layer.contrast_limits = (-150, 300)
    loaded = [('CT', ct_data.shape)]
    for folder in sorted(patient_dir.iterdir()):
        if not folder.is_dir() or not folder.name.startswith(MR_PREFIXES):
            continue
        for j, filenames in enumerate(discover_series(folder), start=1):
            try:
                image, data, ds = load_volume(filenames)
                label = folder.name if j == 1 else f'{folder.name}_{j}'
                mr_for = str(getattr(ds, 'FrameOfReferenceUID', ''))
                if ct_for and mr_for not in ct_for:
                    LOG.warning('%s has different FrameOfReferenceUID from CT. No registration is performed; avoid interpreting overlaps.', label)
                viewer.add_image(data, name=label, colormap='gray',
                                 affine=affine_zyx(image), visible=False)
                loaded.append((label, data.shape))
            except Exception as exc:
                LOG.exception('Could not load %s series %d: %s', folder, j, exc)
    if not no_rtstruct:
        rt_path = find_rtstruct(patient_dir)
        if rt_path is None:
            LOG.warning('No RTSTRUCT found in %s', patient_dir)
        else:
            masks, stats = rasterize_rtstruct(rt_path, ct_image, ct_uids)
            for roi_name, mask in masks.items():
                if not mask.any():
                    LOG.warning('%s: zero rasterized voxels (inspect contour geometry or references)', roi_name)
                    continue
                rgba = COLORS.get(roi_name, '#bbbbbb')
                # Force a deterministic color for label value 1 (foreground).
                # napari >= 0.5 uses `colormap`; older versions use `color`.
                color_kwargs = (
                    {'colormap': {1: rgba, None: '#00000000'}}
                    if 'colormap' in inspect.signature(viewer.add_labels).parameters
                    else {'color': {1: rgba, None: '#00000000'}}
                )
                viewer.add_labels(
                    mask.astype(np.uint8), name=f'ROI | {roi_name}',
                    affine=affine_zyx(ct_image), visible=True, opacity=0.5,
                    **color_kwargs,
                )
            for name, s in stats.items():
                LOG.info('%s: %d contour items; %d skipped; %d mismatched refs',
                         name, s['contours'], s['skipped'], s['unmatched_refs'])
    viewer.layers.selection.active = ct_layer
    viewer.reset_view()
    print(f'Loaded {patient}: ' + ', '.join(f'{n} {shape}' for n, shape in loaded))
    napari.run()


def main():
    parser = argparse.ArgumentParser(description='View a single UVM-RAD patient in Napari (read-only)')
    parser.add_argument('--root', type=Path, default=DEFAULT_ROOT, help='Root of UVM_RAD_DATASET_FINAL')
    parser.add_argument('--patient', required=True, help='Number or ID, e.g. 63 or UVM-063')
    parser.add_argument('--no-rtstruct', action='store_true', help='Load CT/MR without contour rasterization')
    parser.add_argument('--debug', action='store_true', help='Verbose logging')
    args = parser.parse_args()
    # Default: show fatal loading errors only; use --debug for warnings and ROI details.
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.ERROR,
                        format='%(levelname)s: %(message)s')
    if not args.debug:
        logging.getLogger('OpenGL').setLevel(logging.ERROR)
    launch(args.root, normalize_patient(args.patient), no_rtstruct=args.no_rtstruct)


if __name__ == '__main__':
    main()
