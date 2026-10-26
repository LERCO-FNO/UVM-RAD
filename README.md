# UVM-RAD DICOM dataset reader

A read-only Napari viewer for CT, available MR series, and RTSTRUCT contours in
the UVM-RAD DICOM dataset. The script does not modify the dataset or perform
MRI-to-CT registration. RTSTRUCT contours are rasterized on the CT grid for
display.

The dataset itself is not included in this repository. You need access to a
local copy of it before running the viewer.

## Requirements

- Windows, macOS, or Linux
- Python 3.10 or newer (Python 3.12 is recommended)
- [uv](https://docs.astral.sh/uv/)
- A local copy of the UVM-RAD dataset

## Install uv

Follow the [official uv installation guide](https://docs.astral.sh/uv/getting-started/installation/)
if uv is not already installed.

On Windows, install with WinGet:

```powershell
winget install --id=astral-sh.uv -e
```

On macOS or Linux, use the installer from Astral:

```sh
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Open a new terminal after installation and check that uv is available:

```sh
uv --version
```

## Set up the project

In a terminal, change to the directory where you cloned this repository:

```powershell
cd D:\path\to\UVM-RAD
```

Then synchronize the project environment:

```sh
uv sync
```

`uv sync` creates the local `.venv` environment and installs the dependencies
declared in `pyproject.toml`, using the exact versions in `uv.lock`. The
`requirements.txt` file lists the same direct dependencies for tools that use
the standard pip requirements format.

If uv cannot find a compatible Python installation, install one and retry:

```sh
uv python install 3.12
uv sync --python 3.12
```

## Run the viewer

The dataset root directory should contain patient directories such as
`UVM-063`, with the CT, MR, and optional RTSTRUCT folders inside each patient
directory. By default, the viewer uses the values in `config.py`:

- `DATASET_ROOT` is the dataset directory.
- `PATIENT_ID` is the patient to open (currently `1`).

### Option 1: Run in VS Code using `config.py`

Edit `DATASET_ROOT` and `PATIENT_ID` in `config.py` as needed. In VS Code,
select the project's `.venv` Python interpreter, open `uvm_dataset_reader.py`,
and choose **Run Python File**. The configured patient and dataset path are
used when no command-line options are supplied.

### Option 2: Run from the command line

Pass `--root` and `--patient` to override the defaults from `config.py`:

Windows (PowerShell):

```powershell
uv run python .\uvm_dataset_reader.py --patient 28 --root "D:\FNO_image_segmentation\TCIA uveal - anonym\UVM_RAD_DATASET_FINAL"
```

macOS or Linux:

```sh
uv run python uvm_dataset_reader.py --patient 63 --root "/path/to/UVM_RAD_DATASET_FINAL"
```

The patient can be provided as a number (`63`) or an ID (`UVM-063`). If
`--patient` or `--root` is omitted, its value from `config.py` is used. To see
all available options:

```sh
uv run python uvm_dataset_reader.py --help
```

Useful options:

- `--no-rtstruct` loads CT and MR images without searching for or rasterizing
  RTSTRUCT contours.
- `--debug` enables warnings and detailed logging.

The viewer displays images in their DICOM physical coordinates. It does not
register MR images to CT, so images with different frames of reference should
not be interpreted as aligned.
