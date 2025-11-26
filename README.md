**Project**
- **Name**: `kerfer` — small SVG path utilities for offsetting paths to account for laser cutter kerf.
- **Purpose**: Analyze SVG path winding, open/close subpaths, and perform simple linear offsets on path outlines using the excellent `svgelements` library (https://github.com/meerk40t/svgelements).

**Requirements**
- **Python**: 3.8+ recommended.
- **Dependencies**: See `requirements.txt` (uses `svgelements`).

**Install**
- **Create virtual environment and install**:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

**Command Line Interface (CLI)**
- **Overview**: The project includes a small CLI in `kerfer.py` so you can run operations without editing the source. The CLI accepts an input file path, operation flags to open, offset, and close subpaths, and an output file path.
Normal sequence of operation is open, offset, and then close.
- **Basic command**:

```powershell
.\.venv\Scripts\python.exe .\kerfer.py -i "..\svgs\box1.svg" --open --offset 0.15 --close --output "box1_dilated_0.15.svg"
```

- **Flags**:
	- `-i, --input`: Input SVG file path. If none provided, then a small test file is generated.
	- `--open`: Open closed subpaths (replace `Close` segments with `Line`).
	- `--offset`: Floating-point offset distance (same units as the SVG). Example: `--offset 0.15`.
	- `--close`: Close open subpaths (replace final `Line` of each segment with a `Close` when subpath endpoints match).
	- `-o, --output`: Write modified SVG to this path. If omitted the script prints a short summary to stdout.
	- `-v, --verbose`: Enable debug logging output.

- **Operation order**: When multiple operations are given they are applied in this order: `--open` -> `--offset` -> `--close`.

- **Examples**:

```powershell
# Print a summary of the SVG (no file written)
.\.venv\Scripts\python.exe .\kerfer.py -i "D:\Projects\kerfer\box_asploded.svg"

# Open, offset by 0.15 units, and write to output.svg
.\.venv\Scripts\python.exe .\kerfer.py -i "input.svg" -o "output.svg" --open --offset 0.15

# Close only, verbose logging enabled
.\.venv\Scripts\python.exe .\kerfer.py -i "input.svg" -o "closed.svg" --close -v
```

**Development notes**
- **Code structure**: Top-level functions provide the main functionality: `calculate_is_clockwise`, `open_svg`/`close_svg`, and `offset_svg`.
- **Limitations**: Curved segments (`CubicBezier`, `QuadraticBezier`, `Arc`) are currently treated only by their endpoints for area/offset calculations — full-curve handling requires sampling or analytic integration.
- **Testing**: There are no automated tests included. To experiment, create a small test `SVG` and use `svgelements.SVG.parse()` as shown in `kerfer.py`.

**Contributing**
- **Suggestions**: If you want better curve handling, consider sampling bezier/arc segments or using a geometry library to compute exact curve contributions to signed area and offsetting.

**License**
- This project is released under the MIT License — see the `LICENSE` file included in the repository.

**Contact / Attribution**
- This is a small personal utility; adapt freely for local use.

