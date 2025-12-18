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
	- `-o, --output`: Write modified SVG to this path. If omitted, the script generates a unique filename that includes the dilation amount.
	- `-b, --break`: Break apart subpaths into separate paths
	- `-n, --nest`: Nest subpaths into parent paths
	- `-l, --line`: Open closed subpaths (replace each Close with a *Line*)
	- `-z, --zero_cull`: Removes *zero-length* segments from paths
	- `-s, --simplify`: Remove unnecessary points from paths to *simplify* them
	- `-d, --dilate`: Perpendicular offset *dilation* distance (in same units as SVG)
	- `-c, --close`: Close open subpaths (replace final `Line` of each subpath with a `Close` when subpath endpoints match).
	- `-r, --rebreak`: Rebreak subpaths into separate paths for laser cutting
	- `-a, --all`: Default if no other processing specified except dilation: Do *all* the steps - break, nest, line, zero-cull, simplify, dilate, close, rebreak
	- `-v, --verbose`: Enable debug logging output.

- **Operation order**: When multiple operations are given they are applied in this order: break, nest, line, zero-cull, simplify, dilate, close, rebreak.

- **Examples**:

```powershell
# Open, offset by 0.15 units, and write to output.svg
.\.venv\Scripts\python.exe .\kerfer.py -i "input.svg" -a -d 0.15 -o "output.svg"

# Open, offset by 0.2 units, and write to output_d0.2mm.svg (if input has mm for units)
.\.venv\Scripts\python.exe .\kerfer.py -i "input.svg" -a -d 0.15
```

**Development notes**
- **Code structure**: Top-level functions provide the main functionality: e.g., `open_svg`, `dilate_svg`, and `close_svg`.
- **Limitations**: Curved segments (`CubicBezier`, `QuadraticBezier`, `Arc`) are currently treated only by their endpoints for area/offset calculations — full-curve handling requires sampling or analytic integration.
- **Testing**: Automated tests included. To demonstrate, run the app with no imput file specified; a small test SVG will be created and processed, with each intermediate result written to a test output directory.

**Contributing**
- **Suggestions**: For curve handling, consider sampling bezier/arc segments or using a geometry library to compute exact curve contributions to signed area and offsetting.

**License**
- This project is released under the MIT License — see the `LICENSE` file included in the repository.

**Contact / Attribution**
- This is a small personal utility; adapt freely for local use.

