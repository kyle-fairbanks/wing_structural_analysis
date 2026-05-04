# Wing Structural Analysis Tool

A Python-based structural analysis tool for evaluating composite sandwich wing structures using Classical Laminate Theory (CLT), multi-cell torsion analysis, and comprehensive failure mode assessment.

## Overview

This tool performs detailed structural analysis of composite sandwich wings with the following capabilities:

- **Sandwich beam theory** with CLT-derived properties for accurate composite behavior
- **Multi-cell torsion analysis** using Bredt-Batho theory for thin-walled sections
- **Sweep angle effects** (0° to 60°) with proper load projection
- **Comprehensive stress analysis**: bending, shear (web and core), and torsion
- **Factor of Safety (FoS) calculations** for all failure modes
- **Weight estimation** by component (skins, core, webs, rods)
- **Excel and plot outputs** for detailed result inspection

### Material System

**Default configuration uses an aerospace composite material system:**
- **Face sheets**: FCIM255 fabric (HiMax CGL4012) - ±45° biax carbon, 200 gsm, with T700S 12k fiber
- **Core**: ROHACELL 31 IG-F closed-cell foam (32 kg/m³, 36 MPa modulus)
- **Spar caps**: T700S unidirectional carbon/epoxy rods (138 GPa, 1500 kg/m³)

### Analysis Scope

- **Configuration**: Half-wing cantilever beam (root-to-tip)
- **Loading**: Symmetric aerodynamic lift at specified g-limit
- **Geometry**: Uniform chord with customizable airfoil section and web positions
- **Outputs**: Internal loads (V, M, T), stresses, FoS by station and sweep angle

---

## Key Features

✅ **Classical Laminate Theory (CLT)** - CLT-based ABD matrix formulation for laminate stiffness calculation  
✅ **Multi-cell torsion** - Bredt solver with corrected single-cell handling  
✅ **GA shear partition** - GA stiffness model with weighted web split  
✅ **Sandwich wrinkling checks** - Cube-root formula for face-sheet buckling (first-order screen)  
✅ **Transformed section analysis** - Neutral axis and moment of inertia calculation for mixed materials  
✅ **Tsai-Wu failure criterion** - Composite-specific failure assessment for skins  
✅ **Automated sweep studies** - Batch analysis across multiple sweep angles  
✅ **Professional outputs** - Excel workbook + generated PNG plots for result inspection and reporting  
✅ **Built-in validation** - Bredt solver test runs automatically on startup  

---

## Quick Start

### 1. Install Dependencies

```bash
pip install numpy pandas matplotlib scipy shapely xlsxwriter
```

**Requirements:**
- Python 3.7 or later
- numpy >= 1.18
- pandas >= 1.0
- matplotlib >= 3.1
- scipy >= 1.4 (optional, for airfoil smoothing)
- shapely >= 1.7 (required, for geometry operations)
- xlsxwriter >= 1.2 (for Excel output)

### 2. Run Basic Analysis

```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53
```

This runs the default configuration with a single sweep angle (0°). Results are written to `wing_results_final.xlsx` and plot files are generated in the current directory.

### 3. View Results

**Excel Output** (`wing_results_final.xlsx`):
- Open the workbook and review the `0_Final_Summary` sheet for overall FoS
- Check `Weight_Summary` for component mass breakdown
- Inspect `2_All_Sweeps` for detailed station-by-station results

**Plot Output**:
- `airfoil_with_webs.png` - Cross-section geometry with web locations
- `lift_distribution_full_span.png` - Spanwise load distribution verification
- `FoS_governing_vs_span.png` - Minimum FoS envelope across span

---

## Installation

### Standard Installation

```bash
# Clone or download the repository
cd /path/to/wing-analysis

# Install Python dependencies
pip install numpy pandas matplotlib scipy shapely xlsxwriter

# Verify installation
python Wing.py --help
```

### Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv venv

# Activate (Linux/Mac)
source venv/bin/activate

# Activate (Windows)
venv\Scripts\activate

# Install dependencies
pip install numpy pandas matplotlib scipy shapely xlsxwriter
```

### Dependencies Reference

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | >= 1.18 | Numerical computations, arrays |
| pandas | >= 1.0 | Data tables, Excel I/O |
| matplotlib | >= 3.1 | Plotting and visualization |
| scipy | >= 1.4 | Spline interpolation (optional) |
| shapely | >= 1.7 | Geometry operations (required) |
| xlsxwriter | >= 1.2 | Excel file generation |

---

## Basic Usage Examples

### Example 1: Single Sweep Angle

```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53 --sweeps 0
```

Analyzes the wing at 0° sweep with default materials and flight conditions.

### Example 2: Sweep Study (0° to 60° in 5° increments)

```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53 --sweeps "0,5,10,15,20,25,30,35,40,45,50,55,60"
```

Performs parametric sweep study across 13 sweep angles. Results include comparative plots across all angles.

### Example 3: No-Webs Configuration

```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53 --no_webs
```

Analyzes wing with webs neutralized (ultra-thin and weak). Useful for isolating skin/core behavior or exploring web-free designs.

### Example 4: Custom Airfoil

```bash
python Wing.py --chord_m 0.15 --span_ft 5.0 --airfoil my_airfoil.txt
```

Uses a custom airfoil from `my_airfoil.txt` (must be two-column x,y format, normalized or dimensional).

### Example 5: Output to Custom File

```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53 --out my_results.xlsx --no_plots
```

Saves results to `my_results.xlsx` and skips plot generation for faster batch runs.

---

## Output Files

### Excel Workbook Structure

The main output file (default: `wing_results_final.xlsx`) contains multiple sheets:

| Sheet Name | Contents |
|------------|----------|
| `0_Final_Summary` | Worst-case sweep angle and minimum FoS across all analyses |
| `0_Equations` | Theoretical equation reference (sandwich theory, Bredt, etc.) |
| `1_Inputs` | Analysis configuration (geometry, flight conditions, settings) |
| `1_Materials` | Material properties used (E, G, allowables, densities) |
| `1_Airfoil` | Airfoil section coordinates (x, y in meters) |
| `2_All_Sweeps` | Combined station-by-station results for all sweep angles |
| `2_Summaries_All` | Summary row for each sweep (min FoS, governing mode, etc.) |
| `Weight_Summary` | Component mass breakdown (skins, core, webs, rods, total) |
| `Loads_Λ##` | Per-sweep loads (lift, V, M, T) by station |
| `Stresses_Λ##` | Per-sweep stresses (σ, τ) by station |
| `FoS_Λ##` | Per-sweep factors of safety by station |
| `Summary_Λ##` | Per-sweep summary (geometry, EI, min FoS, governing mode) |

*Note: `##` represents sweep angle (e.g., `Loads_Λ00` for 0°, `Loads_Λ30` for 30°)*

### Plot Files Generated

**Geometry:**
- `airfoil_with_webs.png` - Airfoil cross-section showing web locations

**Load Verification:**
- `lift_distribution_full_span.png` - Spanwise lift distribution (full wingspan, symmetric)

**Internal Loads (per sweep angle, overlaid):**
- `bending_moment_vs_span.png` - Bending moment M(x)
- `shear_force_vs_span.png` - Shear force V(x)
- `torque_vs_span.png` - Torsional moment T(x)

**Stresses (per sweep angle, overlaid):**
- `rod_stress_top_vs_span.png`, `rod_stress_bottom_vs_span.png` - Spar cap normal stresses
- `skin_stress_top_vs_span.png`, `skin_stress_bottom_vs_span.png` - Skin normal stresses
- `web_shear_stress_vs_span.png` - Web shear stress (VQ/It)
- `skin_torsion_stress_vs_span.png` - Skin torsion stress (Bredt)
- `core_shear_stress_vs_span.png` - Core through-thickness shear

**Factors of Safety (per sweep angle, log scale, overlaid):**
- `FoS_rod_top_vs_span.png`, `FoS_rod_bottom_vs_span.png` - Rod bending FoS
- `FoS_skin_top_vs_span.png`, `FoS_skin_bottom_vs_span.png` - Skin bending FoS
- `FoS_web_shear_vs_span.png` - Web shear FoS
- `FoS_skin_torsion_vs_span.png` - Skin torsion FoS (Tsai-Wu + simple shear)
- `FoS_core_shear_vs_span.png` - Core shear FoS
- `FoS_governing_vs_span.png` - **Minimum FoS envelope** (most critical at each station)

**Diagnostics:**
- `diagnostic_drivers.png` - Web and torsion stress driver comparison

---

## Configuration

### Editing Analysis Parameters

The main configuration section is in `Wing.py` lines 2786-2810. Key parameters:

```python
# Geometry
SPAN_FT = 4.53                # ft (wingspan)
CHORD_M = 0.10668             # m (uniform chord)
WEB1_POSITION = 0.2           # fraction of chord (front spar at 20% chord)
WEB2_POSITION = 0.65          # fraction of chord (rear spar at 65% chord)
ROD_WIDTH_MM = 3.0            # mm (spar cap width)
ROD_HEIGHT_MM = 3.0           # mm (spar cap height)
ROD_GAP_MM = 0                # mm (clearance from inner skin into foam)

# Flight Conditions
VELOCITY_MS = 15.87           # m/s
AIR_DENSITY = 1.0573          # kg/m³
LIFT_COEFFICIENT = 1.43       # dimensionless
G_LIMIT = 5.0                 # load factor

# Analysis Settings
N_STATIONS = 101              # number of spanwise stations (resolution)
SWEEP_ANGLES = list(np.linspace(0, 60, 13))  # sweep angles to analyze
```

### Modifying Materials

Edit the `Materials` dataclass (lines 41-129) to change material properties:

```python
mat = Materials(
    # Face sheet layup (top-to-bottom, will be mirrored)
    face_stack=[(45.0, 1), (-45.0, 1)],  # one ±45° biax layer
    
    # Web layup
    web_stack=[(45.0, 1), (-45.0, 1)],   # one ±45° biax layer
    
    # Ply thickness (mm)
    ply_t_skin_mm=0.13,  # FCIM255 fabric thickness
    ply_t_web_mm=0.13,   # FCIM255 fabric thickness
    
    # Core properties
    t_core=0.003,         # m (3 mm)
    rho_core=32.0,        # kg/m³
    E_core=36e6,          # Pa
    G_core=13e6,          # Pa
    tau_core_allow=0.4e6, # Pa
    
    # Laminate elastic properties (T700S UD at 60% Vf)
    E1_skin=134e9,        # Pa (longitudinal)
    E2_skin=10e9,         # Pa (transverse)
    G12_skin=5.0e9,       # Pa (in-plane shear)
    
    # Strength allowables (T700S UD at 60% Vf)
    Xt_skin=2.86e9,       # Pa (tensile)
    Xc_skin=1.45e9,       # Pa (compressive)
    tau_skin_allow=136e6, # Pa (shear)
)
```

---

## Command-Line Reference

```bash
python Wing.py [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--airfoil` | str | `Modified_NACA_6_Series.txt` | Path to airfoil coordinate file |
| `--span_ft` | float | `4.53` | Full wingspan in feet (tip to tip) |
| `--chord_m` | float | **Required** | Chord length in meters (uniform) |
| `--sweeps` | str | `0` | Comma-separated sweep angles (e.g., `"0,15,30,45"`) |
| `--out` | str | `wing_results_final.xlsx` | Output Excel filename |
| `--no_plots` | flag | Off | Skip plot generation (faster for batch runs) |
| `--no_webs` | flag | Off | Neutralize webs (analyze without structural webs) |

### Examples with Options

```bash
# Minimal run (requires chord)
python Wing.py --chord_m 0.1

# Full sweep study with custom output
python Wing.py --chord_m 0.15 --span_ft 5.0 --sweeps "0,10,20,30,40,50,60" --out sweep_study.xlsx

# Fast batch run without plots
python Wing.py --chord_m 0.12 --no_plots

# No-webs configuration study
python Wing.py --chord_m 0.10668 --no_webs --out no_webs_results.xlsx
```

---

## Understanding Results

### Factor of Safety (FoS) Interpretation

**FoS > 3.0**: Conservative design, high margin  
**FoS 1.5 - 3.0**: Adequate for preliminary design, verify with detailed analysis  
**FoS < 1.5**: Marginal, requires detailed FEA and testing  
**FoS < 1.0**: Failure predicted at design load  

### Governing Failure Modes

The tool identifies the critical failure mode at each station. Common modes:

- **Skin top/bottom**: Bending tension or compression failure
- **Rod top/bottom**: Spar cap bending failure (typically near root)
- **Web shear**: Web shear failure from combined transverse and torsional shear (RSS combination)
- **Skin torsion**: Skin shear failure from torsional loading (Bredt)
- **Core shear**: Core through-thickness shear failure
- **Skin wrinkling**: Face-sheet wrinkling (sandwich-specific buckling mode, c=0.5 coefficient)
- **Web shear buckling**: Shear buckling of web panels (infinite-strip model)

### Typical Failure Mode Locations

- **Root region**: Rod bending (high moment), web shear (high V)
- **Mid-span**: Often balanced, multiple modes near critical
- **Tip region**: Usually non-critical (low loads)

---

## Assumptions and Limitations

### Model Assumptions

- **Half-wing symmetry**: Assumes symmetric flight (no roll, yaw)
- **Uniform chord**: No taper (rectangular planform)
- **Cantilever BC**: Fully fixed root, free tip
- **Linear elastic**: No plasticity, no damage progression
- **Static loading**: No dynamic, fatigue, or aeroelastic effects

### Known Limitations

1. ~~**No skin buckling analysis**~~ **Sandwich wrinkling now included** (first-order screen)
   - Face-sheet wrinkling checked (cube-root formula, c=0.5)
   - Coefficient range 0.4-0.6 depending on assumptions
   - Does NOT check post-buckling strength or intracell dimpling
2. **Simplified core shear**: Assumes 85% effective width (b_eff) for LE/TE ineffectiveness
3. **No transverse normal stress**: Through-thickness stress not modeled
4. **Rigid airfoil section**: No deformation of cross-section shape
5. **No manufacturing effects**: Assumes perfect fiber alignment, no voids or defects

### When to Use FEA Instead

Consider detailed finite element analysis if:
- FoS < 1.5 (marginal design)
- High compression skins (buckling concern)
- Complex load cases (combined bending, torsion, axial)
- Detailed local stress at cutouts or attachments needed
- Aeroelastic flutter or divergence analysis required

---

## Troubleshooting

### Common Errors

**Error: "Chord (inp.chord_m) must be provided and > 0"**  
→ Solution: Always provide `--chord_m` argument with positive value

**Error: "web2_p must be greater than web1_p"**  
→ Solution: Ensure rear web position > front web position (e.g., web1=0.2, web2=0.65)

**Error: "Airfoil file not found"**  
→ Solution: Check that airfoil file exists in current directory or provide full path

**Warning: "Web height invalid: h_web=0.0mm"**  
→ Solution: Skins may be too thick for airfoil thickness; reduce ply count or use thicker airfoil

### Performance Tips

- Use `--no_plots` for faster batch runs (Excel output only)
- Reduce `N_STATIONS` (in code) from 101 to 51 for coarser/faster analysis
- Run single sweep first to verify configuration before multi-sweep study

---

## Further Documentation

- **[USER_MANUAL.md](USER_MANUAL.md)** - Comprehensive usage guide with detailed examples
- **[TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)** - Theory, equations, and material properties
- **[VALIDATION.md](VALIDATION.md)** - Test cases, verification procedures, and benchmarks

---

## Support and Contributing

### Getting Help

1. Review the [USER_MANUAL.md](USER_MANUAL.md) for detailed usage instructions
2. Check [VALIDATION.md](VALIDATION.md) for verification procedures
3. Inspect terminal output for diagnostic messages and validation results

### Reporting Issues

When reporting issues, please include:
- Full command line used
- Python version (`python --version`)
- Package versions (`pip list`)
- Terminal output (including error messages)
- Input file (airfoil coordinates if custom)

### Contributing

Contributions welcome! Areas for enhancement:
- Tapered wing geometry (variable chord)
- Non-rectangular planform (elliptical, trapezoidal)
- Dynamic load cases (gust, landing impact)
- Skin buckling analysis
- Multi-material optimization

---

## License and Citation

**For academic use**, please cite:
> Wing Structural Analysis Tool (2025). Composite sandwich wing analysis using CLT and multi-cell torsion theory. University of Arizona Senior Design Capstone Team 26012.

**Disclaimer**: This tool is for preliminary design and research purposes. All designs must be validated with detailed FEA and physical testing before manufacturing or flight.

---

## Version History

**Current Version**: 1.0  
- CLT-enabled composite analysis
- Multi-cell Bredt torsion solver
- GA-based shear partition (webs vs. core)
- Transformed section analysis with correct neutral axis
- Tsai-Wu composite failure criterion
- Comprehensive Excel and plot outputs
- Built-in validation tests

---

*For detailed usage instructions, see [USER_MANUAL.md](USER_MANUAL.md)*


