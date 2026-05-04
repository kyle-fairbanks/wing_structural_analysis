# Wing Structural Analysis - User Manual

Comprehensive operational guide for composite sandwich wing structural analysis.

**Version 1.0**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Input Configuration](#2-input-configuration)
3. [Running Analyses](#3-running-analyses)
4. [Interpreting Results](#4-interpreting-results)
5. [Customization Examples](#5-customization-examples)
6. [Troubleshooting](#6-troubleshooting)
7. [Advanced Topics](#7-advanced-topics)

---

## 1. Introduction

### 1.1 Purpose

This manual provides detailed instructions for using the Wing Structural Analysis tool to design and analyze composite sandwich wings. The tool is intended for:

- Preliminary wing structural design
- Parametric studies (sweep angle, layup, geometry)
- Design iteration and optimization
- Research applications

### 1.2 Prerequisites

**Required knowledge:**
- Basic composite mechanics (fiber/matrix, laminate orientation)
- Structural analysis fundamentals (bending, shear, torsion)
- Wing aerodynamics (lift distribution, loading)

**Software requirements:**
- Python 3.7 or later
- See [README.md](README.md) for package dependencies

### 1.3 Analysis Workflow

```
Define Geometry → Configure Materials → Set Flight Conditions → Run Analysis → Review Results → Iterate
```

Typical design iteration takes 5-30 seconds per sweep angle depending on resolution.

---

## 2. Input Configuration

### 2.1 Geometry Parameters

#### 2.1.1 Wing Planform

**Span (`span_ft` or `SPAN_FT` in code)**
- **Units**: feet, full wingspan. The analysis is performed on one semi-span from root to tip
- **Range**: 2 - 20 ft typical for small UAV/RC aircraft
- **Default**: 4.53 ft
- **Command line**: `--span_ft 4.53`

**Chord (`chord_m` or `CHORD_M` in code)**
- **Units**: meters (uniform, no taper)
- **Range**: 0.05 - 0.50 m typical
- **Default**: Must be provided (no default)
- **Command line**: `--chord_m 0.10668`
- **⚠️ REQUIRED**: Analysis will fail if not specified

#### 2.1.2 Web Positions

Webs are internal shear-carrying structures at fixed chordwise locations.

**Front Web (`web1_p` or `WEB1_POSITION`)**
- **Units**: fraction of chord (0 = LE, 1 = TE)
- **Range**: 0.10 - 0.40 typical
- **Default**: 0.2 (20% chord)
- **Location**: Edit in code line 2795

**Rear Web (`web2_p` or `WEB2_POSITION`)**
- **Units**: fraction of chord
- **Range**: 0.50 - 0.80 typical
- **Default**: 0.65 (65% chord)
- **Location**: Edit in code line 2796
- **⚠️ CONSTRAINT**: Must be > web1_p

**Typical configurations:**
```
Light loading:  web1=0.25, web2=0.70 (wider spacing)
Heavy loading:  web1=0.15, web2=0.60 (narrower spacing)
Single spar:    web1=0.30, web2=0.31 (minimal separation, effectively one D-cell)
```

#### 2.1.3 Spar Caps (Rods)

Rods are spanwise unidirectional reinforcements that carry bending loads.

**Rod Width (`rod_width_m` or `ROD_WIDTH_MM`)**
- **Units**: meters (or mm in config constant)
- **Range**: 2 - 10 mm typical
- **Default**: 3.0 mm
- **Location**: Edit in code line 2797

**Rod Height (`rod_height_m` or `ROD_HEIGHT_MM`)**
- **Units**: meters (or mm in config constant)
- **Range**: 2 - 10 mm typical
- **Default**: 3.0 mm
- **Location**: Edit in code line 2798
- **Note**: Must fit inside core thickness

**Rod Gap (`rod_gap_mm` or `ROD_GAP_MM`)**
- **Units**: mm (clearance from inner skin surface into foam)
- **Range**: 0.2 - 2.0 mm typical
- **Default**: 0.5 mm
- **Location**: Edit in code line 2799
- **Purpose**: Manufacturing clearance, prevents rod contact with skin

**Rod Count and Distribution:**
The tool assumes **4 total rods**: 2 at top, 2 at bottom, distributed chordwise.

Default positions (fraction of chord):
```python
rod_positions_xc = [0.2, 0.2, 0.65, 0.65]  # 2 forward @ 20%, 2 aft @ 65%
```

To change positions, edit line 1156 in `Wing.py`.

#### 2.1.4 Airfoil Section

**Airfoil File (`airfoil_file`)**
- **Format**: Two-column text file (x, y coordinates)
- **Coordinate system**: Normalized (0-1) or dimensional (will be normalized)
- **Point order**: Leading edge first, wrap around airfoil back to LE
- **Default**: `Modified_NACA_6_Series.txt`
- **Command line**: `--airfoil my_airfoil.txt`

**File format example:**
```
1.000000  0.000000
0.995000  0.001234
0.990000  0.002456
...
0.000000  0.000500  (leading edge)
...
0.010000 -0.002345
...
1.000000  0.000000
```

**Requirements:**
- At least 20 points for smooth geometry
- Closed loop (first point = last point, or will be closed automatically)
- No duplicate consecutive points
- LE should be near x = 0, TE near x = 1

---

### 2.2 Material Configuration

Materials are defined in the `Materials` dataclass (lines 41-129 in `Wing.py`).

#### 2.2.1 Face Sheet (Skin) Layup

**Face Stack (`face_stack`)**

Defines the ply sequence for ONE face (top or bottom). The stack is applied symmetrically to both faces.

**Format**: List of (angle_deg, n_plies) tuples
```python
face_stack = [(45.0, 1), (-45.0, 1)]  # One ±45° biaxial layer
```

**Interpretation:**
- `(45.0, 1)`: One ply at +45°
- `(-45.0, 1)`: One ply at -45°
- Total face thickness: 2 plies × 0.13 mm/ply = 0.26 mm

**Common layup patterns:**

| Layup | Description | face_stack |
|-------|-------------|------------|
| Single ±45° biax | Light, shear-dominated | `[(45, 1), (-45, 1)]` |
| Double ±45° biax | Moderate thickness | `[(45, 2), (-45, 2)]` |
| 0°/±45° hybrid | Bending-optimized | `[(0, 1), (45, 1), (-45, 1)]` |
| Quasi-isotropic | Balanced properties | `[(0, 1), (45, 1), (-45, 1), (90, 1)]` |

#### 2.2.2 Web Layup

**Web Stack (`web_stack`)**

Defines the ply sequence for internal webs.

**Format**: Same as face_stack
```python
web_stack = [(45.0, 1), (-45.0, 1)]  # One ±45° biax layer (optimal for shear)
```

**Design guidance:**
- ±45° orientation is optimal for shear webs
- 0° plies contribute little to web shear capacity
- Typical: use same layup as skins for manufacturing simplicity

#### 2.2.3 Ply Thickness

**Skin Ply Thickness (`ply_t_skin_mm`)**
- **Units**: mm per ply
- **Default**: 0.13 mm (FCIM255 fabric)
- **Range**: 0.10 - 0.25 mm typical for carbon fabrics

**Web Ply Thickness (`ply_t_web_mm`)**
- **Units**: mm per ply
- **Default**: 0.13 mm (same fabric as skins)
- **Note**: Can differ from skins if using different fabric

#### 2.2.4 Core Properties

**Core Thickness (`t_core`)**
- **Units**: meters
- **Default**: 0.003 m (3.0 mm)
- **Range**: 2 - 10 mm typical for small wings
- **Trade-off**: Thicker = higher bending stiffness, heavier

**Core Density (`rho_core`)**
- **Units**: kg/m³
- **Default**: 32.0 (ROHACELL 31 IG-F)
- **Range**: 25 - 200 typical for foam cores

**Core Elastic Modulus (`E_core`)**
- **Units**: Pa
- **Default**: 36e6 Pa (36 MPa, ROHACELL 31 IG-F)
- **Note**: Tension/compression modulus

**Core Shear Modulus (`G_core`)**
- **Units**: Pa
- **Default**: 13e6 Pa (13 MPa, ROHACELL 31 IG-F)
- **Critical**: Drives core shear capacity

**Core Shear Allowable (`tau_core_allow`)**
- **Units**: Pa
- **Default**: 0.4e6 Pa (0.4 MPa, ROHACELL 31 IG-F)
- **Design driver**: Often critical in highly loaded regions

#### 2.2.5 Laminate Elastic Properties

These properties represent the **cured laminate** behavior (fiber + resin system).

**Skin Properties:**
- `E1_skin`: 134e9 Pa (longitudinal modulus, T700S UD at 60% Vf)
- `E2_skin`: 10e9 Pa (transverse modulus)
- `G12_skin`: 5.0e9 Pa (in-plane shear modulus)
- `nu12_skin`: 0.30 (Poisson's ratio)

**Web Properties:**
- Same as skin by default (same fabric)
- Edit separately if using different material

#### 2.2.6 Strength Allowables

**Skin Allowables:**
- `Xt_skin`: 2.86e9 Pa (2860 MPa, tensile strength T700S UD at 60% Vf)
- `Xc_skin`: 1.45e9 Pa (1450 MPa, compressive strength)
- `tau_skin_allow`: 136e6 Pa (136 MPa, in-plane shear)

**Web Allowables:**
- `tau_web_allow`: 136e6 Pa (same as skin, in-plane shear)

**Rod Allowables:**
- `Xt_rod`: 1.72e9 Pa (1720 MPa, tensile)
- `Xc_rod`: 1.83e9 Pa (1830 MPa, compressive/flexural)

---

### 2.3 Flight Conditions

#### 2.3.1 Aerodynamic Inputs

**Airspeed (`V_ms` or `VELOCITY_MS`)**
- **Units**: m/s
- **Default**: 15.87 m/s (~35 mph, typical small RC)
- **Location**: Edit in code line 2802

**Air Density (`rho` or `AIR_DENSITY`)**
- **Units**: kg/m³
- **Default**: 1.0573 (sea level, 15°C per ISA)
- **Location**: Edit in code line 2803

**Lift Coefficient (`cl` or `LIFT_COEFFICIENT`)**
- **Units**: dimensionless
- **Default**: 1.43 (near max lift for many airfoils)
- **Location**: Edit in code line 2804
- **Note**: Section CL, not 3D wing CL

**Dynamic Pressure (`q_pa`):**
Computed automatically from V and ρ:
```
q = 0.5 × ρ × V²
```
If provided explicitly, overrides V/ρ.

#### 2.3.2 Load Factor

**G-Limit (`g_limit` or `G_LIMIT`)**
- **Units**: dimensionless (multiples of 1g)
- **Default**: 5.0 (typical maneuvering load)
- **Location**: Edit in code line 2805
- **Range**: 2.5 (light loads) to 12 (aerobatic)

**Total Wing Lift:**
```
L_total = AUW × g_limit  (at design condition)
```

The lift distribution is normalized so that integrating over the half-wing gives exactly 0.5 × L_total.

---

### 2.4 Analysis Settings

#### 2.4.1 Spanwise Stations

**Number of Stations (`n_stations` or `N_STATIONS`)**
- **Units**: integer count
- **Default**: 101
- **Location**: Edit in code line 2808
- **Range**: 21 (coarse, fast) to 201 (fine, slow)
- **Trade-off**: More stations = better resolution, slower analysis

#### 2.4.2 Sweep Angles

**Sweep Angles (`sweeps` or `SWEEP_ANGLES`)**
- **Units**: degrees
- **Default**: 0° to 60° in 5° increments (13 angles)
- **Command line**: `--sweeps "0,15,30,45,60"`
- **Location**: Edit in code line 2809

**Single sweep example:**
```bash
python Wing.py --chord_m 0.1 --sweeps 0
```

**Multi-sweep example:**
```bash
python Wing.py --chord_m 0.1 --sweeps "0,10,20,30,40,50"
```

#### 2.4.3 Output Settings

**Excel Output File (`out_excel` or `--out`)**
- **Default**: `wing_results_final.xlsx`
- **Command line**: `--out my_results.xlsx`

**Plot Generation (`make_plots` or `--no_plots`)**
- **Default**: Plots enabled
- **Command line**: `--no_plots` to disable

---

## 3. Running Analyses

### 3.1 Command-Line Interface

#### 3.1.1 Basic Syntax

```bash
python Wing.py [OPTIONS]
```

#### 3.1.2 Complete Option Reference

| Option | Type | Required | Default | Description |
|--------|------|----------|---------|-------------|
| `--airfoil` | string | No | `Modified_NACA_6_Series.txt` | Airfoil coordinate file |
| `--span_ft` | float | No | `4.53` | Wingspan (feet) |
| `--chord_m` | float | **Yes** | None | Chord length (meters) |
| `--sweeps` | string | No | `"0"` | Comma-separated sweep angles |
| `--out` | string | No | `wing_results_final.xlsx` | Output Excel filename |
| `--no_plots` | flag | No | Off | Skip plot generation |
| `--no_webs` | flag | No | Off | Neutralize webs (no structural contribution) |

#### 3.1.3 Running Single Sweep

**Minimum required:**
```bash
python Wing.py --chord_m 0.10668
```

**With all options explicit:**
```bash
python Wing.py --airfoil Modified_NACA_6_Series.txt --span_ft 4.53 --chord_m 0.10668 --sweeps 0 --out results.xlsx
```

#### 3.1.4 Running Multi-Sweep Study

**Example: 0° to 60° in 15° steps**
```bash
python Wing.py --chord_m 0.10668 --sweeps "0,15,30,45,60"
```

**Example: Dense sweep study (5° increments)**
```bash
python Wing.py --chord_m 0.10668 --sweeps "0,5,10,15,20,25,30,35,40,45,50,55,60"
```

Analysis time: ~5-10 seconds per sweep angle @ 101 stations.

---

### 3.2 Configuration Section in Code

For parameters not exposed via command line, edit the configuration section (lines 2786-2810):

```python
# --- Geometry ---
SPAN_FT = 4.53                # ft (default if not overridden by CLI)
CHORD_M = 0.10668             # m
WEB1_POSITION = 0.2           # fraction of chord (x/c)
WEB2_POSITION = 0.65          # fraction of chord (x/c)
ROD_WIDTH_MM = 3.0            # mm
ROD_HEIGHT_MM = 3.0           # mm
ROD_GAP_MM = 0.5              # mm

# --- Flight Conditions ---
VELOCITY_MS = 15.87           # m/s
AIR_DENSITY = 1.0573          # kg/m³
LIFT_COEFFICIENT = 1.43       # dimensionless
G_LIMIT = 5.0                 # load factor

# --- Analysis Settings ---
N_STATIONS = 101              # number of spanwise stations
SWEEP_ANGLES = list(np.linspace(0, 60, 13))  # 0° to 60° in 5° increments
```

**Workflow:**
1. Edit values in configuration section
2. Save `Wing.py`
3. Run: `python Wing.py --chord_m 0.10668`

---

### 3.3 No-Webs Mode

Use `--no_webs` flag to analyze wing without structural contribution from webs.

**Purpose:**
- Isolate skin/core shear behavior
- Study web-free designs (monocoque)
- Verify web effectiveness (compare with vs. without)

**Implementation:**
```bash
python Wing.py --chord_m 0.10668 --no_webs
```

**Effect:**
- Web thickness reduced to 1 μm (negligible)
- Web modulus reduced to 1 MPa (ultra-weak)
- Web FoS set to 1e10 (non-limiting, excluded from governing mode)
- Web results hidden from Excel output

**Typical outcome:** Core shear becomes critical (ROHACELL limit ~0.4 MPa).

---

### 3.4 Batch Processing

For parametric studies, use shell scripting:

**Example: Chord study**
```bash
for chord in 0.08 0.10 0.12 0.14 0.16; do
    python Wing.py --chord_m $chord --out results_c${chord}.xlsx --no_plots
done
```

**Example: Span study**
```bash
for span in 3.0 4.0 5.0 6.0; do
    python Wing.py --chord_m 0.10668 --span_ft $span --out results_span${span}.xlsx --no_plots
done
```

---

## 4. Interpreting Results

### 4.1 Excel Output Structure

The output workbook contains multiple sheets organized hierarchically.

#### 4.1.1 Summary Sheets

**Sheet: `0_Final_Summary`**

Single-row summary of worst-case results across all sweeps.

| Column | Description |
|--------|-------------|
| `Worst Sweep (deg)` | Sweep angle with lowest overall FoS |
| `Min FoS @Design G` | Minimum factor of safety at g_limit across all sweeps and stations |

**Interpretation:**
- If Min FoS > 3.0: Very conservative design
- If Min FoS 1.5 - 3.0: Reasonable preliminary design
- If Min FoS < 1.5: Marginal, needs refinement or detailed FEA

**Sheet: `0_Equations`**

Reference table of governing equations used in analysis (sandwich theory, Bredt, etc.).

#### 4.1.2 Input Documentation Sheets

**Sheet: `1_Inputs`**

Single row with all input parameters (geometry, flight conditions, settings).

**Sheet: `1_Materials`**

Single row with all material properties (E, G, ρ, allowables).

**Sheet: `1_Airfoil`**

Two columns (x, y) with scaled airfoil coordinates in meters.

#### 4.1.3 Weight Summary

**Sheet: `Weight_Summary`**

Component-by-component mass breakdown for the **full wing** (both halves).

| Parameter | Typical Value | Notes |
|-----------|---------------|-------|
| `Planform area (full) [m^2]` | 0.30 | Full wing area |
| `Skins mass [kg]` | 0.12 | Top + bottom faces, both halves |
| `Core mass [kg]` | 0.029 | Foam core |
| `Webs mass [kg]` | 0.008 | Internal webs |
| `Rods mass [kg]` | 0.041 | Spar caps (4 rods spanwise) |
| `TOTAL mass [kg]` | 0.20 | Sum of all components |
| `TOTAL weight [N]` | 1.96 | Mass × 9.80665 |
| `TOTAL weight [lb]` | 0.44 | Weight in pounds |
| `Half-wing mass [kg]` | 0.10 | One half only |

**Typical mass fractions:**
- Skins: 50-65%
- Rods: 20-30%
- Core: 10-20%
- Webs: 3-7%

#### 4.1.4 Combined Multi-Sweep Results

**Sheet: `2_All_Sweeps`**

Station-by-station results for **all sweep angles combined** into one table.

**Key columns:**

| Column | Units | Description |
|--------|-------|-------------|
| `sweep_deg` | deg | Sweep angle for this row |
| `i` | - | Station index (0 = root, n-1 = tip) |
| `x (m)` | m | Spanwise position from root |
| `lift_per_unit_span (N/m)` | N/m | Distributed lift intensity |
| `V (N)` | N | Shear force (flapwise) |
| `M (N·m)` | N·m | Bending moment (flapwise) |
| `T (N·m)` | N·m | Torsional moment |
| `sigma_skin_top (Pa)` | Pa | Top skin normal stress |
| `sigma_skin_bot (Pa)` | Pa | Bottom skin normal stress |
| `sigma_rod_top (Pa)` | Pa | Top rod normal stress |
| `sigma_rod_bot (Pa)` | Pa | Bottom rod normal stress |
| `tau_web_shear (Pa)` | Pa | Web shear stress (VQ/It) |
| `tau_skin_tor (Pa)` | Pa | Skin torsion stress (Bredt) |
| `tau_core (Pa)` | Pa | Core through-thickness shear |
| `FoS_skin_top` | - | Skin top factor of safety |
| `FoS_skin_bot` | - | Skin bottom factor of safety |
| `FoS_rod_top` | - | Rod top factor of safety |
| `FoS_rod_bot` | - | Rod bottom factor of safety |
| `FoS_web_shear` | - | Web shear factor of safety |
| `FoS_skin_tor` | - | Skin torsion factor of safety |
| `FoS_core_shear` | - | Core shear factor of safety |

**Usage:**
- Filter by `sweep_deg` to view single sweep
- Sort by FoS columns to find critical stations
- Plot columns vs. `x (m)` for spanwise trends

**Sheet: `2_Summaries_All`**

One row per sweep angle with summary metrics.

| Column | Description |
|--------|-------------|
| `Sweep (deg)` | Sweep angle |
| `Chord (m)` | Chord length |
| `Span (m)` | Full wing span |
| `Web1 height (m)` | Front web height (after subtracting skins) |
| `Web2 height (m)` | Rear web height (after subtracting skins) |
| `EI_total (N·m^2)` | Total bending stiffness |
| `Min FoS @g_limit` | Minimum FoS at this sweep |
| `Max G capability` | Maximum g-load before failure (= Min FoS × g_limit) |
| `Governing failure mode` | Critical failure mode (e.g., "Rod top", "Web shear") |

**Interpretation:**
- `Min FoS @g_limit` decreasing with sweep → sweep makes structure worse
- `Governing failure mode` changing with sweep → different critical mechanisms

#### 4.1.5 Per-Sweep Detailed Sheets

For each sweep angle, four sheets are generated:

**Sheet: `Loads_Λ##`** (e.g., `Loads_Λ00` for 0°, `Loads_Λ30` for 30°)

Columns: `i`, `x (m)`, `lift_per_unit_span (N/m)`, `V (N)`, `M (N·m)`, `T (N·m)`

**Sheet: `Stresses_Λ##`**

Columns: `i`, `x (m)`, `sigma_*`, `tau_*` (all stress components)

**Sheet: `FoS_Λ##`**

Columns: `i`, `x (m)`, `FoS_*` (all FoS values)

**Sheet: `Summary_Λ##`**

Single-row summary for this sweep (geometry, EI, min FoS, governing mode)

---

### 4.2 Plot Outputs

All plots are saved as PNG files in the current directory.

#### 4.2.1 Geometry Plots

**File: `airfoil_with_webs.png`**

Shows airfoil cross-section with vertical dashed lines at web positions.

**Purpose:** Verify web locations relative to airfoil geometry.

#### 4.2.2 Load Verification Plots

**File: `lift_distribution_full_span.png`**

Spanwise lift distribution w(y) across full wingspan (−b/2 to +b/2).

**Purpose:** Verify load shape and symmetry.

**Expected shape:**
- Broad plateau inboard (tuned elliptic with shoulders)
- Smooth decay to zero at tip
- Symmetric about y = 0

**Check:** Integrate area under curve should equal L_total = AUW × g_limit.

#### 4.2.3 Internal Load Plots

All load plots show multiple curves (one per sweep angle) on the same axes.

**File: `bending_moment_vs_span.png`**

- **X-axis:** Spanwise position (m, root to tip)
- **Y-axis:** Bending moment M (N·m)
- **Trend:** Maximum at root, decreases to zero at tip

**File: `shear_force_vs_span.png`**

- **Y-axis:** Shear force V (N)
- **Trend:** Maximum near root, decreases toward tip

**File: `torque_vs_span.png`**

- **Y-axis:** Torsional moment T (N·m)
- **Trend:** Accumulates from tip to root if lift is offset from structural axis

**Interpretation:**
- Higher sweep → lower M, V (cos Λ projection effect)
- T should be continuous and smooth

#### 4.2.4 Stress Plots

**Normal Stresses (Bending):**
- `rod_stress_top_vs_span.png`: Spar cap top (tension for positive M)
- `rod_stress_bottom_vs_span.png`: Spar cap bottom (compression)
- `skin_stress_top_vs_span.png`: Skin top face
- `skin_stress_bottom_vs_span.png`: Skin bottom face

**Shear Stresses:**
- `web_shear_stress_vs_span.png`: Web shear from transverse shear force (VQ/It)
- `skin_torsion_stress_vs_span.png`: Skin shear from torsion (Bredt)
- `core_shear_stress_vs_span.png`: Core through-thickness shear

**Typical trends:**
- Normal stresses: Maximum near root (high M)
- Web shear: Maximum near root (high V)
- Skin torsion: Depends on moment arm and T distribution

#### 4.2.5 Factor of Safety Plots

**All FoS plots use log scale** to show wide range of values clearly.

**Files:**
- `FoS_rod_top_vs_span.png`, `FoS_rod_bottom_vs_span.png`
- `FoS_skin_top_vs_span.png`, `FoS_skin_bottom_vs_span.png`
- `FoS_web_shear_vs_span.png`
- `FoS_skin_torsion_vs_span.png`
- `FoS_core_shear_vs_span.png`

**Most important: `FoS_governing_vs_span.png`**

Shows the **minimum FoS envelope** (worst of all modes at each station).

**Interpretation:**
- Flat regions → margin, not critical
- Dips → critical locations
- Values < 3 → needs attention
- Values < 1.5 → marginal design
- Values < 1.0 → predicted failure

---

### 4.3 Understanding Factor of Safety Results

#### 4.3.1 FoS Definition

```
FoS = Allowable Stress / Applied Stress
```

For normal stresses (tension/compression), sign matters:
```
FoS_tension = Xt / σ      (if σ > 0)
FoS_compression = Xc / |σ|  (if σ < 0)
```

For shear stresses (always positive magnitude):
```
FoS_shear = τ_allow / |τ|
```

#### 4.3.2 Typical FoS Values by Mode

| Failure Mode | Typical FoS Range | Governing Region |
|--------------|-------------------|------------------|
| Rod top | 3 - 10 | Root (high moment, tension) |
| Rod bottom | 2 - 8 | Root (high moment, compression) |
| Skin top | 5 - 20 | Root to mid-span |
| Skin bottom | 5 - 20 | Root to mid-span |
| Web shear | 2 - 15 | Root (high shear force) |
| Skin torsion | 10 - 50 | Mid-span to tip |
| Core shear | 1.5 - 5 | **Often critical**, root region |
| Skin wrinkling top/bot | 2 - 10 | Compression regions, depends on core stiffness (±20% uncertainty in c coefficient) |
| Web shear buckling | 8 - 40 | High shear regions, infinite-strip model (conservative for continuous beam) |

#### 4.3.3 Design Criteria

**Conservative design:** FoS ≥ 3.0 everywhere

**Typical design targets:**
- Ultimate load (g_limit): FoS ≥ 1.5
- Limit load (2/3 × g_limit): FoS ≥ 2.25
- Proof load (1/2 × g_limit): FoS ≥ 3.0

**Buckling hierarchy for sandwich panels:**
- Face-sheet wrinkling typically governs for foam cores (E_core < 200 MPa)
- Material strength may govern for very stiff cores (honeycomb, metal)
- Wrinkling stress scales as ∛(E_face · E_core · G_core)
- First-order screen: expect ±20% uncertainty from coefficient choice

**Stiffness-weighted web split:**
- Shear force split between webs by their shear stiffness GA = G·t·h
- Handles asymmetric webs robustly (different heights or materials)

**Failure hierarchy preference:**
1. Core shear (predictable, benign failure mode)
2. Web shear (local, repairable)
3. Skin/rod bending (structural failure, avoid)

#### 4.3.4 Governing Failure Mode Analysis

The `Governing failure mode` column in summary sheets identifies the critical mode.

**Common scenarios:**

**"Rod bottom" or "Rod top":**
- Bending-critical design
- Typical for slender, lightly loaded wings
- **Action:** Increase rod size, add more rods, or stiffen skins

**"Web shear":**
- Shear-critical in heavily loaded, stubby wings
- **Action:** Thicken webs, add plies, or add more webs

**"Core shear":**
- Very common, especially in no-webs mode
- **Action:** Use denser core (ROHACELL 71, 110), thicken core, or add webs

**"Skin torsion":**
- Torsion-critical, often at high sweep or with large moment arm
- **Action:** Optimize web spacing, reduce lever arm, or thicken skins

#### 4.3.5 Combined Failure Modes

**Web total shear stress:**
- Combines transverse shear (from GA partition) and torsion (from Bredt) using RSS:
  ```
  τ_total = √(τ_trans² + τ_tor²)
  ```
- Both τ_trans and τ_tor are arrays over span (per station)
- Transverse shear split between webs by stiffness: V_i = V_web · (GA_i / Σ GA_j)
- Ensures both load paths checked against single allowable and buckling limit

---

### 4.4 Terminal Output Diagnostics

The tool prints extensive diagnostics during analysis. Key messages:

#### 4.4.1 Load Normalization

```
[LOAD-NORM] Target half-wing lift @ g_limit: 155.324 N
[LOAD-NORM] Integrated half-wing lift (after scale): 155.324 N
[LOAD-NORM] Applied scale factor: 1.234567
```

**Verify:**
- Target lift should equal (AUW × g_limit / 2)
- Integrated lift should match target to within 0.01 N
- Scale factor is applied to make the match exact

#### 4.4.2 Geometry Checks

```
DEBUG: chord=0.106680 m, span=1.3807 m, sweep=0.000 deg
DEBUG: web positions x1=0.021336 m, x2=0.069342 m
DEBUG: t_face=0.00026 m, t_core=0.003 m, t_web=0.00026 m
```

**Verify:**
- Chord matches input
- Web positions reasonable (x1 < x2, both within chord)
- Thicknesses physical (face ~ 0.1-1 mm, core ~ 2-10 mm)

#### 4.4.3 Neutral Axis and Section Properties

```
=== TRANSFORMED SECTION PROPERTIES ===
Neutral axis position: z_NA = 0.0123 mm (from geometric midplane)
Total transformed I: I_total_trans = 1.23456e-08 m^4
Implied E_avg = 89.2 GPa
```

**Verify:**
- z_NA near zero for symmetric sandwich (offset indicates asymmetry)
- I_total_trans scales with chord^4 roughly
- E_avg between core E (36 MPa) and fiber E (134 GPa), closer to fiber

#### 4.4.4 FoS Statistics

```
=== FACTOR OF SAFETY DIAGNOSTICS ===
FoS_skin_top: min=8.234, mean=15.456
FoS_rod_top: min=3.123, mean=7.890
FoS_web_shear: min=2.456, mean=5.678
FoS_core_shear: min=1.234, mean=3.456
```

**Verify:**
- All min values > 1.0 (no predicted failure)
- Core shear often lowest (common)
- Mean >> min indicates localized critical region (usually root)

#### 4.4.5 Governing Failure Analysis

```
=== GOVERNING FAILURE ANALYSIS ===
  Rod bottom          : FoS = 2.456
  Web shear           : FoS = 3.123
  Core shear          : FoS = 1.234  ← MINIMUM
  Skin torsion        : FoS = 8.901
Minimum FoS: 1.234 (mode: Core shear)
```

**Interpret:**
- Lowest FoS mode is governing
- If < 1.5, requires design changes
- If mode changes with sweep, indicates complex load interaction

---

## 5. Customization Examples

### 5.1 Changing Material Properties

#### Example 5.1.1: Higher Modulus Fiber (M55J instead of T700S)

**Edit `Materials` dataclass:**

```python
mat = Materials(
    # M55J fiber properties (540 GPa fiber, ~300 GPa UD at 60% Vf)
    E1_skin=300e9,     # Pa (was 134e9 for T700S)
    E2_skin=10e9,      # Pa (transverse unchanged)
    G12_skin=8.0e9,    # Pa (higher shear modulus)
    nu12_skin=0.30,
    
    # Strength slightly lower than T700S for high modulus fiber
    Xt_skin=2.0e9,     # Pa (was 2.86e9)
    Xc_skin=1.2e9,     # Pa (was 1.45e9)
    tau_skin_allow=100e6,  # Pa (was 136e6)
    
    # Rod properties (same fiber)
    E_rod=300e9,       # Pa
    Xt_rod=2.0e9,      # Pa
    Xc_rod=1.2e9,      # Pa
)
```

**Expected effect:**
- Higher stiffness → lower strains, lower curvature
- Lower strength → lower FoS if stress is similar
- Net effect depends on load case (stiffness often dominates)

#### Example 5.1.2: Denser Core (ROHACELL 110 instead of 31)

```python
mat = Materials(
    t_core=0.003,          # m (same thickness)
    rho_core=110.0,        # kg/m³ (was 32)
    E_core=180e6,          # Pa, 180 MPa (was 36 MPa)
    G_core=70e6,           # Pa, 70 MPa (was 13 MPa)
    tau_core_allow=2.3e6,  # Pa, 2.3 MPa (was 0.4 MPa)
)
```

**Expected effect:**
- Much higher core shear capacity (2.3 vs 0.4 MPa)
- Core shear FoS increases ~5.75×
- Weight increases ~3.4× (but core is small fraction of total)
- Core shear unlikely to govern

---

### 5.2 Adjusting Layup Stacks

#### Example 5.2.1: Thicker Skins (2 layers of ±45° biax)

```python
mat = Materials(
    face_stack=[(45.0, 2), (-45.0, 2)],  # was [(45, 1), (-45, 1)]
    # Total skin thickness: 4 plies × 0.13 mm = 0.52 mm (was 0.26 mm)
)
```

**Expected effect:**
- Bending stiffness EI increases significantly (t_skin enters squared in z²)
- Skin stress decreases (more area to share load)
- Weight increases (more material)
- Core shear may become more critical (relative to bending strength)

#### Example 5.2.2: 0°/±45° Hybrid Layup

```python
mat = Materials(
    face_stack=[(0.0, 1), (45.0, 1), (-45.0, 1)],
    # 0° plies provide axial stiffness, ±45° provide shear resistance
)
```

**Expected effect:**
- Axial (0°) stiffness very high in 0° plies
- Bending stiffness increases due to high E1 contribution
- CLT correctly accounts for orientation effects

#### Example 5.2.3: Quasi-Isotropic Layup

```python
mat = Materials(
    face_stack=[(0.0, 1), (45.0, 1), (-45.0, 1), (90.0, 1)],
    # Balanced properties in all directions
)
```

**Expected effect:**
- More isotropic behavior (less orientation-dependent)
- Heavier (more plies)
- Good for uncertain load directions

---

### 5.3 Custom Airfoil Input

#### Example 5.3.1: Creating Custom Airfoil File

**NACA 2412 example (`naca2412.txt`):**

Use an airfoil coordinate generator (e.g., XFOIL, online tools) to export x,y coordinates.

**Format:**
```
1.0000    0.0000
0.9900    0.0036
0.9800    0.0071
...
0.0000    0.0000  ← leading edge
...
0.9900   -0.0012
1.0000    0.0000
```

**Requirements:**
- At least 50 points for smooth outline
- LE near x=0, TE at x=1
- Closed loop (or will be closed automatically)

**Usage:**
```bash
python Wing.py --chord_m 0.12 --airfoil naca2412.txt
```

#### Example 5.3.2: Dimensional vs. Normalized Coordinates

**The tool accepts both formats:**

**Normalized (0-1):**
```
0.0000  0.0000
0.0100  0.0145
0.0200  0.0205
...
```

**Dimensional (e.g., 100mm chord):**
```
0.0000  0.0000
1.0000  1.4500
2.0000  2.0500
...
```

The tool normalizes automatically by detecting chord = xmax − xmin.

---

### 5.4 Modifying Flight Envelope

#### Example 5.4.1: High-Speed Configuration

**Edit configuration section:**

```python
VELOCITY_MS = 30.0          # m/s (was 15.87, double speed)
AIR_DENSITY = 1.0573        # kg/m³ (unchanged)
LIFT_COEFFICIENT = 0.8      # dimensionless (lower CL at higher speed)
G_LIMIT = 3.0               # load factor (reduced for high-speed case)
```

**Effect on dynamic pressure:**
```
q_old = 0.5 × 1.0573 × 15.87² = 133 Pa
q_new = 0.5 × 1.0573 × 30.0² = 476 Pa  (3.6× higher)
```

**Net lift change depends on CL and g_limit:**
```
L ∝ q × CL × g_limit
L_new / L_old = (476/133) × (0.8/1.43) × (3.0/5.0) = 1.20
```

So loads increase 20% despite lower g-limit.

#### Example 5.4.2: High-G Aerobatic Configuration

```python
VELOCITY_MS = 25.0
AIR_DENSITY = 1.0573
LIFT_COEFFICIENT = 1.6      # Near max CL
G_LIMIT = 12.0              # Aerobatic pull-up
```

**Expected result:**
- Very high loads (12g)
- Low FoS everywhere, likely < 1.5 in critical regions
- Indicates need for heavier structure

---

## 6. Troubleshooting

### 6.1 Common Errors and Fixes

#### Error 6.1.1: Missing Chord

**Error message:**
```
ValueError: Chord (inp.chord_m) must be provided and > 0.
```

**Cause:** `--chord_m` not specified on command line.

**Fix:**
```bash
python Wing.py --chord_m 0.10668
```

#### Error 6.1.2: Web Order Violation

**Error message:**
```
ValueError: web2_p must be greater than web1_p (webs must be ordered root->tip in x/c).
```

**Cause:** Rear web position ≤ front web position.

**Fix:** Edit configuration to ensure web2 > web1:
```python
WEB1_POSITION = 0.2
WEB2_POSITION = 0.65  # must be > 0.2
```

#### Error 6.1.3: Airfoil File Not Found

**Error message:**
```
FileNotFoundError: Airfoil file not found: my_airfoil.txt
```

**Cause:** File doesn't exist in current directory.

**Fix:**
- Verify file exists: `ls my_airfoil.txt` (Linux/Mac) or `dir my_airfoil.txt` (Windows)
- Provide full path: `--airfoil /full/path/to/my_airfoil.txt`
- Or copy file to current directory

#### Error 6.1.4: Invalid Web Height

**Warning message:**
```
ValueError: Web height invalid: h1_web=0.000mm, h2_web=0.000mm.
Skins (t_face=0.500mm) may be too thick for airfoil.
```

**Cause:** Skin thickness ≥ half of airfoil thickness at web location.

**Diagnosis:**
- Airfoil too thin at web positions
- Skins too thick
- Web positions in very thin region (near LE or TE)

**Fix (choose one):**
1. Reduce skin thickness: Fewer plies or thinner fabric
2. Move webs to thicker region: Adjust WEB1_POSITION, WEB2_POSITION
3. Use thicker airfoil: Scale chord or use fatter airfoil section
4. Reduce core thickness: Makes room for skins

#### Error 6.1.5: Shapely Import Error

**Error message:**
```
ImportError: shapely is required. Install with: pip install shapely
```

**Fix:**
```bash
pip install shapely
```

If installation fails, try:
```bash
pip install shapely --upgrade
```

or use conda:
```bash
conda install shapely
```

---

### 6.2 Validation Checks

#### Check 6.2.1: Load Integration

**What to verify:**

Terminal output shows:
```
[LOAD-NORM] Target half-wing lift @ g_limit: 155.324 N
[LOAD-NORM] Integrated half-wing lift (after scale): 155.324 N
```

**Acceptance criteria:**
- Integrated value matches target to < 0.1 N

**If they don't match:**
- Likely bug in code (report issue)
- Or numerical integration tolerance too loose

#### Check 6.2.2: Bredt Solver Test

**What to verify:**

Terminal output at startup:
```
Bredt single-cell test status: bredt_solved cond: 124.5
  tau analytic: 500000.0
  tau solver: [500000. 500000. 500000. 500000.]
```

**Acceptance criteria:**
- Status: `bredt_solved` (not `ill_conditioned`)
- Condition number < 1e6
- tau_solver matches tau_analytic to within 1%

**If test fails:**
- Likely code error or numerical issue
- Check for updates or report bug

#### Check 6.2.3: FoS Sanity

**What to verify:**

All minimum FoS values should be positive and finite.

```
FoS_skin_top: min=5.234, mean=12.456
FoS_web_shear: min=2.345, mean=8.901
```

**Red flags:**
- Any min FoS < 0 → Sign error in code
- Any min FoS = inf → Zero stress somewhere (unusual)
- Any min FoS = NaN → Numerical issue

---

### 6.3 Numerical Issues

#### Issue 6.3.1: Ill-Conditioned Bredt Matrix

**Symptom:**

Terminal output shows:
```
Bredt multi-cell status: ill_conditioned, cond: 1.23e15
```

**Cause:**
- Very disparate wall thicknesses or shear moduli
- Near-zero cell areas
- Degenerate geometry (webs too close together)

**Fix:**
1. Check web positions are reasonable (not too close)
2. Verify material properties are physical
3. Increase web thickness if very thin
4. Consider single-cell simplification if appropriate

#### Issue 6.3.2: Zero or Negative FoS

**Symptom:**

Excel shows FoS = 0 or negative values.

**Cause:**
- Allowable = 0 (check material properties)
- Applied stress has wrong sign
- Division by zero in FoS calculation

**Fix:**
- Review material allowables (all should be positive)
- Check for recent code changes
- Report as potential bug

---

### 6.4 Performance Issues

#### Issue 6.4.1: Slow Analysis

**Symptom:**

Analysis takes > 1 minute per sweep.

**Diagnosis:**

Check `N_STATIONS` in configuration:
```python
N_STATIONS = 101  # default
```

**Solutions:**

1. **Reduce stations** (faster, less resolution):
```python
N_STATIONS = 51  # 2× faster
```

2. **Skip plots** (saves ~30% time):
```bash
python Wing.py --chord_m 0.1 --no_plots
```

3. **Single sweep first** (verify before multi-sweep):
```bash
python Wing.py --chord_m 0.1 --sweeps 0
```

#### Issue 6.4.2: Large Output Files

**Symptom:**

Excel file > 50 MB.

**Cause:**

Many sweeps × many stations = huge `2_All_Sweeps` sheet.

**Solution:**

Reduce data volume:
- Fewer sweep angles
- Fewer stations (N_STATIONS = 51 instead of 101)
- Delete per-sweep sheets after viewing (keep `2_All_Sweeps` only)

---

## 7. Advanced Topics

### 7.1 Rod Position Optimization

**Current implementation:** 4 rods at fixed positions (2 forward, 2 aft).

**To change rod positions, edit line 1156:**

```python
rod_positions_xc = [0.2, 0.2, 0.65, 0.65]  # fraction of chord
```

**Example: Concentrate rods near front spar:**
```python
rod_positions_xc = [0.18, 0.22, 0.18, 0.22]  # all near 20% chord
```

**Effect:**
- Rods closer to webs → less rod moment arm → potentially lower web shear Q
- But concentrating area may increase local stresses

**Optimization approach:**
1. Run baseline (current positions)
2. Try concentrated layout (rods near single location)
3. Try distributed layout (rods spread across chord)
4. Compare FoS results and choose best

---

### 7.2 Multi-Cell Topology Understanding

**Current cell structure:**

3 cells created by 2 webs:

```
LE -------|web1|-------|web2|------- TE
   Cell 0        Cell 1       Cell 2
```

**Walls in Bredt system:**

1. Top skin (shared between cells 0-1 and 1-2)
2. Front web (internal to cell 1)
3. Middle skin (shared between cells 1-2)
4. Rear web (internal to cell 2)
5. Bottom skin (shared between cells 0 and 2)

**Effect of web spacing:**
- Narrow cells → higher torsional stiffness (smaller A)
- Wide cells → lower torsional stiffness
- Optimal spacing depends on loading

---

### 7.3 Parametric Study Design

**Example: Chord-Span Trade Study**

**Goal:** Find minimum weight configuration for fixed g-limit.

**Approach:**

1. **Define parameter ranges:**
```python
chords = [0.08, 0.10, 0.12, 0.14, 0.16]  # m
spans = [3.0, 4.0, 5.0, 6.0]             # ft
```

2. **Run matrix:**
```bash
for chord in 0.08 0.10 0.12 0.14 0.16; do
  for span in 3.0 4.0 5.0 6.0; do
    python Wing.py --chord_m $chord --span_ft $span --out "results_c${chord}_s${span}.xlsx" --no_plots
  done
done
```

3. **Extract results:**
Open each Excel file, record:
- Weight (from `Weight_Summary` sheet)
- Min FoS (from `0_Final_Summary`)
- Governing mode (from `0_Final_Summary`)

4. **Filter feasible designs:**
Keep only designs with Min FoS ≥ 1.5.

5. **Select minimum weight:**
From feasible set, choose lowest weight.

**Plot results:**
- Weight vs. chord (for each span)
- FoS vs. chord
- Identify optimal region

---

### 7.4 Interpreting Sweep Effects

**Physical mechanism:**

Sweep angle Λ projects the lift force:
```
L_flapwise = L_total × cos(Λ)
L_inplane = L_total × sin(Λ)
```

**Effect on loads:**
- Flapwise bending M decreases with sweep (∝ cos Λ)
- Torsion T decreases if moment arm also projected
- In-plane force increases (not currently modeled)

**FoS trends with sweep:**
- Bending-critical designs: FoS increases with sweep (good)
- Torsion-critical designs: FoS may decrease if torsion mechanism changes
- Shear-critical designs: FoS increases with sweep (V decreases)

**Typical result:**
- Min FoS increases with sweep for most configurations
- But high sweep (> 50°) may have other issues (not modeled: in-plane bending, aeroelasticity)

---

### 7.5 CLT Implementation Details

**ABD Matrix Construction:**

For a laminate stack, the tool computes:

```
A = ∫ Qbar dz     (extensional stiffness)
B = ∫ Qbar z dz   (bending-extension coupling)
D = ∫ Qbar z² dz  (bending stiffness)
```

Where Qbar is the transformed ply stiffness matrix.

**Effective Bending Stiffness:**

```python
E_eff_bending = 12 × D11 / t³
```

Used for skin stress scaling in bending.

**Membrane Shear Modulus:**

```python
G_membrane = A66 / t
```

Used for web and skin shear stiffness in torsion.

**Advantage of CLT:**

- Correctly accounts for ply orientations (0°, ±45°, 90°)
- Handles unsymmetric laminates (B matrix)
- Provides physical basis for composite properties

---

## Appendix A: Quick Reference Tables

### A.1 Typical Material Properties

| Material | E (GPa) | ρ (kg/m³) | Cost/Perf | Application |
|----------|---------|-----------|-----------|-------------|
| T700S UD (60% Vf) | 134 | 1500 | Moderate | Standard high-perf |
| M55J UD (60% Vf) | 300 | 1600 | High | High stiffness |
| IM7 UD (60% Vf) | 160 | 1550 | Moderate | High strength |
| ROHACELL 31 IG-F | 0.036 | 32 | Low | Light load |
| ROHACELL 71 | 0.105 | 75 | Moderate | Medium load |
| ROHACELL 110 | 0.180 | 110 | High | Heavy load |

### A.2 Typical Layup Sequences

| Layup | Orientation | Purpose | EI Benefit | Weight |
|-------|-------------|---------|------------|--------|
| ±45° biax | [45/-45] | Shear, light | Low | Low |
| ±45° biax double | [45/-45]₂ | Shear, moderate | Moderate | Moderate |
| 0°/±45° | [0/45/-45] | Bending-optimized | High | Moderate |
| ±45°/0° | [45/-45/0] | Alternative | High | Moderate |
| Quasi-isotropic | [0/45/-45/90] | All-around | Moderate | High |

### A.3 Design Guidelines

| Parameter | Typical Range | Conservative | Aggressive |
|-----------|---------------|--------------|------------|
| Chord (m) | 0.08 - 0.20 | > 0.12 | < 0.10 |
| Span/Chord | 10 - 30 | < 15 | > 25 |
| Core (mm) | 2 - 8 | > 5 | < 3 |
| Skin plies | 2 - 8 | > 4 | < 3 |
| Rod count | 2 - 8 | > 4 | < 4 |
| FoS target | 1.5 - 3.0 | > 2.5 | < 2.0 |

---

## Appendix B: File Formats

### B.1 Airfoil Coordinate File

**Format:** Plain text, two columns (space or tab delimited)

**Example:**
```
# Optional comment lines starting with #
1.00000  0.00000
0.99500  0.00123
0.99000  0.00245
...
0.00000  0.00000
...
1.00000  0.00000
```

**Requirements:**
- Columns: x (chordwise), y (thickness direction)
- Points: At least 20, recommended 100-200
- Order: Start at TE, go around upper surface, through LE, back along lower surface to TE
- Closure: First point should equal last point (or will be closed automatically)
- Normalization: Will be auto-normalized to chord = 1

---

## Appendix C: Glossary

| Term | Definition |
|------|------------|
| ABD matrix | Laminate stiffness matrix (A=extensional, B=coupling, D=bending) |
| Biax | Fabric with two fiber directions (e.g., ±45°) |
| Bredt-Batho | Theory for torsion of thin-walled multi-cell sections |
| CLT | Classical Laminate Theory |
| FoS | Factor of Safety (allowable / applied) |
| GA | Shear stiffness (G × Area) |
| g-limit | Maximum design load factor (multiples of 1g) |
| LE / TE | Leading Edge / Trailing Edge |
| NA | Neutral Axis (zero bending strain location) |
| Sandwich | Structure with stiff skins and lightweight core |
| Spar cap | Spanwise rod or beam flange carrying bending |
| UD | Unidirectional (all fibers parallel) |
| VQ/It | Classical beam shear stress formula |
| Web | Internal vertical wall carrying shear |

---

*For theoretical background, see [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)*

*For validation procedures, see [VALIDATION.md](VALIDATION.md)*


