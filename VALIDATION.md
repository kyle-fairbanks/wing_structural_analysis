# Wing Structural Analysis - Validation Guide

Test cases, verification procedures, and quality assurance for composite wing analysis.

**Version 1.0**

---

## Table of Contents

1. [Built-In Validation](#1-built-in-validation)
2. [Recommended Verification Workflow](#2-recommended-verification-workflow)
3. [Example Test Case](#3-example-test-case)
4. [Known Limitations](#4-known-limitations)
5. [When to Seek Further Analysis](#5-when-to-seek-further-analysis)
6. [Quality Assurance Checklist](#6-quality-assurance-checklist)

---

## 1. Built-In Validation

### 1.1 Automatic Validation Tests

The tool runs validation tests automatically on startup. Review terminal output to verify correct operation.

#### 1.1.1 Bredt Torsion Solver Test

**Location:** Runs at program start (line 2864 in `Wing.py`)

**Test case:**
- Rectangular box section: 200 mm × 100 mm × 1 mm wall thickness
- Uniform shear modulus: G = 1 GPa
- Applied torque: T = 10 N·m

**Expected output (after Fix #1):**
```
Bredt single-cell test status: bredt_solved cond: 1.0

tau analytic: 250000.0

tau solver: [250000. 250000. 250000. 250000.]

```

**Acceptance criteria:**
- Status: `bredt_solved`
- **cond: 1.0 (dummy value for nc==1, not a conditioning measure)**
- Accuracy: **exact match (not just 1% tolerance)**

**Note**: For nc==1, no matrix is solved (direct calculation), so conditioning is not applicable.

**If test fails:**
- Check for code changes in Bredt solver functions
- Verify numpy/scipy versions are compatible
- Report issue with full terminal output

---

### 1.2 Diagnostic Output Interpretation

#### 1.2.1 Load Normalization Check

**What to look for:**

```
[LOAD-NORM] Target half-wing lift @ g_limit: 155.324 N
[LOAD-NORM] Integrated half-wing lift (after scale): 155.324 N
[LOAD-NORM] Applied scale factor: 1.234567
```

**Verification:**

| Item | Expected | Actual | Status |
|------|----------|--------|--------|
| Target lift | = 0.5 × AUW × g_limit | ___ N | ✓ / ✗ |
| Integrated lift | Matches target ± 0.1 N | ___ N | ✓ / ✗ |
| Scale factor | 0.5 - 2.0 (reasonable) | ___ | ✓ / ✗ |

**Interpretation:**
- Scale factor ≈ 1.0: Initial lift shape well-matched to target
- Scale factor >> 1.0: Initial shape too weak, scaled up
- Scale factor << 1.0: Initial shape too strong, scaled down
- Scale factor outside [0.1, 10]: Potential issue, investigate

#### 1.2.2 Geometry Validation

**What to look for:**

```
DEBUG: chord=0.106680 m, span=1.3807 m, sweep=0.000 deg
DEBUG: web positions x1=0.021336 m, x2=0.069342 m
DEBUG: t_face=0.00026 m, t_core=0.003 m, t_web=0.00026 m
```

**Verification checklist:**

- [ ] Chord matches input (within rounding)
- [ ] Span = input_ft × 0.3048 (for full wing, multiply by 2)
- [ ] x_web1 = web1_p × chord
- [ ] x_web2 = web2_p × chord
- [ ] x_web2 > x_web1 (rear web aft of front web)
- [ ] t_face = n_plies × ply_thickness (check calculation)
- [ ] t_core matches input
- [ ] All thicknesses in reasonable range (0.1 - 10 mm)

#### 1.2.3 Section Properties Check

**What to look for:**

```
=== TRANSFORMED SECTION PROPERTIES ===
Neutral axis position: z_NA = 0.0123 mm (from geometric midplane)
Total transformed I: I_total_trans = 1.23456e-08 m^4
Implied E_avg = 89.2 GPa
```

**Verification:**

| Property | Typical Range | Your Value | Status |
|----------|---------------|------------|--------|
| z_NA | ±2 mm (symmetric sandwich) | ___ mm | ✓ / ✗ |
| I_total_trans | 1e-9 to 1e-7 m⁴ (small wings) | ___ m⁴ | ✓ / ✗ |
| E_avg | 50 - 120 GPa (fiber-dominated) | ___ GPa | ✓ / ✗ |

**Red flags:**
- |z_NA| > 5 mm: Highly asymmetric section (unusual)
- E_avg < 10 GPa: Core-dominated (check if rod/skin properties correct)
- E_avg > 150 GPa: Suspiciously high (check for input errors)

#### 1.2.4 FoS Summary Statistics

**What to look for:**

```
=== FACTOR OF SAFETY DIAGNOSTICS ===
FoS_skin_top: min=8.234, mean=15.456
FoS_rod_top: min=3.123, mean=7.890
FoS_web_shear: min=2.456, mean=5.678
FoS_core_shear: min=1.234, mean=3.456
```

**Sanity checks:**

- [ ] All minimum FoS > 0 (no negative values)
- [ ] All minimum FoS < 1e5 (no infinities)
- [ ] Mean FoS > Min FoS (expected, loads highest at root)
- [ ] Mean / Min ratio: 1.5 - 5× typical
- [ ] No NaN values

**Typical patterns:**
- Core shear: Often lowest FoS (ROHACELL 31 is weak in shear)
- Rod stresses: Often second-lowest (high bending moment at root)
- Skin torsion: Usually high FoS (skins are strong in shear)

---

## 2. Recommended Verification Workflow

### 2.1 Pre-Analysis Checks

Before running analysis, verify inputs:

#### Step 1: Verify Geometry Inputs

```python
# Check configuration section (lines 2786-2810)
CHORD_M = 0.10668  # ✓ Must be provided
WEB1_POSITION = 0.2  # ✓ 0 < web1 < 1
WEB2_POSITION = 0.65  # ✓ web2 > web1
ROD_WIDTH_MM = 3.0  # ✓ Reasonable (1-10 mm)
ROD_HEIGHT_MM = 3.0  # ✓ Fits in core (< t_core)
```

**Command-line check:**
```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53 --sweeps 0
```

Expected: No errors, analysis starts.

#### Step 2: Verify Material Properties

```python
# Check Materials dataclass
mat = Materials(
    face_stack=[(45, 1), (-45, 1)],  # ✓ Valid ply angles
    t_core=0.003,  # ✓ Physical (2-10 mm)
    E1_skin=134e9,  # ✓ Reasonable for T700S
    tau_core_allow=0.4e6,  # ✓ Matches ROHACELL 31
)
```

**Red flags:**
- E values in wrong units (should be Pa, not GPa)
- Allowables < moduli (unphysical)
- Density = 0 (causes mass = 0)

#### Step 3: Verify Flight Conditions

```python
VELOCITY_MS = 15.87  # ✓ Reasonable for RC/UAV
AIR_DENSITY = 1.0573  # ✓ Sea level ISA
LIFT_COEFFICIENT = 1.43  # ✓ Near max CL
G_LIMIT = 5.0  # ✓ Typical maneuvering load
```

**Sanity check dynamic pressure:**
```
q = 0.5 × ρ × V² = 0.5 × 1.0573 × 15.87² = 133 Pa
```

Typical range: 50 - 500 Pa for small RC aircraft.

---

### 2.2 Post-Analysis Verification

After analysis completes, verify results are physically reasonable.

#### Step 1: Check Load Integration

**Terminal output:**
```
[LOAD-NORM] Target half-wing lift @ g_limit: 155.324 N
[LOAD-NORM] Integrated half-wing lift (after scale): 155.324 N
```

**Manual check:**
```
AUW = 6.33 lb × 4.448 N/lb = 28.15 N
Target = 0.5 × 28.15 × 5.0 = 70.4 N  (half-wing at 5g)
```

✓ Values should match. If not, investigate load calculation.

#### Step 2: Verify Root Loads (Equilibrium)

**Excel:** Open `2_All_Sweeps` sheet, filter sweep_deg = 0, look at i = 0 (root).

**Expected:**
- `V(0)` ≈ Target half-wing lift (within 1%)
- `M(0)` ≈ V × (span/3 to span/4) for typical distributions
- `T(0)` > 0 if moment arm ≠ 0

**Example:**
```
V_root = 155 N
span_half = 0.69 m
M_root ≈ 155 × 0.69/3 ≈ 35 N·m  (rough estimate)
```

Actual M depends on load distribution shape.

#### Step 3: Check FoS Trends

**Excel:** Plot FoS vs. x (spanwise position).

**Expected trends:**
- FoS **decreases** from tip to root (loads increase)
- Minimum FoS near root (0 < x < span/4)
- FoS >> minimum at tip (low loads)

**Red flags:**
- FoS increases toward root (load direction wrong?)
- Minimum FoS at mid-span (unusual, investigate)
- FoS constant across span (load not varying?)

#### Step 4: Compare Weight Estimate

**Excel:** Open `Weight_Summary` sheet.

**Typical values for 4.5 ft wing, 0.1 m chord:**

| Component | Typical Mass | Your Value | Ratio |
|-----------|--------------|------------|-------|
| Skins | 60-120 g | ___ g | ___ |
| Core | 20-40 g | ___ g | ___ |
| Webs | 5-15 g | ___ g | ___ |
| Rods | 30-60 g | ___ g | ___ |
| **Total** | **120-240 g** | **___ g** | **___** |

**If weight is way off (> 2× expected):**
- Check density values (wrong units?)
- Check planform area calculation
- Check ply counts (too many layers?)

#### Step 5: Review Governing Failure Mode

**Excel:** `0_Final_Summary` sheet, note `Governing failure mode`.

**Common modes and interpretation:**

| Mode | Interpretation | Action |
|------|----------------|--------|
| Core shear | Foam is limiting factor | Use denser core or add webs |
| Rod top/bottom | Bending-critical | Add rods, increase size, or thicken skins |
| Web shear | Shear-critical | Thicken webs or add more webs |
| Skin torsion | Torsion-critical | Optimize web spacing, reduce moment arm |

**Red flags:**
- Mode changes drastically between similar cases (check for errors)
- Unexpected mode (e.g., "Skin torsion" for no-sweep case with zero torque)

---

### 2.3 Consistency Checks Across Sweeps

If running multi-sweep study:

#### Check 1: Load Projection

**Theory:** Flapwise loads ∝ cos(Λ)

**Test:**
```
M_root(30°) / M_root(0°) ≈ cos(30°) = 0.866
```

**Excel:** Compare root moments from `Summaries_All_Sweeps`.

**Tolerance:** Within 5% (small deviations due to discretization).

#### Check 2: FoS Trends with Sweep

**Expected:**
- Bending FoS increases with sweep (M decreases)
- Shear FoS increases with sweep (V decreases)
- Overall min FoS typically increases with sweep

**If min FoS decreases with sweep:**
- Check if torsion is increasing (may offset bending benefit)
- Verify sweep projection is applied correctly

#### Check 3: Weight Independence

**Check:** Weight should be **identical** across all sweeps (geometry doesn't change).

**Excel:** `Weight_Summary` is computed once (sweep-independent).

If weight varies between runs, likely a bug (report).

---

## 3. Example Test Case

### 3.1 Baseline Configuration

**Objective:** Verify tool produces reasonable results for typical small RC wing.

#### 3.1.1 Input Parameters

**Geometry:**
```python
SPAN_FT = 4.53           # 4.53 ft full wingspan
CHORD_M = 0.10668        # 106.68 mm chord
WEB1_POSITION = 0.2      # 20% chord
WEB2_POSITION = 0.65     # 65% chord
ROD_WIDTH_MM = 3.0       # 3 mm × 3 mm rods
ROD_HEIGHT_MM = 3.0
ROD_GAP_MM = 0
```

**Materials:**
```python
face_stack = [(45, 1), (-45, 1)]  # One FCIM255 biax layer per face
web_stack = [(45, 1), (-45, 1)]   # One FCIM255 biax layer per web
t_core = 0.003 m  # 3 mm ROHACELL 31 IG-F
```

**Flight Conditions:**
```python
VELOCITY_MS = 15.87      # 35.5 mph
AIR_DENSITY = 1.0573     # Sea level ISA
LIFT_COEFFICIENT = 1.43  # Near max CL
G_LIMIT = 5.0            # 5g maneuvering load
```

**Run command:**
```bash
python Wing.py --chord_m 0.10668 --span_ft 4.53 --sweeps 0
```

---

#### 3.1.2 Expected Outputs

**Load integration:**
```
[LOAD-NORM] Target half-wing lift @ g_limit: ~155 N
[LOAD-NORM] Integrated half-wing lift (after scale): ~155 N
```

**Root loads (i = 0 in Excel):**

| Load | Expected Value | Tolerance |
|------|----------------|-----------|
| V_root | 150 - 160 N | ±10 N |
| M_root | 30 - 40 N·m | ±5 N·m |
| T_root | 0.5 - 2.0 N·m | ±1 N·m |

**Section properties:**
```
z_NA: 0.0 ± 1.0 mm  (symmetric)
EI_total: 5e-8 to 1e-7 N·m²  (depends on exact layup)
E_avg: 80 - 100 GPa  (fiber-dominated)
```

**Weight summary (full wing):**

| Component | Expected Mass | Tolerance |
|-----------|---------------|-----------|
| Skins | 80 - 100 g | ±20 g |
| Core | 25 - 35 g | ±5 g |
| Webs | 8 - 12 g | ±3 g |
| Rods | 40 - 50 g | ±10 g |
| **TOTAL** | **160 - 200 g** | **±30 g** |

**Factor of Safety:**

| Mode | Expected FoS | Status |
|------|--------------|--------|
| Skin top | 8 - 15 | Non-critical |
| Skin bottom | 8 - 15 | Non-critical |
| Rod top | 3 - 6 | Moderate |
| Rod bottom | 2 - 5 | Moderate |
| Web shear | 3 - 8 | Moderate |
| **Core shear** | **1.5 - 3.0** | **Often governs** |
| Skin torsion | 15 - 50 | Non-critical |

**Governing mode:** Likely **Core shear** (ROHACELL 31 is weak at 0.4 MPa).

**Min FoS @g_limit:** 1.8 - 2.5 typical (adequate for preliminary design, marginal for production).

---

#### 3.1.3 Interpretation Guidance

**If min FoS < 1.5:**
- Considered marginal
- **Actions:** Upgrade to ROHACELL 71 (0.7 MPa → 2.9 MPa shear), thicken core, or add webs

**If min FoS 1.5 - 3.0:**
- Adequate for preliminary design
- **Actions:** Proceed with detailed FEA, consider fatigue analysis, prototype testing

**If min FoS > 3.0:**
- Very conservative, overdesigned
- **Actions:** Reduce material (lighter/cheaper), or keep margin for safety/unknown loads

**Weight assessment:**
- 160-200 g for this wing is reasonable (5.6-7.1 oz)
- As fraction of AUW: 160 g / (6.33 lb × 454 g/lb) = 5.6% (typical for composite wings)

---

### 3.2 Parametric Variation Test

**Objective:** Verify tool responds correctly to parameter changes.

#### Test 3.2.1: Chord Scaling

**Vary chord:** 0.08, 0.10, 0.12, 0.14, 0.16 m (span fixed)

**Expected trends:**
- Weight ∝ chord (linear, area increases)
- EI ∝ chord⁴ (bending stiffness)
- Stresses ∝ 1/chord² (for same g-limit, M ∝ chord, z ∝ chord, EI ∝ chord⁴ → σ ∝ M×z/EI ∝ 1/chord²)
- FoS ∝ chord² (inverse of stress)

**Acceptance:**
- Doubling chord → ~4× higher FoS ✓
- Weight doubles ✓

#### Test 3.2.2: G-Limit Scaling

**Vary g_limit:** 2.5, 5.0, 7.5, 10.0 (all else fixed)

**Expected trends:**
- All loads (V, M, T) ∝ g_limit (linear)
- All stresses ∝ g_limit (linear)
- FoS ∝ 1/g_limit (inverse)
- Weight unchanged (geometry doesn't change)

**Acceptance:**
- Doubling g_limit → halves all FoS ✓
- Root V, M double ✓
- Weight identical ✓

#### Test 3.2.3: Core Upgrade

**Test:** Change from ROHACELL 31 to ROHACELL 110

```python
# Before (ROHACELL 31)
mat = Materials(
    rho_core=32.0,
    E_core=36e6,
    G_core=13e6,
    tau_core_allow=0.4e6,
)

# After (ROHACELL 110)
mat = Materials(
    rho_core=110.0,
    E_core=180e6,
    G_core=70e6,
    tau_core_allow=2.3e6,
)
```

**Expected changes:**
- Core shear FoS increases ~5.75× (2.3 / 0.4)
- Core mass increases ~3.4× (110 / 32)
- Total weight increases ~10-15% (core is ~15-20% of total)
- Bending stiffness EI increases slightly (~5%, core E matters little)
- If core shear was governing, new governing mode likely different

**Acceptance:**
- Core shear no longer critical ✓
- Weight increase is small fraction of total ✓
- Overall min FoS improves if core shear was limiting ✓

---

## 4. Known Limitations

### 4.1 Model Assumptions

**Understand these limitations before interpreting results:**

#### 4.1.1 Half-Wing Symmetry

**Assumption:** Symmetric flight (no roll, yaw, or side loads).

**Implication:**
- Valid for straight-and-level flight, symmetric pull-up
- **Invalid** for:
  - Aileron deflection (asymmetric loading)
  - Sideslip (side loads on vertical tail, induced roll)
  - Spin or spiral (complex 3D loads)

**Recommendation:** For asymmetric cases, analyze worst-case side (higher lift) with additional safety factor.

#### 4.1.2 Uniform Chord (No Taper)

**Assumption:** Constant chord from root to tip.

**Implication:**
- Real wings often have taper (tip chord < root chord)
- Taper affects:
  - Weight distribution (lighter toward tip)
  - Stiffness distribution (EI decreases toward tip)
  - Stress distribution (stress relief at tip due to lower loads)

**Recommendation:**
- Use mean chord for preliminary sizing
- For tapered wing, model conservatively (use minimum chord) or extend code to handle taper; irrelevant for current application

#### 4.1.3 Static Loading Only

**Assumption:** Loads are static (applied slowly, no dynamics).

**Implication:**
- **Invalid** for:
  - Impact loads (landing, gust, bird strike)
  - Vibration (flutter, buffet, propeller-induced)
  - Fatigue (repeated loading)

**Recommendation:**
- For dynamic loads, use dynamic amplification factors (DAF ≈ 1.2 - 1.5)
- For fatigue, perform separate life analysis with load spectra

#### 4.1.4 Perfect Bond (No Debonding)

**Assumption:** Skins perfectly bonded to core, no delamination.

**Implication:**
- Real wings can experience skin-core debond (especially at high peel stresses)
- Tool does not check peel/through-thickness normal stress

**Recommendation:**
- Use proven manufacturing methods (vacuum bag, autoclave)
- Consider peel stress analysis for critical regions (high curvature, load introduction)

---

### 4.2 Analysis Limitations

#### 4.2.1 ~~No Skin Buckling Analysis~~ Sandwich Wrinkling Included (First-Order Screen)

**Status**: Sandwich face-sheet wrinkling now checked (as of Fix #6).

**Current implementation:**
- Cube-root formula: `σ_wr = c·(E_face·E_core·G_core/(1-ν²))^(1/3)`
- Uses CLT in-plane A11/t for E_face (accounts for ±45° orientation)
- Coefficient c = 0.5 (typical range 0.4-0.6, ±20% uncertainty)
- Proper physics for continuous core support (not Euler column buckling)

**What it captures:**
- Local face-sheet wrinkling (correct mode for sandwich with foam core)
- Core stiffness contribution (E_core, G_core both matter)
- Approximate effective modulus from CLT

**Limitations:**
- First-order screen, not definitive stability allowable
- Assumes symmetric wrinkling (top/bottom faces same)
- Coefficient uncertainty: results may vary ±20% with literature formulas
- Does NOT check:
  - Core shear crimping (separate mode, usually non-critical)
  - Post-wrinkling strength (conservative: assumes failure at wrinkling)
  - Intracell dimpling (very short wavelength, higher critical stress)

**When detailed buckling FEA needed:**
- FoS_wrinkling < 1.5 (within uncertainty band)
- Very thin skins (t < 0.2 mm) with weak core
- High compression (σ > 500 MPa)
- Asymmetric sandwich (different top/bottom faces)

#### 4.2.2 Simplified Core Shear Model

**Limitation:** Core shear area assumed as `b_eff × t_core` with b_eff = 0.85 × chord.

**Implication:**
- Actual core shear distribution is non-uniform (higher near webs)
- 15% reduction is empirical, not derived from first principles

**Sensitivity:**
- ±10% variation in b_eff → ±10% variation in core shear FoS
- If core shear is critical, verify with FEA

**Recommendation:**
- For heavily-loaded wings (core shear FoS < 2), validate with FEA
- Consider local reinforcement at root (doublers, thicker core)

#### 4.2.3 No Aeroelastic Effects

**Limitation:** Wing shape is rigid (no deformation under load).

**Implication:**
- Real wings deflect, twist, and deform → changes aerodynamic loads
- Flutter (dynamic instability) not analyzed
- Divergence (static instability at high speed) not checked

**When aeroelastic effects matter:**
- Flexible wings (low EI, high span)
- High speeds (V > 30 m/s for small wings)
- High aspect ratio (span/chord > 15)

**Recommendation:**
- For flexible wings, use aeroelastic analysis tools (NASTRAN, ZAERO)
- Check flutter speed with simplified methods (Goland wing model)

#### 4.2.4 Linear Elastic Material

**Limitation:** No plasticity, no damage progression, no nonlinear effects.

**Implication:**
- Real composites show nonlinear behavior near failure (matrix cracking, fiber breakage)
- Tool predicts failure at first-ply failure (conservative)
- No progressive failure analysis (redistribution after local failure)

**Recommendation:**
- FoS > 1.5 provides margin for nonlinear effects
- For critical designs, use progressive failure codes (ANSYS ACP, Abaqus/Explicit)

#### 4.2.5 Shear Model Choice: GA Partition vs. VQ/It

**Current implementation**: GA stiffness-based partition ("parallel springs") with stiffness-weighted web split

**What it means:**
- Transverse shear V splits between webs and core by stiffness ratio: GA_i / Σ(GA_j)
- Webs split by individual stiffness: V_web_i = V_web_total · (GA_web_i / Σ GA_web)
- Web stress computed directly: τ_web_i = V_web_i / (h_web_i · t_web_i)
- Core stress computed directly: τ_core = V_core / (b_eff · t_core)
- Does NOT use VQ/It (beam theory shear flow)

**Why this choice:**
- Correct for sandwich with vastly different web/core shear moduli (G_web ~ 5 GPa, G_core ~ 13 MPa)
- Avoids circular logic of extracting forces from stress fields
- Computationally simple and physically transparent
- Handles asymmetric webs robustly (different heights, materials, damage)

**Alternative (not implemented):**
- Full thin-wall shear flow distribution: q(s) = ∫(VQ/I)ds around closed section
- Would give similar results (webs much stiffer → attract most flow)
- More complex to implement correctly
- Primarily useful for detailed stress distribution around perimeter

**Impact**: GA model is conservative for core, slightly unconservative for webs (but margin is small).

---

### 4.3 Geometry Limitations

#### 4.3.1 Rectangular Planform Only

**Limitation:** Span and chord are constant (no sweep-back, taper, twist, or dihedral in planform).

**Note:** Sweep angle affects **loading** (cos Λ projection) but not **planform shape** (wing is still rectangular in plan view).

**Implication:**
- Real wings often have:
  - Sweep-back (LE/TE swept, not just load projection)
  - Taper (tip chord < root chord)
  - Twist (washout, tip at lower angle of attack)
  - Dihedral (tip higher than root)

**Recommendation:**
- For complex planforms, use 3D FEA or extend code

#### 4.3.2 Fixed Airfoil Section

**Limitation:** Airfoil shape is constant along span.

**Implication:**
- Real wings may have varying airfoil (different section at root vs. tip)

**Recommendation:**
- Use representative airfoil (e.g., root section for strength-critical analysis)

---

## 5. When to Seek Further Analysis

### 5.1 FoS Thresholds

**Use these guidelines to decide when detailed FEA or testing is required:**

| FoS Range | Status | Action |
|-----------|--------|--------|
| **> 3.0** | Conservative | Preliminary design adequate, consider optimization |
| **2.0 - 3.0** | Adequate | Acceptable for production with testing |
| **1.5 - 2.0** | Marginal | Proceed with caution, FEA recommended |
| **1.0 - 1.5** | Inadequate | Redesign or detailed FEA **required** |
| **< 1.0** | Failure predicted | **Do not build**, redesign immediately |

### 5.2 When FEA is Recommended

**Consider detailed finite element analysis if:**

1. **Low FoS:** Min FoS < 1.5 anywhere
2. **Thin skins:** Face thickness < 0.5 mm with high compression
3. **Complex geometry:** Load introduction, cutouts, attachments
4. **High loads:** g_limit > 8 or dynamic impact loads
5. **Flexible wing:** Large deflections (δ_tip > span/20)
6. **Critical application:** Human-rated, high consequence of failure

**What FEA provides:**
- 3D stress distribution (not just beam theory)
- Skin buckling analysis (stability check)
- Local stress concentrations (holes, corners, load points)
- Validated with industry-standard tools (ANSYS, Abaqus, Nastran)

### 5.3 When Physical Testing is Required

**Test requirements (regulatory or best practice):**

1. **Ultimate load test:**
   - Apply 1.5× limit load (or 1.0× ultimate load = g_limit)
   - Hold for 3 seconds without failure
   - **Required** for manned aircraft (FAA Part 23, EASA CS-23)

2. **Fatigue test:**
   - Apply load spectrum (representative flight profile)
   - Typically 2× design life in cycles
   - **Recommended** for high-utilization vehicles

3. **Material coupon tests:**
   - Verify assumed properties (E, strength, etc.)
   - Tension, compression, shear tests per ASTM standards
   - **Required** if using non-standard materials

4. **Ground vibration test (GVT):**
   - Measure natural frequencies, mode shapes
   - Validate flutter analysis
   - **Required** for high-speed aircraft

### 5.4 When Expert Consultation Needed

**Seek aerospace structural engineer if:**

- Designing manned aircraft (legal/safety requirements)
- FoS < 1.5 and unsure how to proceed
- Novel materials or construction methods
- Atypical load cases (crash, blast, bird strike)
- Certification or regulatory compliance required

---

## 6. Quality Assurance Checklist

### 6.1 Pre-Run Checklist

**Before running analysis, verify:**

- [ ] **Geometry**
  - [ ] Chord specified (required parameter)
  - [ ] Span reasonable (2 - 20 ft for typical wings)
  - [ ] Web positions ordered (web2 > web1)
  - [ ] Airfoil file exists and is readable

- [ ] **Materials**
  - [ ] Ply thickness in meters (not mm)
  - [ ] Moduli in Pascals (not GPa)
  - [ ] Allowables in Pascals (not MPa)
  - [ ] Densities in kg/m³
  - [ ] Face/web stacks have at least one ply

- [ ] **Flight Conditions**
  - [ ] V, ρ, CL are physical values
  - [ ] g_limit > 0 (typically 3 - 10)
  - [ ] Dynamic pressure q = 0.5ρV² is reasonable (50-500 Pa)

- [ ] **Output Settings**
  - [ ] Output filename is writable (no permission errors)
  - [ ] Sufficient disk space for Excel file (typically < 10 MB)

### 6.2 Post-Run Checklist

**After analysis completes, verify:**

- [ ] **Bredt Test**
  - [ ] Status: `bredt_solved`
  - [ ] Condition number < 1000
  - [ ] Accuracy < 1% error

- [ ] **Load Normalization**
  - [ ] Target lift = 0.5 × AUW × g_limit
  - [ ] Integrated lift matches target ± 0.1 N
  - [ ] Scale factor between 0.1 and 10

- [ ] **Geometry**
  - [ ] Neutral axis near zero (|z_NA| < 2 mm)
  - [ ] EI in reasonable range (1e-9 to 1e-6 N·m² for small wings)
  - [ ] Web heights positive (h_web > 0)

- [ ] **Results**
  - [ ] All FoS values positive and finite
  - [ ] Minimum FoS > 0 (no predicted failure if design is adequate)
  - [ ] Weight in expected range (5-10% of AUW typical)
  - [ ] Governing failure mode makes physical sense

- [ ] **Trends**
  - [ ] FoS decreases from tip to root (expected)
  - [ ] Root V approximately equals integrated lift
  - [ ] Plots show smooth curves (no jagged discontinuities)

### 6.3 Acceptance Criteria

**For a valid analysis, all of the following must be true:**

| Check | Criterion | Status |
|-------|-----------|--------|
| Bredt test | Passed (status = bredt_solved, error < 1%) | ✓ / ✗ |
| Load integration | ∫w dx = Target ± 0.1 N | ✓ / ✗ |
| No NaN/Inf | All results finite | ✓ / ✗ |
| Physical FoS | All FoS > 0 | ✓ / ✗ |
| Weight sanity | 100 - 500 g typical (check vs. similar designs) | ✓ / ✗ |
| Trend check | FoS decreases toward root | ✓ / ✗ |

**If any check fails, do not use results. Debug and re-run.**

---

## 7. Reporting and Documentation

### 7.1 Minimum Documentation for Design Review

**When presenting analysis results, include:**

1. **Input Summary**
   - Geometry (span, chord, web positions)
   - Material properties (fabric type, core grade, layup)
   - Flight conditions (V, ρ, CL, g_limit)

2. **Analysis Results**
   - Minimum FoS and location
   - Governing failure mode
   - Weight breakdown
   - Key plots (FoS_governing, root loads)

3. **Validation**
   - Bredt test result
   - Load integration check
   - Comparison to similar designs (if available)

4. **Limitations**
   - State assumptions (half-wing, no buckling, etc.)
   - Identify areas needing further analysis

5. **Recommendations**
   - Design changes if FoS < 1.5
   - Next steps (FEA, testing, etc.)

### 7.2 Archival

**Save for future reference:**
- Input files (airfoil coordinates, config parameters)
- Output files (Excel workbook, plots)
- Analysis log (terminal output, copy to text file)
- This validation checklist with filled values

**Version control:**
- Tag code version used (if using Git)
- Record date of analysis
- Note any custom modifications to code

---

## 8. Troubleshooting Guide

### 8.1 Results Don't Match Expectations

**Symptom:** FoS much higher or lower than expected.

**Possible causes:**

1. **Wrong units:**
   - Check: E in Pa (not GPa), t in m (not mm)
   - Fix: Convert all inputs to SI base units

2. **Wrong g_limit:**
   - Check: g_limit set to design value (not 1.0)
   - Fix: Update G_LIMIT in config

3. **Wrong airfoil:**
   - Check: Airfoil is correct file
   - Fix: Verify file contents, check thickness distribution

4. **Missing load:**
   - Check: Load integration equals target
   - Fix: Review load calculation, check for scale factor issues

### 8.2 Analysis Crashes or Hangs

**Symptom:** Program stops unexpectedly or runs forever.

**Possible causes:**

1. **Invalid geometry:**
   - Singular airfoil (zero thickness somewhere)
   - Self-intersecting airfoil
   - Fix: Check airfoil file, plot to visualize

2. **Ill-conditioned matrix:**
   - Bredt matrix singular (cells have zero area?)
   - Fix: Check web positions, cell areas in debug output

3. **Infinite loop:**
   - Rare, likely bug
   - Fix: Ctrl+C to interrupt, report issue

### 8.3 Plots Look Wrong

**Symptom:** Jagged curves, discontinuities, or unexpected shapes.

**Possible causes:**

1. **Too few stations:**
   - N_STATIONS = 11 → coarse, jagged
   - Fix: Increase to 51 or 101

2. **Numerical noise:**
   - Division by small numbers near tip (V → 0)
   - Fix: Normal, plots should still be interpretable

3. **Sign error:**
   - Loads negative when should be positive
   - Fix: Review sign convention, check if plot is upside-down

---

## Appendix: Validation Test Record Template

**Use this template to document validation for each configuration:**

```
=== WING ANALYSIS VALIDATION RECORD ===

Date: _______________
Analyst: _______________
Configuration: _______________

--- INPUT SUMMARY ---
Span: _______ ft
Chord: _______ m
Web positions: _______, _______
Core: ROHACELL ___ (t = ___ mm)
Fabric: _______ (n_plies per face = ___)
Rods: ___ mm × ___ mm
g_limit: _______

--- VALIDATION CHECKS ---
[ ] Bredt test passed (status = _______, cond = _______, error = _____%)
[ ] Load integration: Target = _____ N, Integrated = _____ N, Δ = _____ N
[ ] Neutral axis: z_NA = _____ mm (< 2 mm for symmetric?)
[ ] FoS all positive: min = _____, max = _____
[ ] Weight: Total = _____ g (reasonable?)
[ ] Trends: FoS decreases toward root? (yes/no)

--- RESULTS SUMMARY ---
Min FoS: _____ (location: x = _____ m)
Governing mode: _______________
Weight breakdown:
  Skins: _____ g
  Core: _____ g
  Webs: _____ g
  Rods: _____ g
  TOTAL: _____ g

--- ACCEPTANCE ---
[ ] All validation checks passed
[ ] Results physically reasonable
[ ] Documentation complete

Approved by: _______________ Date: _______________
```

---

*For usage instructions, see [USER_MANUAL.md](USER_MANUAL.md)*

*For theoretical background, see [TECHNICAL_REFERENCE.md](TECHNICAL_REFERENCE.md)*


