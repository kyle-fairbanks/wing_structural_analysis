# Wing Structural Analysis - Technical Reference

Engineering theory, equations, and material property reference.

**Version 1.0**

---

## Table of Contents

1. [Theoretical Basis](#1-theoretical-basis)
2. [Material Properties Reference](#2-material-properties-reference)
3. [Load Model](#3-load-model)
4. [Stress Calculations](#4-stress-calculations)
5. [Validation Cases](#5-validation-cases)
6. [References and Standards](#6-references-and-standards)

---

## 1. Theoretical Basis

### 1.1 Sandwich Beam Theory

#### 1.1.1 Fundamental Assumptions

**Sandwich structure:**
- Two thin, stiff face sheets (skins) separated by thick, lightweight core
- Faces carry bending (normal stresses)
- Core carries transverse shear and stabilizes faces against buckling
- Perfect bond between faces and core (no slip)

**Beam theory assumptions:**
- Plane sections remain plane (Euler-Bernoulli for slender beams)
- Small deflections (linear elastic)
- No transverse normal stress (σz = 0)
- Shear deformation in core allowed (core is soft in shear)

#### 1.1.2 Parallel-Axis Theorem

For composite sections with multiple materials, the total bending stiffness is:

```
EI_total = Σ E_i × (I_i + A_i × z_i²)
         = Σ E_i × I_i_about_NA
```

Where:
- `E_i` = elastic modulus of material i
- `I_i` = self moment of inertia of component i about its own centroid
- `A_i` = cross-sectional area of component i
- `z_i` = distance from neutral axis to centroid of component i

**Components in wing sandwich:**

1. **Top skin:**
   ```
   EI_skin_top = E_skin × (I_skin_self + A_skin × z_top²)
   ```

2. **Bottom skin:**
   ```
   EI_skin_bot = E_skin × (I_skin_self + A_skin × z_bot²)
   ```

3. **Core:**
   ```
   EI_core = E_core × (b × t_core³) / 12
   ```

4. **Rods (spar caps):**
   ```
   EI_rod_top = n_top × E_rod × (I_rod_self + A_rod × z_rod_top²)
   EI_rod_bot = n_bot × E_rod × (I_rod_self + A_rod × z_rod_bot²)
   ```

**Total:**
```
EI_total = EI_skin_top + EI_skin_bot + EI_core + EI_rod_top + EI_rod_bot
```

#### 1.1.3 Transformed Section Analysis

When materials have different moduli, use **transformed section** method:

1. **Choose reference modulus:** `E_ref` (typically use skin or fiber modulus)

2. **Compute modular ratios:**
   ```
   n_skin = E_skin / E_ref
   n_rod  = E_rod / E_ref
   n_core = E_core / E_ref
   ```

3. **Transformed areas:**
   ```
   A_skin_trans = n_skin × A_skin
   A_rod_trans  = n_rod × A_rod
   A_core_trans = n_core × A_core
   ```

4. **Neutral axis position:**
   ```
   z_NA = Σ(n_i × A_i × z_i) / Σ(n_i × A_i)
   ```

5. **Transformed moment of inertia:**
   ```
   I_trans = Σ n_i × (I_i_self + A_i × (z_i - z_NA)²)
   ```

6. **Stress recovery:**
   ```
   σ_i = (M × (z - z_NA) / I_trans) × (E_i / E_ref)
   ```

#### 1.1.4 Skin Placement

**Skin centroid locations (measured from core midplane):**

```
z_skin_top = +(t_core / 2) + (t_skin / 2)
z_skin_bot = -(t_core / 2) - (t_skin / 2)
```

**Rod placement (inside core, flush with inner skin surface):**

```
z_rod_desired = ±[(t_core / 2) - (rod_height / 2) - gap]
```

Where `gap` is clearance from inner skin surface into foam.

**Constraint:** If rod is taller than available core space, clamp to fit inside:
```
z_rod_actual = min(z_rod_desired, t_core/2 - margin)
```

---

### 1.2 Classical Laminate Theory (CLT)

#### 1.2.1 Ply Stiffness Matrix

For a single orthotropic ply with principal axes aligned with fiber direction:

```
       ⎡ σ₁ ⎤       ⎡ ε₁ ⎤
       ⎢ σ₂ ⎥ = [Q]⎢ ε₂ ⎥
       ⎣ τ₁₂⎦       ⎣γ₁₂⎦
```

**Reduced stiffness matrix Q (in principal axes):**

```
      ⎡ Q₁₁  Q₁₂   0  ⎤
[Q] = ⎢ Q₁₂  Q₂₂   0  ⎥
      ⎣  0    0   Q₆₆ ⎦
```

Where:
```
Q₁₁ = E₁ / (1 - ν₁₂ν₂₁)
Q₂₂ = E₂ / (1 - ν₁₂ν₂₁)
Q₁₂ = (ν₁₂ × E₂) / (1 - ν₁₂ν₂₁)
Q₆₆ = G₁₂

ν₂₁ = (E₂ / E₁) × ν₁₂  (reciprocity)
```

#### 1.2.2 Transformed Stiffness Matrix

For a ply at angle θ from the reference axis:

```
[Q̄] = [T]^T × [Q] × [T]
```

**Transformation (simplified for key terms):**

```
c = cos(θ), s = sin(θ)

Q̄₁₁ = Q₁₁c⁴ + 2(Q₁₂ + 2Q₆₆)c²s² + Q₂₂s⁴
Q̄₁₂ = (Q₁₁ + Q₂₂ - 4Q₆₆)c²s² + Q₁₂(c⁴ + s⁴)
Q̄₂₂ = Q₁₁s⁴ + 2(Q₁₂ + 2Q₆₆)c²s² + Q₂₂c⁴
Q̄₆₆ = (Q₁₁ + Q₂₂ - 2Q₁₂ - 2Q₆₆)c²s² + Q₆₆(c² - s²)²
```

**Special cases:**
- θ = 0°: Q̄ = Q (aligned)
- θ = 90°: Q̄₁₁ = Q₂₂, Q̄₂₂ = Q₁₁
- θ = ±45°: Maximum Q̄₆₆ (optimal for shear)

#### 1.2.3 Laminate ABD Matrix

For a laminate with N plies, each at height z_k to z_{k+1}:

**Extensional stiffness A:**
```
A_ij = Σ Q̄ᵢⱼ^(k) × (z_k - z_{k+1})
```

**Bending-extension coupling B:**
```
B_ij = (1/2) × Σ Q̄ᵢⱼ^(k) × (z_k² - z_{k+1}²)
```

**Bending stiffness D:**
```
D_ij = (1/3) × Σ Q̄ᵢⱼ^(k) × (z_k³ - z_{k+1}³)
```

**Constitutive equation:**
```
⎡ N ⎤   ⎡ A  B ⎤ ⎡ ε⁰ ⎤
⎢   ⎥ = ⎢     ⎥ ⎢    ⎥
⎣ M ⎦   ⎣ B  D ⎦ ⎣ κ  ⎦
```

Where:
- N = force resultants (N/m)
- M = moment resultants (N·m/m)
- ε⁰ = midplane strains
- κ = curvatures

#### 1.2.4 Effective Properties from CLT

**Effective bending modulus (per unit width):**

```
E_eff_bending = 12 × D₁₁ / t³
```

Where t is total laminate thickness.

**Membrane shear modulus:**

```
G_membrane = A₆₆ / t
```

**Implementation in code:**

See functions:
- `ply_stiffness()` (line 569)
- `transform_Q()` (line 580)
- `laminate_ABD_from_stack()` (line 604)
- `effective_bending_EI_from_D()` (line 633)
- `laminate_membrane_G()` (line 637)

---

### 1.3 Multi-Cell Torsion Theory (Bredt-Batho)

#### 1.3.1 Thin-Walled Section Assumptions

**Assumptions:**
- Wall thickness << cell dimensions
- Shear stress constant through wall thickness
- Shear flow `q = τ × t` is continuous around cells

**Variables:**
- `τ_j` = shear stress in wall j (Pa)
- `t_j` = wall thickness (m)
- `q_j` = shear flow in wall j (N/m)
- `A_c` = enclosed area of cell c (m²)
- `T` = applied torque (N·m)

#### 1.3.2 Governing Equations

**Torque equilibrium:**

For multi-cell section:
```
T = Σ 2 × A_c × q_c_net
```

Where `q_c_net` is the net shear flow around cell c.

**Twist compatibility:**

All cells must twist at the same rate θ'. For cells sharing walls:
```
θ'_c = (1 / 2A_c) × Σ (q_j × L_j) / (G_j × t_j)
```

Must be equal for all cells → system of equations.

#### 1.3.3 Matrix Formulation (Multi-Cell)

For `nc` cells and `nw` walls:

**Unknowns:** Shear flows `q_j` in each wall (length nw)

**Equations:**
1. **(nc - 1) compatibility equations:** Enforce equal twist rates between adjacent cells
2. **1 equilibrium equation:** Torque balance

**Matrix form:**
```
[M] {q} = {b}
```

Where:
- M is (nc × nw) matrix
- q is (nw × 1) vector of shear flows
- b is (nc × 1) vector with T in last row, zeros elsewhere

**Compatibility rows (i = 1 to nc-1):**
```
M[i, j] = (L_j / (G_j × t_j)) × (loop[i, j] - loop[0, j])
```

Where `loop[c, j] = 1` if wall j bounds cell c, else 0.

**Equilibrium row (last):**
```
M[nc, j] = 2 × Σ A_c_contrib_j
b[nc] = T
```

**Solution:**
```
{q} = [M]^(-1) {b}    (or least-squares if overdetermined)
```

**Stress recovery:**
```
τ_j = |q_j| / t_j
```

**Implementation:**

See functions:
- `assemble_bredt_matrix()` (line 760)
- `multicell_bredt_solver()` (line 801)

#### 1.3.4 Special Case: Single Cell

For single-cell closed section:

**Torque balance:**
```
T = 2 × A × q
```

**Twist rate:**
```
θ' = q × Σ(L_j / (G_j × t_j)) / (2A)
```

**Shear stress in wall j:**
```
τ_j = T / (2 × A × t_j)
```

**Analytic validation:**

Rectangular box (width B, height H, uniform thickness t, shear modulus G):
```
A = B × H
τ_analytic = T / (2 × B × H × t)
```

---

### 1.4 Shear Partition Model

#### 1.4.1 Physical Basis

In sandwich wing with webs and core, transverse shear force V is carried by **parallel load paths:**

1. **Webs:** Vertical walls carry shear via VQ/It mechanism
2. **Core:** Through-thickness shear across full width

**Key insight:** Webs and core act as **parallel springs** → load splits by stiffness ratio.

#### 1.4.2 Stiffness-Based Partition

**Shear stiffness (force per unit deflection):**
```
GA_web  = n_webs × G_web × t_web × h_web
GA_core = G_core × t_core × b_eff
```

Where:
- `n_webs` = number of webs (typically 2)
- `G_web` = web shear modulus (Pa, from CLT)
- `t_web` = web thickness (m)
- `h_web` = web height (m, rod-to-rod distance)
- `G_core` = core shear modulus (Pa, datasheet)
- `t_core` = core thickness (m)
- `b_eff` = effective core width (m, < chord to account for LE/TE ineffectiveness)

**Total stiffness:**
```
GA_total = GA_web + GA_core
```

**Force partition:**
```
V_web_total = V × (GA_web / GA_total)
V_core = V × (GA_core / GA_total)
```

**Check:** V_web_total + V_core = V ✓

#### 1.4.3 Effective Core Width

**Physical justification:**

Core near leading and trailing edges is less effective due to:
- Low thickness (airfoil tapers)
- Edge boundary effects
- Manufacturing quality (edges may have resin-rich regions)

**Model:**
```
b_eff = 0.85 × chord
```

(15% reduction accounts for LE/TE ineffectiveness)

**Sensitivity:** Results not highly sensitive to exact value (0.80 - 0.90 reasonable range).

#### 1.4.4 Advantages Over Traditional VQ/It-Only

**Traditional approach problems:**
1. Uses full V in VQ/It for webs → overestimates web stress
2. Uses full V for core → overestimates core stress
3. Double-counts load (V counted twice)

**GA partition approach:**
1. Physically consistent (parallel springs)
2. No circular logic
3. Each mechanism carries only its share of V
4. VQ/It used only for stress distribution within webs (with V_web, not full V)

**Implementation:**

See function `shear_partition_GA()` (line 438)

---

## 2. Material Properties Reference

### 2.1 Fabric: Hexcel FCIM255 (HiMax CGL4012)

**Type:** Carbon fiber biaxial fabric (±45°)

**Fiber:** T700S 12k tow, standard modulus

**Construction:**
- Biaxial (two directions: +45° and -45°)
- Stitched (not woven), minimal crimp
- Balanced (equal fiber in both directions)

**Properties (per datasheet):**

| Property | Value | Units | Notes |
|----------|-------|-------|-------|
| Areal weight | 200 | g/m² | Total for ±45° pair |
| Fabric thickness | 0.13 | mm | Dry fabric |
| Fiber direction | ±45° | deg | Relative to fabric edge |
| Fiber type | T700S 12k | - | Standard modulus |

**Usage in tool:**
- `ply_t_skin_mm = 0.13` (total fabric thickness)
- `face_stack = [(45, 1), (-45, 1)]` (one fabric layer)
- Areal weight used for weight estimation (200 gsm)

---

### 2.2 Fiber: T700S Carbon

**Manufacturer:** Toray

**Type:** Standard modulus, high strength carbon fiber

**Fiber Properties (from Toray datasheet):**

| Property | Value | Units |
|----------|-------|-------|
| Tensile modulus | 230 | GPa |
| Tensile strength | 4900 | MPa |
| Density | 1.80 | g/cm³ |
| Elongation at break | 2.1 | % |

**Composite Properties (UD, 60% fiber volume fraction, epoxy matrix):**

| Property | Symbol | Value | Units | Notes |
|----------|--------|-------|-------|-------|
| Longitudinal modulus | E₁ | 134 | GPa | 0° direction |
| Transverse modulus | E₂ | 10 | GPa | 90° direction |
| In-plane shear modulus | G₁₂ | 5.0 | GPa | Typical for carbon/epoxy |
| Poisson's ratio | ν₁₂ | 0.30 | - | Typical |
| Tensile strength (0°) | Xt | 2860 | MPa | From datasheet |
| Compressive strength (0°) | Xc | 1450 | MPa | From datasheet |
| In-plane shear strength | S | 136 | MPa | From datasheet |
| Density (cured UD) | ρ | 1500 | kg/m³ | At 60% Vf |

**Micromechanics check (Rule of Mixtures for E₁):**
```
E₁ = Vf × Ef + (1 - Vf) × Em
   = 0.60 × 230 + 0.40 × 3.5 = 139.4 GPa ≈ 134 GPa ✓
```

**Usage in tool:**
```python
E1_skin = 134e9    # Pa
E2_skin = 10e9     # Pa
G12_skin = 5.0e9   # Pa
nu12_skin = 0.30
Xt_skin = 2.86e9   # Pa
Xc_skin = 1.45e9   # Pa
tau_skin_allow = 136e6  # Pa
```

---

### 2.3 Core: ROHACELL 31 IG-F

**Manufacturer:** Evonik

**Type:** Closed-cell rigid polymethacrylimide (PMI) foam

**Grade:** 31 IG-F (31 nominal density, IG = Integral Foam, F = Fine cell structure)

**Properties (from Evonik datasheet, 23°C, dry):**

| Property | Symbol | Value | Tolerance | Units |
|----------|--------|-------|-----------|-------|
| Density | ρ | 32 | ±7 | kg/m³ |
| Compressive modulus | E_c | 36 | - | MPa |
| Compressive strength | σ_c | 0.4 | - | MPa |
| Tensile modulus | E_t | 36 | - | MPa |
| Tensile strength | σ_t | 0.7 | - | MPa |
| Shear modulus | G | 13 | - | MPa |
| Shear strength | τ | 0.4 | - | MPa |
| Elongation at break | ε | 4 | - | % |

**Temperature effects:**
- Tg (glass transition) > 180°C
- Service temperature: −200 to +180°C
- Properties stable at typical flight temperatures (−50 to +70°C)

**Moisture effects:**
- Closed-cell → minimal water absorption
- < 2% weight gain at 100% RH

**Usage in tool:**
```python
t_core = 0.003          # m (3.0 mm)
rho_core = 32.0         # kg/m³
E_core = 36e6           # Pa (36 MPa)
G_core = 13e6           # Pa (13 MPa)
tau_core_allow = 0.4e6  # Pa (0.4 MPa)
```

**Design considerations:**
- Shear strength (0.4 MPa) often limiting
- For higher loads, consider ROHACELL 51 (0.7 MPa), 71 (1.2 MPa), or 110 (2.3 MPa)
- Denser grades have proportionally higher E, G, and strengths

---

### 2.4 Rods / Spar Caps: T700S UD Carbon/Epoxy

**Type:** Pultruded or filament-wound unidirectional carbon/epoxy

**Fiber:** T700S, aligned in spanwise direction (0° orientation)

**Properties (from representative UD pultruded rod datasheet):**

| Property | Symbol | Value | Units | Notes |
|----------|--------|-------|-------|-------|
| Tensile modulus | E | 138 | GPa | 0° direction, ~60% Vf |
| Tensile strength | X_t | 1720 | MPa | Ultimate |
| Flexural strength | σ_f | 1830 | MPa | Used as proxy for compression |
| Poisson's ratio | ν | 0.30 | - | Typical |
| Density | ρ | 1500 | kg/m³ | 1.5 g/cm³ |

**Usage in tool:**
```python
E_rod = 138e9           # Pa
nu_rod = 0.30
Xt_rod = 1.72e9         # Pa (tension)
Xc_rod = 1.83e9         # Pa (flexural as compression proxy)
rho_rod = 1500.0        # kg/m³
```

**Geometry:**
- Rectangular cross-section (typical for small wings)
- Default: 3 mm × 3 mm
- Can be changed via `ROD_WIDTH_MM`, `ROD_HEIGHT_MM` in config

**Placement:**
- Inside foam core
- Flush with inner skin surface (gap = 0 mm clearance)
- 4 total: 2 top, 2 bottom

---

### 2.5 Comparison Table

| Material | E (GPa) | G (GPa) | ρ (kg/m³) | σ_allow (MPa) | τ_allow (MPa) | Cost |
|----------|---------|---------|-----------|---------------|---------------|------|
| T700S UD (skin) | 134 | 5.0 | 1500 | 1450 - 2860 | 136 | $$ |
| T700S UD (rod) | 138 | - | 1500 | 1720 - 1830 | - | $$ |
| ROHACELL 31 | 0.036 | 0.013 | 32 | 0.4 - 0.7 | 0.4 | $ |
| ROHACELL 71 | 0.105 | 0.042 | 75 | 1.0 - 1.6 | 1.2 | $$ |
| ROHACELL 110 | 0.180 | 0.070 | 110 | 1.7 - 2.9 | 2.3 | $$$ |

**Ratios (dimensionless):**

| Ratio | Value | Implication |
|-------|-------|-------------|
| E_fiber / E_core | 3700× | Fiber dominates bending stiffness |
| G_fiber / G_core | 385× | Fiber dominates torsional stiffness |
| ρ_fiber / ρ_core | 47× | Core is weight-efficient for stiffness |
| Xt_fiber / Xc_core | 7150× | Fiber carries bending, core carries shear |

---

## 3. Load Model

### 3.1 Half-Wing Cantilever Model

**Geometry:**
- Half-wing from root (fuselage attachment) to tip
- Span = distance from root to tip (user input, converted to meters)
- Uniform chord (no taper in current version)

**Boundary conditions:**
- Root: Fully fixed (all DOF constrained)
- Tip: Free

**Coordinate system:**
- x: spanwise (root = 0, tip = span/2)
- y: chordwise (LE = 0, TE = chord)
- z: vertical (midplane = 0, up = positive)

### 3.2 Lift Distribution

#### 3.2.1 Spanwise Load Shape

**Goal:** Realistic lift distribution that:
1. Peaks slightly inboard (not at root, not at mid-span)
2. Has broad plateau inboard
3. Decays smoothly to zero at tip
4. Is symmetric about root for full wing

**Analytic model:**

```python
eta = x / (b/2)   # normalized span position (0 at root, 1 at tip)
base = 1 - eta²   # elliptic base (zero at tip)
shape_raw = base^n × (1 + a × eta² × base)
```

**Parameters (tuned):**
- `n = 0.20` (very gentle falloff, creates plateau)
- `a = 0.80` (mid-span "shoulder" bump)

**Normalization:**

```python
shape = shape_raw / max(shape_raw)  # peak value = 1
```

**Result:**
- Root: w(0) ≈ 0.90 (90% of peak)
- Peak: x ≈ 0.15 to 0.25 (slightly inboard)
- Tip: w(b/2) = 0 (enforced)

#### 3.2.2 Load Magnitude Normalization

**Target total lift:**

```
L_total = AUW × g_limit
```

For half-wing:
```
L_half = 0.5 × AUW × g_limit
```

**Distributed load:**

```python
w(x) = w_shape(x) × scale_factor
```

Where `scale_factor` is chosen such that:

```
∫[0 to b/2] w(x) dx = L_half
```

**Implementation:**

```python
L_half_integrated = trapezoid(w_initial, xspan)
scale = L_half_target / L_half_integrated
w_final = w_initial × scale
```

Verification printed during analysis:
```
[LOAD-NORM] Target half-wing lift @ g_limit: 155.324 N
[LOAD-NORM] Integrated half-wing lift (after scale): 155.324 N
```

#### 3.2.3 Aerodynamic vs. Weight-Based Loading

**Two paths to determine lift magnitude:**

**Path 1: Aerodynamic (if V, ρ, CL provided)**
```
q = 0.5 × ρ × V²
L_full_1g = q × CL × chord × span_full
L_half_target = 0.5 × L_full_1g × g_limit
```

**Path 2: Wing loading (if aerodynamic data unavailable)**
```
WL = wing_loading_lbft2 × 47.88 N/m² per lb/ft²
A_full = span_full × chord
L_full_1g = WL × A_full
L_half_target = 0.5 × L_full_1g × g_limit
```

**Both paths produce identical wing loading at design condition.**

---

### 3.3 Sweep Projection

#### 3.3.1 Load Components

For sweep angle Λ, the lift force projects into:

**Flapwise (vertical plane perpendicular to spar):**
```
L_flap = L_total × cos(Λ)
```

**In-plane (horizontal, parallel to spar):**
```
L_inplane = L_total × sin(Λ)
```

**Current model:**
- Flapwise loads drive bending V, M
- In-plane loads currently neglected (small effect for Λ < 30°, conservative)

#### 3.3.2 Bending and Shear Integration

**Shear force V(x):**

```
V(x) = ∫[x to b/2] w_flap(ξ) dξ
```

Integrated tip-to-root using trapezoidal rule.

**Bending moment M(x):**

```
M(x) = ∫[x to b/2] V(ξ) dξ
```

Also integrated tip-to-root.

**Sign convention:**
- Positive V: upward shear force (tends to shear section upward)
- Positive M: causes tension on top surface, compression on bottom (sagging)

#### 3.3.3 Torsional Loading

**Mechanism:**

Lift acts at chordwise location `x_lift` (typically near quarter-chord).  
Structural neutral axis at `x_NA` (depends on section properties).

**Moment arm:**
```
d = x_lift - x_NA
```

**Distributed torque:**
```
m_t(x) = w_flap(x) × d
```

Projected by sweep:
```
m_t_proj(x) = m_t(x) × cos(Λ)
```

**Integrated torsional moment T(x):**

```
T(x) = ∫[x to b/2] m_t_proj(ξ) dξ
```

**Sign convention:**
- Positive T: nose-up twist (trailing edge twists upward relative to leading edge)

---

## 4. Stress Calculations

### 4.1 Bending Stresses

#### 4.1.1 Physical Approach (Curvature-Based)

**Beam curvature:**
```
κ = M / EI_total
```

**Strain at distance z from neutral axis:**
```
ε(z) = κ × z
```

**Stress in material with modulus E:**
```
σ(z) = E × ε(z) = E × κ × z
```

#### 4.1.2 Skin Bending Stress

**CLT-derived effective bending modulus:**

```
E_eff_bending = 12 × D₁₁ / t_skin³
```

Where D₁₁ is from ABD matrix.

**Skin stress:**

```
σ_skin_top = E_eff_bending × κ × z_skin_top
σ_skin_bot = E_eff_bending × κ × z_skin_bot
```

Where:
```
z_skin_top = +(t_core/2 + t_skin/2)
z_skin_bot = -(t_core/2 + t_skin/2)
```

#### 4.1.3 Rod Bending Stress

**Rod stress (isotropic material, simple):**

```
σ_rod_top = E_rod × κ × z_rod_top
σ_rod_bot = E_rod × κ × z_rod_bot
```

Where z_rod positions are determined by placement geometry (see Section 1.1.4).

#### 4.1.4 Sign Convention

- Positive σ: tension
- Negative σ: compression
- Positive M (sagging): tension on top, compression on bottom
- Negative M (hogging): compression on top, tension on bottom

---

### 4.2 Web Shear Stress

#### 4.2.1 VQ/It Formula (Classical Beam Shear)

For vertical webs carrying transverse shear:

```
τ_web = V_web × Q / (I_trans × t_web)
```

**Key:** Use `V_web` (web portion from GA partition), not total V.

**First moment of area Q:**

For vertical web cut at neutral axis:

```
Q = ∫[forward of web] y × dA_trans
```

Where:
- y = distance from neutral axis
- dA_trans = transformed area element
- Integration over all material forward of the web

**Transformed moment of inertia:**

```
I_trans = ∫ (n_i × y²) dA_i
```

Integrated over entire cross-section about neutral axis.

#### 4.2.2 Implementation Details

**Q calculation (continuous integration):**

```python
# Skin contribution (top)
Q_skin_top = ∫[x=0 to x_web] A_strip × (y_upper - z_NA) dx

# Skin contribution (bottom)
Q_skin_bot = ∫[x=0 to x_web] A_strip × (y_lower - z_NA) dx

# Rod contribution (discrete)
Q_rods = Σ A_rod_trans × (z_rod - z_NA)  for rods forward of web

# Total
Q_total = Q_skin_top + Q_skin_bot + Q_rods
```

**I_trans calculation (global section property):**

Computed once using transformed section method (Section 1.1.3).

**Consistent usage:**
- Q_1 for web 1 (uses material forward of web 1)
- Q_2 for web 2 (uses material forward of web 2)
- I_trans is same for both (global property)

---

### 4.3 Core Shear Stress

#### 4.3.1 Through-Thickness Shear Model

**Mechanism:**

Core carries shear via through-thickness distortion (like shear in homogeneous beam).

**Stress:**

```
τ_core = V_core / A_core_shear
```

Where:
```
A_core_shear = b_eff × t_core
```

And `b_eff = 0.85 × chord` (effective width, accounting for LE/TE regions).

#### 4.3.2 Physical Justification

**Why not VQ/It for core?**

Core shear is **not** beam shear (Q/I mechanism) because:
1. Core is thick (t/h ≈ 0.5, not thin)
2. Core is very soft (E_core << E_skin)
3. Core strain is dominated by shear deformation, not bending-induced shear

**Correct model:** Direct shear stress from shear force:
```
τ = V / A
```

This is analogous to Timoshenko beam shear correction.

---

### 4.4 Torsional Stresses

#### 4.4.1 Bredt Shear Flow Solution

**Process:**

1. Solve Bredt system for shear flows q_j in each wall
2. Convert to stress: τ_j = |q_j| / t_j

**Skin torsion stress:**

Maximum of perimeter wall stresses (top, mid-chord, bottom skins):
```
τ_skin_tor = max(τ_top_skin, τ_mid_skin, τ_bot_skin)
```

**Web torsion stress:**

From internal web walls:
```
τ_web1_tor = |q_web1| / t_web
τ_web2_tor = |q_web2| / t_web
```

#### 4.4.2 Total Web Stress

**Superposition:**

Webs experience both VQ/It shear (from transverse V) and Bredt shear (from torsion T):

```
τ_web_VQIt = V_web × Q / (I_trans × t_web)
τ_web_Bredt = |q_web| / t_web
```

**Conservative combination (magnitudes sum):**

```
τ_web_total = τ_web_VQIt + τ_web_Bredt
```

**In code:**

- `tau_web1_arr` = VQ/It component (from shear force)
- `tau_web1` = Bredt component (from torsion)
- These are stored separately for diagnostics
- FoS uses worst-case of each component

---

### 4.5 Failure Criteria

#### 4.5.1 Normal Stress (Tension/Compression)

**Simple max-stress criterion:**

```
σ_applied ≥ 0 (tension):   FoS = Xt / σ_applied
σ_applied < 0 (compression): FoS = Xc / |σ_applied|
```

**Applied to:**
- Skin top/bottom
- Rod top/bottom

#### 4.5.2 Shear Stress (Simple)

**Max shear stress criterion:**

```
FoS = τ_allow / |τ_applied|
```

**Applied to:**
- Web shear (VQ/It + Bredt)
- Core shear
- Skin membrane shear (if used)

#### 4.5.3 Tsai-Wu Composite Criterion (2D)

**General form:**

```
F₁σₓ + F₂σᵧ + F₁₁σₓ² + F₂₂σᵧ² + F₆₆τₓᵧ² + 2F₁₂σₓσᵧ < 1
```

**Coefficients:**

```
F₁ = 1/Xt - 1/Xc
F₂ = 1/Yt - 1/Yc
F₁₁ = 1/(Xt × Xc)
F₂₂ = 1/(Yt × Yc)
F₆₆ = 1/S²
F₁₂ = -0.5 × √(F₁₁ × F₂₂)  (typical assumption)
```

**Factor of Safety:**

If Tsai-Wu index `TW > 0`, then:
```
FoS_TW = 1 / √TW
```

**Applied to:**
- Skin torsion (combines normal stress from bending + shear stress from torsion)

**Conservative approach in code:**

For skins, torsion stress is Bredt wall stress (not laminate in-plane stress), so:
```
TW_index = tsai_wu(σₓ = σ_skin_bending, σᵧ = 0, τₓᵧ = 0, ...)
```

Setting τₓᵧ = 0 is conservative (ignores beneficial shear capacity).

---

## 5. Validation Cases

### 5.1 Built-In Validation

#### 5.1.1 Single-Cell Bredt Test

**Function:** `validate_bredt_single_cell()` (lines 2502-2518)

**Geometry:**
- Rectangular box: B = 0.2 m, H = 0.1 m
- Uniform thickness: t = 1 mm
- Uniform shear modulus: G = 1 GPa
- Applied torque: T = 10 N·m

**Analytic solution:**

```
A_cell = B × H = 0.02 m²
τ_analytic = T / (2 × A_cell × t)
          = 10 / (2 × 0.02 × 0.001)
          = 250,000 Pa
```

**Numerical solution:**

Bredt solver assembles matrix and solves for q in each wall.

**Expected result:**
```
τ_solver = [250000, 250000, 250000, 250000] Pa
```

(all four walls have same stress for rectangular section with uniform properties)

**Acceptance criteria:**
- Status: `bredt_solved`
- Condition number < 1000
- |τ_solver - τ_analytic| / τ_analytic < 1%

**Printed at startup:**
```
Bredt single-cell test status: bredt_solved cond: 124.5
  tau analytic: 250000.0
  tau solver: [250000. 250000. 250000. 250000.]
```

#### 5.1.2 Load Integration Check

**Test:** Verify spanwise lift distribution integrates to target value.

**Printed during analysis:**
```
[LOAD-NORM] Target half-wing lift @ g_limit: 155.324 N
[LOAD-NORM] Integrated half-wing lift (after scale): 155.324 N
[LOAD-NORM] Applied scale factor: 1.234567
```

**Acceptance:**
- Integrated value matches target to < 0.1 N (numerical integration tolerance)

---

### 5.2 Material Property Consistency Checks

#### 5.2.1 CLT-Derived vs. Datasheet G

**Test:** Compare G from CLT A₆₆ with datasheet value.

**CLT:**
```
G_membrane = A₆₆ / t_laminate
```

For ±45° biax (one layer each direction):
```
G_membrane ≈ 30-40 GPa  (high due to fiber orientation)
```

**Datasheet:**
```
G₁₂ = 5.0 GPa  (for 0° UD)
```

**Note:** These are **not** the same property!
- G_membrane is **in-plane** shear modulus for the laminate (dominated by fiber tension/compression at ±45°)
- G₁₂ is **ply** shear modulus (matrix-dominated, between fibers)

**For ±45° laminates:**
```
G_membrane ≈ (E₁ + E₂) / 4 + G₁₂ / 2 ≈ 38 GPa
```

This is **correct** and much higher than G₁₂ due to Fiber rotation effect.

**Validation:** No inconsistency; different properties for different purposes.

#### 5.2.2 Neutral Axis Check

**Test:** For symmetric sandwich, z_NA should be near zero.

**Expected:**
```
|z_NA| < 1 mm  (small offset acceptable due to discrete rods)
```

**Printed:**
```
Neutral axis position: z_NA = 0.0123 mm (from geometric midplane)
```

**Interpretation:**
- z_NA ≈ 0: Symmetric section ✓
- z_NA >> 0: Asymmetry (more material on top)
- z_NA << 0: Asymmetry (more material on bottom)

---

### 5.3 External Validation Recommendations

#### 5.3.1 Simple Beam Comparison

**Test case:**

Cantilever beam with uniform load:
- Length L = 0.5 m
- Uniform load w = 1000 N/m
- EI = known value (e.g., 100 N·m²)

**Analytic solution:**

```
V(x) = w × (L - x)
M(x) = w × (L - x)² / 2
δ_tip = w × L⁴ / (8 × EI)
```

**Test:**
- Run tool with this configuration
- Compare root V, M to analytic
- Compare tip deflection (if code extended to include deflections)

#### 5.3.2 Torsion FEA Benchmark

**Test case:**

Model simple multi-cell box beam in FEA (ANSYS, Abaqus, etc.):
- Shell elements for walls
- Apply torque at one end
- Extract shear stress in walls

**Compare:**
- Wall shear stresses τ_j from FEA vs. Bredt solver
- Expect agreement within 5-10% (FEA includes 3D effects, stress concentrations)

#### 5.3.3 CLT Software Cross-Check

**Test case:**

Use dedicated CLT software (e.g., ESAComp, LAMINATOR):
- Input same ply stack: [(45, 1), (-45, 1)]
- Extract ABD matrix
- Compare D₁₁, A₆₆ values

**Expected:** Agreement to within 1% (numerical precision).

---

### 5.4 Sensitivity Studies

#### 5.4.1 Mesh Convergence (Spanwise Stations)

**Test:**

Run same case with varying `N_STATIONS`:
```
N = [21, 51, 101, 201]
```

**Monitor:**
- Root M, V (should converge)
- Min FoS (should converge)

**Expected:**
- N ≥ 51: Results within 1%
- N = 101: Baseline (good balance of speed/accuracy)

#### 5.4.2 Core Effective Width Sensitivity

**Test:**

Vary `b_eff` factor:
```
b_eff = [0.70, 0.80, 0.85, 0.90, 1.00] × chord
```

**Monitor:**
- Core shear FoS (inversely proportional to b_eff)
- Overall min FoS (may be sensitive if core shear governs)

**Expected:**
- 15-20% variation in core shear FoS
- If other modes govern, minimal sensitivity

---

## 6. References and Standards

### 6.1 Composite Mechanics References

1. **Daniel, I.M., Ishai, O.** (2006). *Engineering Mechanics of Composite Materials*, 2nd ed. Oxford University Press.
   - Chapter 6: Classical Laminate Theory

2. **Jones, R.M.** (1999). *Mechanics of Composite Materials*, 2nd ed. Taylor & Francis.
   - Comprehensive treatment of CLT and failure theories

3. **Tsai, S.W., Hahn, H.T.** (1980). *Introduction to Composite Materials*. Technomic Publishing.
   - Tsai-Wu failure criterion derivation

### 6.2 Structural Analysis References

4. **Megson, T.H.G.** (2013). *Aircraft Structures for Engineering Students*, 5th ed. Butterworth-Heinemann.
   - Chapter 19: Thin-Walled Beams (Bredt-Batho theory)

5. **Bruhn, E.F.** (1973). *Analysis and Design of Flight Vehicle Structures*. Tri-State Offset Company.
   - Classic reference for aircraft structures

6. **Niu, M.C.Y.** (1988). *Airframe Structural Design*. Conmilit Press.
   - Practical design methods for composite wings

### 6.3 Material Property Sources

7. **Toray Carbon Fibers America, Inc.** T700S Data Sheet.
   - https://www.toray.com (Technical Data Sheets section)

8. **Evonik Industries AG.** ROHACELL IG/IG-F Product Information.
   - https://www.rohacell.com (ROHACELL PMI Foams section)

9. **HiMax Composites.** FCIM255 (CGL4012) Fabric Specifications.
   - Manufacturer technical data sheets

### 6.4 Relevant Standards

10. **ASTM D3039** - Standard Test Method for Tensile Properties of Polymer Matrix Composite Materials

11. **ASTM D3410** - Standard Test Method for Compressive Properties of Polymer Matrix Composite Materials

12. **ASTM D5379** - Standard Test Method for Shear Properties of Composite Materials (Iosipescu Shear)

13. **MIL-HDBK-17-1F** - Composite Materials Handbook, Volume 1: Polymer Matrix Composites Guidelines for Characterization of Structural Materials
    - Now CMH-17 (managed by SAE)

14. **FAA Advisory Circular 20-107B** - Composite Aircraft Structure
    - Guidelines for certification of composite aircraft

### 6.5 Software Validation References

15. **ANSYS Composite PrepPost (ACP)** - FEA validation for CLT
    - User's guide and verification manual

16. **ESAComp** - Commercial CLT and structural analysis software
    - Benchmark database for validation

---

## Appendix A: Notation and Symbols

### A.1 Geometry

| Symbol | Description | Units |
|--------|-------------|-------|
| b | Wing span (full or half, context-dependent) | m |
| c | Chord length | m |
| x | Spanwise coordinate (root = 0) | m |
| y | Chordwise coordinate (LE = 0) | m |
| z | Vertical coordinate (midplane = 0) | m |
| t_skin | Skin thickness (one face) | m |
| t_core | Core thickness | m |
| t_web | Web thickness | m |
| A | Area | m² |
| I | Moment of inertia | m⁴ |

### A.2 Material Properties

| Symbol | Description | Units |
|--------|-------------|-------|
| E | Elastic modulus | Pa |
| G | Shear modulus | Pa |
| ν | Poisson's ratio | - |
| ρ | Density | kg/m³ |
| Xt | Tensile strength | Pa |
| Xc | Compressive strength | Pa |
| τ_allow | Shear strength | Pa |

### A.3 Loads and Stresses

| Symbol | Description | Units |
|--------|-------------|-------|
| V | Shear force | N |
| M | Bending moment | N·m |
| T | Torsional moment | N·m |
| w | Distributed load | N/m |
| σ | Normal stress | Pa |
| τ | Shear stress | Pa |
| q | Shear flow | N/m |
| κ | Curvature | 1/m |

### A.4 CLT Notation

| Symbol | Description | Units |
|--------|-------------|-------|
| [Q] | Reduced stiffness matrix (ply) | Pa |
| [Q̄] | Transformed stiffness matrix | Pa |
| [A] | Extensional stiffness matrix | N/m |
| [B] | Bending-extension coupling matrix | N |
| [D] | Bending stiffness matrix | N·m |
| N | Force resultant | N/m |
| M | Moment resultant | N·m/m |
| ε⁰ | Midplane strain | - |
| κ | Curvature | 1/m |

---

*For usage instructions, see [USER_MANUAL.md](USER_MANUAL.md)*

*For validation procedures, see [VALIDATION.md](VALIDATION.md)*




