from __future__ import annotations
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wing_structural_analysis_merged_final_clt.py

Revised to include Classical Laminate Theory for skins and webs,
composite Tsai-Wu checks, and CLT-derived shear moduli for Bredt.
Geometry handling left unchanged.
"""

import argparse
import math
import os
from dataclasses import dataclass, asdict
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# shapely required
try:
    from shapely.geometry import Polygon, LineString, LinearRing, MultiLineString
    from shapely.ops import unary_union
    from shapely import ops
    import shapely.validation as shval
except Exception as e:
    raise ImportError("shapely is required. Install with: pip install shapely") from e

# scipy optional for smoothing
try:
    from scipy.interpolate import splprep, splev, splrep
    from scipy.signal import savgol_filter
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

# ---------------------- MATERIALS / INPUTS ----------------------

@dataclass
class Materials:
    """
    Material properties for wing structural analysis.
    
    Fabric: FCIM255 (HiMax CGL4012) - ±45° biax carbon fabric
    Fiber: T700S 12k standard modulus carbon fiber
    Core: ROHACELL 31 IG-F foam
    Rods: T700S UD carbon/epoxy
    """
    # Face sheet: list of (theta_deg, n_plies) top-to-bottom for a face (will be mirrored)
    face_stack: Optional[List[Tuple[float, int]]] = None   # e.g. [(45, 1), (-45, 1)]
    # web stack (list of (theta_deg, n_plies))
    web_stack: Optional[List[Tuple[float, int]]] = None

    # Per-ply thickness (mm) - FCIM255 fabric properties
    ply_t_skin_mm: float = 0.13    # FCIM255 total fabric thickness (200 gsm biax ±45°)
    ply_t_web_mm: float = 0.13     # FCIM255 total fabric thickness (200 gsm biax ±45°)

    # Core (Rohacell 31 IG-F, SI units per datasheet)
    t_core: float = 0.003          # m (3.0 mm)
    rho_core: float = 32.0         # kg/m³ (32 ± 7 per datasheet)
    E_core: float = 36e6           # Pa, 36 MPa (tensile modulus per datasheet)
    G_core: float = 13e6           # Pa, 13 MPa (shear modulus per datasheet)
    tau_core_allow: float = 0.4e6  # Pa, 0.4 MPa (shear strength per datasheet)

    # Rods / spar caps (0° UD carbon/epoxy per provided datasheet)
    E_rod: float = 138e9           # Pa, 138 GPa (20.0 msi, tensile modulus)
    nu_rod: float = 0.30           # typical for carbon/epoxy
    Xt_rod: float = 1.72e9         # Pa, 1.72 GPa (250 ksi, tensile strength)
    Xc_rod: float = 1.83e9         # Pa, 1.83 GPa (265 ksi, flexural strength as proxy)
    rho_rod: float = 1500.0        # kg/m³ (1.5 g/cm³, 0.54 lbs/in³)

    # Representative ply (T700S fiber-dominated) elastic properties
    E1_fiber: float = 134e9        # Pa, 134 GPa (T700S UD composite, 60% Vf)
    E2_fiber: float = 10e9         # Pa, 10 GPa (transverse, typical)
    G12_fiber: float = 5.0e9       # Pa, 5 GPa (typical for carbon/epoxy)
    nu12_fiber: float = 0.30       # typical

    # Skin laminate properties (FCIM255 ±45° biax with T700S fiber)
    E1_skin: float = 134e9         # Pa, 134 GPa (T700S at 60% Vf)
    E2_skin: float = 10e9          # Pa, 10 GPa (transverse)
    G12_skin: float = 5.0e9        # Pa, 5 GPa (in-plane shear modulus)
    nu12_skin: float = 0.30        # typical

    # Web laminate properties (FCIM255 ±45° biax with T700S fiber)
    E1_web: float = 134e9          # Pa, 134 GPa (T700S at 60% Vf)
    E2_web: float = 10e9           # Pa, 10 GPa (transverse)
    G12_web: float = 5.0e9         # Pa, 5 GPa (in-plane shear modulus)
    nu12_web: float = 0.30         # typical

    # Allowables (T700S composite at 60% Vf per datasheet)
    Xt_skin: float = 2.86e9        # Pa, 2,860 MPa (tensile, T700S datasheet)
    Xc_skin: float = 1.45e9        # Pa, 1,450 MPa (compressive, T700S datasheet)
    tau_skin_allow: float = 136e6  # Pa, 136 MPa (in-plane shear, T700S datasheet)

    # Web shear allowable (T700S composite at 60% Vf per datasheet)
    tau_web_allow: float = 136e6   # Pa, 136 MPa (in-plane shear, T700S datasheet)

    # legacy/fallback
    G_web_45: float = 5.0e9        # Pa, 5 GPa

    def __post_init__(self):
        if self.face_stack is None:
            self.face_stack = [(45.0, 1), (-45.0, 1)]
        if self.web_stack is None:
            self.web_stack = [(45.0, 1), (-45.0, 1)]

    @property
    def t_ply_skin(self) -> float:
        return max(1e-6, float(self.ply_t_skin_mm) * 1e-3)

    @property
    def t_ply_web(self) -> float:
        return max(1e-6, float(self.ply_t_web_mm) * 1e-3)

    @property
    def t_skin_total(self) -> float:
        t = 0.0
        for theta, n in self.face_stack:
            t += int(n) * self.t_ply_skin
        return t

    @property
    def t_web_total(self) -> float:
        t = 0.0
        for theta, n in self.web_stack:
            t += int(n) * self.t_ply_web
        return t

@dataclass
class Inputs:
    airfoil_file: str = "Modified_NACA_6_Series.txt"
    span_ft: float = 4.53
    sweep_deg: float = 0.0
    chord_m: Optional[float] = None          # MUST be provided or raise
    web1_p: float = 0.2
    web2_p: float = 0.65
    rod_width_m: float = 3.0e-3
    rod_height_m: float = 3.0e-3
    rod_gap_mm: float = 0
    wing_loading_lbft2: float = 3.96
    AUW_lb: float = 6.33
    g_limit: float = 5.0
    n_stations: int = 101
    out_excel: str = "wing_results_final.xlsx"
    make_plots: bool = True
    lift_offset_c: float = 0.25
    x_sc_c: Optional[float] = None
    cm_ac: float = -0.06
    q_pa: Optional[float] = None
    cl: float = 1.34
    rho: Optional[float] = 1.0573
    V_ms: Optional[float] = 15.87

# ---------------------- AIRFOIL IO + SMOOTHING ----------------------

def load_airfoil_points(path: str) -> np.ndarray:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Airfoil file not found: {path}")
    try:
        data = np.loadtxt(p, usecols=(0, 1))
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data[:, :2].astype(float)
    except Exception:
        df = pd.read_csv(p, header=None)
        if df.shape[1] >= 2:
            return df.iloc[:, :2].to_numpy(dtype=float)
        else:
            raise ValueError("Airfoil file could not be parsed as two columns.")

def normalize_and_order_loop(points: np.ndarray) -> np.ndarray:
    x = points[:, 0].astype(float).copy()
    y = points[:, 1].astype(float).copy()
    i_le = int(np.argmin(x))
    i_te = int(np.argmax(x))
    chord = x[i_te] - x[i_le]
    if chord <= 0:
        raise ValueError("Invalid chord derived from airfoil points.")
    x = (x - x[i_le]) / chord
    order = np.concatenate((np.arange(i_le, len(x)), np.arange(0, i_le)))
    x = x[order]; y = y[order]
    imax = int(np.argmax(x))
    upper = np.column_stack([x[:imax+1], y[:imax+1]])
    lower = np.column_stack([x[imax:], y[imax:]])

    def branch_sort(branch):
        xb = branch[:, 0]; yb = branch[:, 1]
        idx = np.argsort(xb)
        xb = xb[idx]; yb = yb[idx]
        ux, inv, cnt = np.unique(xb, return_inverse=True, return_counts=True)
        uy = np.zeros_like(ux)
        np.add.at(uy, inv, yb)
        uy = uy / np.maximum(cnt, 1)
        return np.column_stack([ux, uy])

    upper = branch_sort(upper)
    lower = branch_sort(lower)
    if lower[0, 0] > lower[-1, 0]:
        lower = lower[::-1]
    loop = np.vstack([upper, lower[1:-1][::-1], upper[0:1]])
    return loop

def build_smooth_envelopes(loop_xy: np.ndarray, ngrid: int = 301, smooth: float = 0.0):
    xs = loop_xy[:, 0]; ys = loop_xy[:, 1]
    i_te = int(np.argmax(xs))
    xu = xs[:i_te+1]; yu = ys[:i_te+1]
    xl = xs[i_te:]; yl = ys[i_te:]

    def unique_xy(xb, yb):
        idx = np.argsort(xb)
        xb = xb[idx]; yb = yb[idx]
        ux, inv, cnt = np.unique(xb, return_inverse=True, return_counts=True)
        uy = np.zeros_like(ux)
        np.add.at(uy, inv, yb)
        uy /= np.maximum(cnt, 1)
        return ux, uy

    xu, yu = unique_xy(xu, yu)
    xl, yl = unique_xy(xl, yl)

    xu[0], xu[-1] = 0.0, 1.0
    xl[0], xl[-1] = 0.0, 1.0

    y_le_avg = 0.5 * (yu[0] + yl[0])
    yu[0] = y_le_avg
    yl[0] = y_le_avg

    def safe_tck(xb, yb):
        m = len(xb)
        k = min(3, m-1)
        if k < 1:
            return None
        return splrep(xb, yb, s=smooth, k=k) if HAS_SCIPY else None

    tck_u = safe_tck(xu, yu)
    tck_l = safe_tck(xl, yl)
    xg = np.linspace(0.0, 1.0, ngrid)

    if tck_u is None:
        yu_g = np.interp(xg, xu, yu)
    else:
        yu_g = splev(xg, tck_u)

    if tck_l is None:
        yl_g = np.interp(xg, xl, yl)
    else:
        yl_g = splev(xg, tck_l)

    yu_g = np.asarray(yu_g, dtype=float)
    yl_g = np.asarray(yl_g, dtype=float)

    y_le_grid = 0.5 * (yu_g[0] + yl_g[0])
    yu_g[0] = y_le_grid
    yl_g[0] = y_le_grid

    swap = yu_g < yl_g
    if np.any(swap):
        tmp = yu_g[swap].copy()
        yu_g[swap] = yl_g[swap]
        yl_g[swap] = tmp

    return xg, yu_g, yl_g

def make_closed_outline_from_envelopes(xg, yu, yl, te_flat=2e-4, scale_chord=1.0, close_le=True, le_radius_pts=40):
    # Handle trailing edge flattening first (in normalized 0-1 space)
    if te_flat > 0:
        x0 = 1.0 - te_flat
        mask = xg <= x0
        x2 = np.r_[xg[mask], x0, 1.0]
        avg_te = 0.5 * (float(yu[-1]) + float(yl[-1]))
        yu2 = np.r_[yu[mask], float(yu[mask][-1]), avg_te]
        yl2 = np.r_[yl[mask], float(yl[mask][-1]), avg_te]
        xu = x2; yu = yu2; yl = yl2
    else:
        xu = xg.copy()
        yu = yu.copy()
        yl = yl.copy()

    # Force LE points to be identical
    if len(xu) > 0:
        y_le_avg = 0.5 * (yu[0] + yl[0])
        yu[0] = y_le_avg
        yl[0] = y_le_avg
        
        # Create smooth LE closure with circular arc
        if close_le and len(xu) > 3:
            # Estimate leading edge radius from the airfoil thickness
            dy_upper = abs(yu[1] - yu[0])  # upper surface deviation
            dy_lower = abs(yl[0] - yl[1])  # lower surface deviation
            dx = xu[1] - xu[0]
            
            thickness_at_1 = abs(yu[1] - yl[1])
            
            if thickness_at_1 > 1e-9 and dx > 1e-9:
                # Estimate circle radius that will smoothly connect both surfaces
                # The radius should be proportional to the local curvature
                r_est = max(thickness_at_1 * 0.5, dx * 0.5)
                
                # Calculate the center of the circle (behind the LE)
                x_center = xu[0] - r_est
                y_center = y_le_avg
                
                # Calculate angles to connect to the first real points on upper and lower surfaces
                angle_to_upper = np.arctan2(yu[1] - y_center, xu[1] - x_center)
                angle_to_lower = np.arctan2(yl[1] - y_center, xu[1] - x_center)
                
                # Create arc that goes from lower surface around to upper surface
                # The arc should span from the lower angle, around the nose, to the upper angle
                n_arc = max(le_radius_pts, 30)
                
                # Ensure we go the correct direction (counterclockwise from lower to upper)
                if angle_to_lower > angle_to_upper:
                    angles = np.linspace(angle_to_lower, angle_to_upper + 2*np.pi, n_arc)
                else:
                    angles = np.linspace(angle_to_lower, angle_to_upper, n_arc)
                
                x_arc = x_center + r_est * np.cos(angles)
                y_arc = y_center + r_est * np.sin(angles)
                
                # Build complete outline: 
                # Start with lower surface (reversed, from TE to point 1)
                # Then arc (from lower point 1 to upper point 1)
                # Then upper surface (from point 1 to TE)
                X_norm = np.concatenate([
                    xu[-1:0:-1],  # lower surface reversed (TE to just after LE)
                    x_arc,         # arc connecting lower to upper at LE
                    xu[1:]         # upper surface (just after LE to TE)
                ])
                Y_norm = np.concatenate([
                    yl[-1:0:-1],  # lower surface reversed
                    y_arc,         # arc
                    yu[1:]         # upper surface
                ])
            else:
                # Fallback: simple closure
                X_norm = np.concatenate([xu, xu[-1:0:-1]])
                Y_norm = np.concatenate([yu, yl[-1:0:-1]])
        else:
            # Simple closure without arc
            X_norm = np.concatenate([xu, xu[-1:0:-1]])
            Y_norm = np.concatenate([yu, yl[-1:0:-1]])
    else:
        X_norm = np.array([])
        Y_norm = np.array([])
    
    # Scale to actual chord
    X = X_norm * scale_chord
    Y = Y_norm * scale_chord
    
    # Remove any duplicate consecutive points
    if len(X) > 1:
        diffs = np.sqrt(np.diff(X)**2 + np.diff(Y)**2)
        keep = np.concatenate([[True], diffs > 1e-12])
        X = X[keep]
        Y = Y[keep]
    
    # Force exact closure
    if len(X) > 0:
        if not np.allclose([X[0], Y[0]], [X[-1], Y[-1]], atol=1e-12):
            X = np.append(X, X[0])
            Y = np.append(Y, Y[0])
    
    return np.column_stack([X, Y])

def sanitize_polygon_xy(xy: np.ndarray, tol_dedup: float = 1e-9, simplify_tol: float = 1e-9, max_attempts: int = 6):
    pts = [(float(x), float(y)) for x, y in xy]

    def dedup(pts_in):
        out = []
        prev = None
        for x, y in pts_in:
            if prev is None or (abs(x - prev[0]) > tol_dedup or abs(y - prev[1]) > tol_dedup):
                out.append((x, y))
                prev = (x, y)
        if out and out[0] != out[-1]:
            out.append(out[0])
        return out

    pts = dedup(pts)
    if len(pts) < 4:
        raise ValueError("Airfoil polygon has too few unique points after dedup.")

    for attempt in range(1, max_attempts+1):
        ring = LinearRing(pts)
        poly = Polygon(ring)
        if poly.is_valid and isinstance(poly, Polygon):
            coords = np.array(poly.exterior.coords)
            if coords.shape[0] >= 4:
                return coords
        try:
            fixed = ops.make_valid(poly)
            if isinstance(fixed, Polygon):
                poly = fixed
            else:
                polys = [g for g in getattr(fixed, "geoms", []) if isinstance(g, Polygon)]
                if polys:
                    poly = max(polys, key=lambda p: p.area)
            if poly.is_valid and isinstance(poly, Polygon):
                coords = np.array(poly.exterior.coords)
                if coords.shape[0] >= 4:
                    return coords
        except Exception:
            pass
        try:
            poly2 = poly.buffer(0)
            if poly2.is_valid and isinstance(poly2, Polygon):
                coords = np.array(poly2.exterior.coords)
                if coords.shape[0] >= 4:
                    return coords
        except Exception:
            pass

    raise ValueError(f"Failed to sanitize airfoil polygon after {max_attempts} attempts")

# ---------------------- GEOMETRY UTILITIES ----------------------

def thickness_at(xval: float, upper: np.ndarray, lower: np.ndarray) -> float:
    y_u = np.interp(xval, upper[:, 0], upper[:, 1])
    y_l = np.interp(xval, lower[:, 0], lower[:, 1])
    return max(0.0, y_u - y_l)

def fos_normal(sig: np.ndarray, Xt: float, Xc: float, max_fos: float = 1e5) -> np.ndarray:
    sig = np.asarray(sig, dtype=float)
    out = np.full(sig.shape, np.nan, dtype=float)
    finite = np.isfinite(sig)
    if not finite.any():
        return out
    pos = (sig >= 0) & finite
    neg = (sig < 0) & finite
    out[pos] = Xt / np.maximum(sig[pos], 1e-12)
    out[neg] = Xc / np.maximum(-sig[neg], 1e-12)
    out = np.clip(out, 0.0, max_fos)
    return out

def shear_partition_GA(V, n_webs, h_web, t_web, G_web, 
                       b_eff, t_core, G_core,
                       webs_active=True):
    """
    Partition transverse shear V between webs and core using stiffness ratios.
    
    Uses parallel spring model: each mechanism carries load proportional to its
    shear stiffness GA = G * Area. This avoids circular logic of extracting forces
    from VQ/It stress fields computed with full V.
    
    Args:
        V: Total transverse shear force (N)
        n_webs: Number of webs
        h_web: Web height (m)
        t_web: Web thickness (m)
        G_web: Web shear modulus (Pa)
        b_eff: Core effective width (m)
        t_core: Core thickness (m)
        G_core: Core shear modulus (Pa)
        webs_active: Whether webs exist/are structural
    
    Returns:
        (V_web_total, V_core): Forces in N
    """
    V_abs = abs(float(V))
    
    if not webs_active or n_webs <= 0:
        return 0.0, V_abs
    
    # Shear stiffness: GA = G * Area
    GA_web = n_webs * G_web * t_web * h_web
    GA_core = G_core * t_core * b_eff
    GA_total = GA_web + GA_core
    
    if GA_total <= 0:
        return 0.0, V_abs
    
    # Partition by stiffness ratio (parallel springs)
    V_web_total = V_abs * (GA_web / GA_total)
    V_core = V_abs - V_web_total
    
    # Safety clamps
    V_web_total = max(0.0, min(V_web_total, V_abs))
    V_core = max(0.0, V_abs - V_web_total)
    
    return V_web_total, V_core

# ---------------------- CLT / COMPOSITE UTILITIES ----------------------

import numpy.linalg as la

def compute_wing_weight(span_m, chord_m, mat, inp, h1_web, h2_web, sb,
                        skin_gsm_per_ply=200.0,
                        web_dry_gsm=200.0,
                        web_aw_is_biax_total=True,
                        rho_f=1780.0, rho_r=1150.0, Vf=0.55):
    """
    Returns a DataFrame 'Weight_Summary' for the FULL wing (both halves), using:
      - Skins by areal mass (skin_gsm_per_ply, plies from mat.face_stack)
      - Webs by FCIM255 dry areal weight (web_dry_gsm) -> cured via Vf model
      - Rods by cross-section and density
      - Core by volume and density
    """
    import pandas as pd
    # --- planform / materials ---
    A_planform = float(span_m) * float(chord_m)          # full wing area [m^2]
    t_core     = float(mat.t_core)
    rho_core   = float(mat.rho_core)

    # --- skins (areal mass per ply, impregnated) ---
    n_skin_plies_per_face = int(sum(int(n) for _, n in mat.face_stack))
    skin_areal_kg_m2_per_face = (skin_gsm_per_ply * n_skin_plies_per_face) / 1e3
    mass_skins_kg = (2.0 * skin_areal_kg_m2_per_face) * A_planform

    # --- core ---
    vol_core_m3  = A_planform * t_core
    mass_core_kg = vol_core_m3 * rho_core

    # cured gsm for ONE biax piece (± pair) at Vf
    resin_over_fiber = ((1 - Vf) * rho_r) / (Vf * rho_f)
    web_dry_each = web_dry_gsm if web_aw_is_biax_total else 2.0 * web_dry_gsm
    web_cured_gsm_per_layer = web_dry_each * (1.0 + resin_over_fiber)
    
    # area for both webs, one side
    A_webs_m2 = float(span_m) * (float(h1_web) + float(h2_web))  # m^2
    
    # layers in your stack per piece (usually 1 unless you stack multiples)
    n_web_layers = max(1, len(mat.web_stack))
    
    # you said TWO separate pieces per web (left/right faces), so factor 2
    n_pieces_per_web = 2
    
    mass_webs_kg = n_pieces_per_web * n_web_layers * (web_cured_gsm_per_layer / 1e3) * A_webs_m2

    # --- rods (spanwise) ---
    n_rods_top   = int(sb.get('n_rods_top', 2))
    n_rods_bot   = int(sb.get('n_rods_bot', 2))
    n_rods_total = n_rods_top + n_rods_bot
    A_rod_cs     = float(inp.rod_width_m) * float(inp.rod_height_m)
    vol_rods_m3  = n_rods_total * float(span_m) * A_rod_cs
    mass_rods_kg = vol_rods_m3 * float(mat.rho_rod)

    # --- totals ---
    mass_total_kg = mass_skins_kg + mass_core_kg + mass_webs_kg + mass_rods_kg
    g0 = 9.80665
    weight_total_N  = mass_total_kg * g0
    weight_total_lb = weight_total_N / 4.4482216153

    rows = [
        ("Planform area (full) [m^2]", A_planform),
        ("Skins mass [kg]",            mass_skins_kg),
        ("Core mass [kg]",             mass_core_kg),
        ("Webs mass [kg]",             mass_webs_kg),
        ("Rods mass [kg]",             mass_rods_kg),
        ("TOTAL mass [kg]",            mass_total_kg),
        ("TOTAL weight [N]",           weight_total_N),
        ("TOTAL weight [lb]",          weight_total_lb),
        ("Half-wing mass [kg]",        0.5 * mass_total_kg),
        ("Skin gsm/ply (impregnated)", skin_gsm_per_ply),
        ("Web dry gsm label",          web_dry_gsm),
        ("AW is total for biax?",      web_aw_is_biax_total),
        ("Vf (web)",                   Vf),
        ("Web cured gsm/layer",        web_cured_gsm_per_layer),
        ("Web layers (from stack)",    n_web_layers),
        ("Rod count (top+bot)",        n_rods_total),
        ("Rod CS [mm x mm]",           f"{inp.rod_width_m*1e3:.1f} x {inp.rod_height_m*1e3:.1f}"),
        ("Core t [mm], ρ [kg/m^3]",    f"{t_core*1e3:.2f}, {rho_core:.1f}"),
    ]
    df = pd.DataFrame(rows, columns=["Parameter", "Value"])
    return df

def ply_stiffness(E1, E2, G12, nu12):
    nu21 = (E2 * nu12) / E1 if E1 != 0 else 0.0
    Q = np.zeros((3,3), dtype=float)
    denom = 1.0 - nu12 * nu21 if abs(1.0 - nu12 * nu21) > 1e-12 else 1e-12
    Q[0,0] = E1 / denom
    Q[1,1] = E2 / denom
    Q[0,1] = (nu12 * E2) / denom
    Q[1,0] = Q[0,1]
    Q[2,2] = G12
    return Q

def transform_Q(Q, theta_deg):
    th = math.radians(theta_deg)
    c = math.cos(th); s = math.sin(th)
    c2 = c*c; s2 = s*s; cs = c*s
    
    Q11 = float(Q[0,0]); Q12 = float(Q[0,1])
    Q22 = float(Q[1,1]); Q66 = float(Q[2,2])
    
    Qbar = np.zeros((3,3), dtype=float)
    Qbar[0,0] = Q11*c2*c2 + 2.0*(Q12 + 2.0*Q66)*c2*s2 + Q22*s2*s2
    Qbar[0,1] = (Q11 + Q22 - 4.0*Q66)*c2*s2 + Q12*(c2*c2 + s2*s2)
    Qbar[1,1] = Q11*s2*s2 + 2.0*(Q12 + 2.0*Q66)*c2*s2 + Q22*c2*c2
    Qbar[0,2] = (Q11 - Q12 - 2.0*Q66)*c*c2*s + (Q12 - Q22 + 2.0*Q66)*s*s2*c
    Qbar[1,2] = (Q11 - Q12 - 2.0*Q66)*c*s2*s + (Q12 - Q22 + 2.0*Q66)*c*c2*s
    # CORRECTED formula:
    Qbar[2,2] = (Q11 + Q22 - 2.0*Q12 - 2.0*Q66)*c2*s2 + Q66*(c2 - s2)**2
    
    # Symmetry
    Qbar[1,0] = Qbar[0,1]
    Qbar[2,0] = Qbar[0,2]
    Qbar[2,1] = Qbar[1,2]
    
    return Qbar

def laminate_ABD_from_stack(stack: List[Tuple[float,int]], ply_props: Dict[str,float], t_ply):
    # stack: list of (theta_deg, nplies) for a single face (order will be used top->bottom)
    ply_list = []
    for theta, n in stack:
        for _ in range(int(n)):
            ply_list.append(float(theta))
    N = len(ply_list)
    total_t = N * t_ply
    if N == 0:
        A = np.zeros((3,3)); B = np.zeros((3,3)); D = np.zeros((3,3)); return A,B,D,total_t
    # z coordinates top-to-bottom (midplane at 0)
    z = np.linspace(total_t/2.0, -total_t/2.0, N+1)
    A = np.zeros((3,3), dtype=float)
    B = np.zeros((3,3), dtype=float)
    D = np.zeros((3,3), dtype=float)
    E1 = ply_props.get('E1', ply_props.get('E1_fiber', 1e9))
    E2 = ply_props.get('E2', ply_props.get('E2_fiber', 1e9))
    G12 = ply_props.get('G12', ply_props.get('G12_fiber', 1e9))
    nu12 = ply_props.get('nu12', ply_props.get('nu12_fiber', 0.3))
    Qply = ply_stiffness(E1, E2, G12, nu12)
    for k in range(N):
        theta = ply_list[k]
        Qbar = transform_Q(Qply, theta)
        z_k = z[k]; z_k1 = z[k+1]
        A += Qbar * (z_k - z_k1)
        B += 0.5 * Qbar * (z_k**2 - z_k1**2)
        D += (1.0/3.0) * Qbar * (z_k**3 - z_k1**3)
    return A, B, D, total_t

def effective_bending_EI_from_D(D_mat):
    EI_per_width = float(D_mat[0,0])
    return EI_per_width

def laminate_membrane_G(A_mat, total_t):
    """Extract membrane shear modulus from ABD matrix"""
    t = max(total_t, 1e-9)
    A66 = float(A_mat[2,2])
    if A66 <= 0:
        print(f"WARNING: Non-positive A66={A66}, using fallback G=1e6")
        return 1e6
    G_membrane = A66 / t
    return max(G_membrane, 1e6)

def laminate_inplane_axial_modulus(A_mat, total_t):
    """
    Extract in-plane axial modulus from ABD matrix A.
    
    E_x = A11 / t  (membrane stiffness per unit width / thickness)
    
    Args:
        A_mat: 3x3 membrane stiffness matrix from CLT
        total_t: Total laminate thickness (m)
    
    Returns:
        E_axial: In-plane axial modulus (Pa)
    """
    t = max(total_t, 1e-9)
    A11 = float(A_mat[0,0])
    E_axial = A11 / t
    return max(E_axial, 1e6)  # Safety floor

def tsai_wu_index(sx, sy, txy, Xt, Xc, Yt=None, Yc=None, S=None):
    # 2D Tsai-Wu conservative wrapper. If Yt/Yc or S missing, use proxies.
    if Yt is None: Yt = Xt * 0.6
    if Yc is None: Yc = Xc * 0.6
    if S is None: S = min(abs(Xt), abs(Xc)) * 0.1
    F1 = 1.0/Xt - 1.0/Xc
    F2 = 1.0/Yt - 1.0/Yc
    F11 = 1.0/(Xt*Xc)
    F22 = 1.0/(Yt*Yc)
    F66 = 1.0/(S*S)
    F12 = -0.5 * math.sqrt(abs(F11 * F22)) if F11*F22 >= 0 else 0.0
    val = F1 * sx + F2 * sy + F11 * sx * sx + F22 * sy * sy + 2.0 * F12 * sx * sy + F66 * txy * txy
    return val

def sandwich_wrinkling_stress(E_face_axial, E_core, G_core, nu_face=0.3, c_coeff=0.5):
    """
    Face-sheet wrinkling stress for sandwich panel (continuous core support).
    
    Simplified cube-root formula (first-order screen):
    σ_wr = c_coeff * (E_face * E_core * G_core / (1 - nu_face²))^(1/3)
    
    Reference: Allen (1969), Zenkert (1995)
    Typical c_coeff range: 0.4-0.6 depending on assumptions
    (symmetric vs antisymmetric wrinkling, isotropic vs orthotropic face)
    
    Args:
        E_face_axial: Face sheet axial modulus in compression direction (Pa)
        E_core: Core compression modulus (Pa)
        G_core: Core shear modulus (Pa)
        nu_face: Face sheet Poisson's ratio
        c_coeff: Wrinkling coefficient (0.5 default, range 0.4-0.6)
    
    Returns:
        sigma_wr: Wrinkling stress (Pa), first-order estimate
    """
    E_eff_face = E_face_axial / (1.0 - nu_face**2)
    sigma_wr = c_coeff * (E_eff_face * E_core * G_core) ** (1.0/3.0)
    return sigma_wr

def web_shear_buckling_stress_infinite_strip(E_web, h_web, t_web, nu_web=0.3):
    """
    Web panel shear buckling stress (infinite-strip model, no spanwise restraint).
    
    For continuous web without discrete ribs/bulkheads, use infinite-strip coefficient.
    Formula: τ_cr = k_s * π² * E / (12(1-ν²)) * (t/h)²
    
    Infinite strip (a/b → ∞): k_s = 5.35 + 4*(b/a)² → 5.35 as a → ∞
    
    Args:
        E_web: Web modulus (Pa)
        h_web: Web clear height between supports (m)
        t_web: Web thickness (m)
        nu_web: Poisson's ratio
    
    Returns:
        tau_cr: Critical shear buckling stress (Pa)
    
    Note: Conservative for continuous sandwich (no discrete panels).
    If ribs/bulkheads added, use finite panel k_s and aspect ratio.
    """
    if h_web <= 0 or t_web <= 0:
        return np.inf
    
    k_s = 5.35  # Infinite-strip coefficient (conservative lower bound)
    E_eff = E_web / (12.0 * (1.0 - nu_web**2))
    tau_cr = k_s * np.pi**2 * E_eff * (t_web / h_web)**2
    return tau_cr

# ---------------------- STRUCTURAL HELPERS (CLT-enabled) ----------------------

def sandwich_bending_properties(chord: float,
                                mat: Materials,
                                rod_w: float, rod_h: float, gap: float,
                                n_rods_top: int = 2, n_rods_bot: int = 2,
                                ply_props: Optional[Dict]=None) -> Dict:
    if ply_props is None:
        ply_props = {'E1': mat.E1_fiber, 'E2': mat.E2_fiber, 'G12': mat.G12_fiber, 'nu12': mat.nu12_fiber}

    # Build CLT for a face (face_stack)
    ply_props_skin = {'E1': mat.E1_skin, 'E2': mat.E2_skin, 'G12': mat.G12_skin, 'nu12': mat.nu12_skin}
    A_face, B_face, D_face, t_face_total = laminate_ABD_from_stack(mat.face_stack, ply_props_skin, mat.t_ply_skin)
    EI_face = effective_bending_EI_from_D(D_face)
    G_skin_eff = laminate_membrane_G(A_face, t_face_total)

    # core self inertia
    I_core_self = (chord * mat.t_core ** 3) / 12.0
    EI_core = mat.E_core * I_core_self

    # rods
    A_rod = rod_w * rod_h
    I_rod_self = (rod_w * rod_h ** 3) / 12.0
    # skin mid-surface location measured from midplane (positive above midplane)
    # skins sit outside the foam: skin centroid is core half-thickness + half skin thickness
    z_face = (mat.t_core / 2.0) + (t_face_total / 2.0)

    # Proper parallel-axis for skins
    EI_skins = 2.0 * (EI_face + chord * t_face_total * z_face ** 2)

    # Place rods inside the foam core, not outside the skins.
    # rod centers should be located inside the foam at a distance from midplane:
    #   z_rod = (t_core/2) - (rod_h/2) - gap
    # where 'gap' is the clearance from the inner skin surface into the foam (m).
    # Clamp to a small value if the geometry would push the rod outside the foam.
    # Place rods inside the foam core: rod centroid measured from midplane
    # desired position: just inside the inner skin surface by 'gap' (gap is clearance from inner skin into foam)
    # --- rod geometry and placement (force rods to be inside the foam, flush with inner skin) ---
    core_half = mat.t_core / 2.0
    rod_half = float(rod_h) / 2.0
    gap_val = max(0.0, float(gap))
    # tiny margin so rods do not touch the inner skin
    margin = min(core_half * 0.01, 1e-6)

    # If user-specified rod is taller than the foam, clamp the effective rod_half so the rod fits inside.
    if rod_half >= core_half - 1e-12:
        # warn and clamp rod height to sit inside foam (preserve area ratio by scaling geometric dimensions if desired)
        rod_half_eff = max(1e-6, core_half - margin)
        rod_h_eff = 2.0 * rod_half_eff
        # keep original rod_w but reduce rod_h behaviorally for centroid/I calculations
        A_rod = float(rod_w) * rod_h_eff
        I_rod_self = (float(rod_w) * rod_h_eff ** 3) / 12.0
        print(f"DEBUG: requested rod_h ({rod_h:.6g} m) >= core thickness; using rod_h_eff={rod_h_eff:.6g} m for geometry")
    else:
        rod_half_eff = rod_half
        rod_h_eff = 2.0 * rod_half_eff
        A_rod = float(rod_w) * rod_h_eff
        I_rod_self = (float(rod_w) * rod_h_eff ** 3) / 12.0

    # Desired flush centroid: just inside inner skin surface (positive above midplane)
    z_rod_flush = core_half - rod_half_eff - margin
    if z_rod_flush <= 0.0:
        # fallback: place rod centroid half-way into the core (safe internal position)
        z_rod = max(1e-6, core_half * 0.5)
    else:
        z_rod = z_rod_flush

    # symmetric top/bottom
    z_rod_top = +z_rod
    z_rod_bot = -z_rod

    # rod contributions (parallel-axis) and include into EI
    EI_rods_top = float(n_rods_top) * float(mat.E_rod) * (I_rod_self + A_rod * z_rod_top**2)
    EI_rods_bot = float(n_rods_bot) * float(mat.E_rod) * (I_rod_self + A_rod * z_rod_bot**2)
    EI_rods = EI_rods_top + EI_rods_bot

    # update EI_total using rods
    EI_total = max(1e-12, EI_skins + EI_core + EI_rods)
    
    # Debug output to show stiffness breakdown
    print(f"\n=== BENDING STIFFNESS BREAKDOWN ===")
    print(f"EI_skins = {EI_skins:.6e} N·m² ({100*EI_skins/EI_total:.1f}%)")
    print(f"EI_core  = {EI_core:.6e} N·m² ({100*EI_core/EI_total:.1f}%)")
    print(f"EI_rods  = {EI_rods:.6e} N·m² ({100*EI_rods/EI_total:.1f}%)")
    print(f"EI_total = {EI_total:.6e} N·m²")
    print(f"====================================\n")

    return {
        'EI_total': EI_total,
        'EI_skins': EI_skins,
        'EI_core': EI_core,
        'EI_rods': EI_rods,
        'z_skin_top': z_face,
        'z_skin_bot': -z_face,
        'z_rod_top': z_rod_top,
        'z_rod_bot': z_rod_bot,
        'A_skin': chord * t_face_total,
        'A_rod': A_rod,
        'd_rod': 2 * z_rod,
        'n_rods_top': n_rods_top,
        'n_rods_bot': n_rods_bot,
        't_skin_total': t_face_total,
        'G_skin_eff': G_skin_eff,
        'D_face': D_face,
        'A_face': A_face
    }

def assemble_bredt_matrix(A_cells: List[float], walls: List[Dict], G_vals: List[float]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Assemble the Bredt system for multi-cell torsion.
    Unknowns: q_j (shear flow on each wall)
    Equations: (nc-1) compatibility rows + 1 torque equilibrium row
    Returns matrix M and rhs vector b where M q = b.
    """
    nw = len(walls)
    nc = len(A_cells)
    if nc == 0 or nw == 0:
        return np.zeros((0, nw)), np.zeros((0,))

    # s_j = L_j / (G_j * t_j)
    s = np.zeros(nw)
    for j, w in enumerate(walls):
        G = G_vals[j] if j < len(G_vals) else G_vals[-1]
        s[j] = float(w['length']) / (max(float(G), 1e-12) * max(float(w['thickness']), 1e-12))

    # loop matrix: loop[c, j] = 1 if wall j bounds cell c, else 0
    loop = np.zeros((nc, nw))
    for j, w in enumerate(walls):
        for c in w['cell_ids']:
            if 0 <= c < nc:
                loop[c, j] = 1.0

    if nc == 1:
        M = np.zeros((1, nw))
        torque_coeffs = np.zeros((nw,))
        for j, w in enumerate(walls):
            torque_coeffs[j] = 2.0 * A_cells[0] * (1.0 if loop[0, j] else 0.0)
        M[0, :] = torque_coeffs
        b = np.array([0.0], dtype=float)
        return M, b

    Mcomp = np.zeros((nc - 1, nw))
    for i in range(1, nc):
        Mcomp[i - 1, :] = loop[i, :] * s - loop[0, :] * s

    # torque row constructed outside
    return Mcomp, np.zeros((Mcomp.shape[0],))

def multicell_bredt_solver(A_cells: List[float], walls: List[Dict], G_vals: List[float], T: float) -> Dict:
    """
    Solve for wall shear flows q_j given torque T using Bredt multi-cell approach.
    Returns q (array length nw), tau_wall (abs shear stress), and status info.
    """
    nw = len(walls)
    nc = len(A_cells)
    if abs(T) < 1e-12 or nc == 0 or nw == 0:
        return {'q': np.zeros(nw), 'tau_wall': np.zeros(nw), 'status': 'zero_T_or_no_cells'}

    # s_j values:
    s = np.zeros(nw)
    for j, w in enumerate(walls):
        G = G_vals[j] if j < len(G_vals) else G_vals[-1]
        s[j] = float(w['length']) / (max(float(G), 1e-12) * max(float(w['thickness']), 1e-12))

    # loop matrix
    loop = np.zeros((nc, nw))
    for j, w in enumerate(walls):
        for c in w['cell_ids']:
            if 0 <= c < nc:
                loop[c, j] = 1.0

    # Build compatibility rows
    if nc == 1:
        # Single closed cell: constant shear flow q0 = T/(2*A_cell)
        q0 = float(T) / (2.0 * float(A_cells[0]))
        q = np.zeros(nw, dtype=float)
        
        # Assign q0 to all walls on the loop (loop membership is 0/1, not signed)
        for j in range(nw):
            q[j] = q0 if loop[0, j] else 0.0  # Defensive: only if wall bounds cell 0
        
        # Torsion shear stress per wall
        tau_wall = np.zeros(nw, dtype=float)
        for j, w in enumerate(walls):
            t = max(float(w['thickness']), 1e-6)
            tau_wall[j] = abs(q[j]) / t
        
        # Return cond=1.0 as dummy (direct calculation, not matrix solve)
        return {'q': q, 'tau_wall': tau_wall, 'status': 'bredt_solved', 'cond': 1.0}
    else:
        Mcomp = np.zeros((nc - 1, nw))
        for i in range(1, nc):
            Mcomp[i - 1, :] = loop[i, :] * s - loop[0, :] * s

        torque_coeffs = np.zeros((nw,))
        for j in range(nw):
            acc = 0.0
            for c in walls[j]['cell_ids']:
                if 0 <= c < nc:
                    acc += A_cells[c]
            torque_coeffs[j] = 2.0 * acc
        M = np.vstack([Mcomp, torque_coeffs.reshape(1, -1)])
        b = np.zeros((M.shape[0],))
        b[-1] = T

    try:
        q, *_ = np.linalg.lstsq(M, b, rcond=None)
    except Exception:
        reg = 1e-12
        MtM = M.T.dot(M) + reg * np.eye(nw)
        q = np.linalg.solve(MtM, M.T.dot(b))

    tau_wall = np.zeros(nw)
    for j, w in enumerate(walls):
        t = max(float(w['thickness']), 1e-6)
        tau_wall[j] = abs(q[j]) / t

    cond = np.nan
    try:
        cond = np.linalg.cond(M)
    except Exception:
        cond = np.nan

    status = 'bredt_solved'
    if cond is not None and not np.isnan(cond) and cond > 1e12:
        status = 'ill_conditioned'

    return {'q': q, 'tau_wall': tau_wall, 'status': status, 'cond': cond}

def convert_to_single_cell_topology(A_cells: List[float], walls: List[Dict], G_vals: List[float]) -> Tuple[List[float], List[Dict], List[float]]:
    """
    Convert multi-cell topology to single-cell by merging all cells
    and removing internal walls (keep only outer perimeter walls).
    
    Internal webs have cell_ids like [1] or [2] (single cell only).
    Outer walls have cell_ids like [0,1] or [1,2] (shared between cells).
    
    Returns:
        A_single: List with one merged cell area
        walls_single: List of walls forming outer perimeter only
        G_single: List of G values corresponding to walls_single
    """
    # Merge all cells into one
    A_single = [sum(A_cells)]
    
    # Filter walls: keep only those that bound multiple cells (outer perimeter)
    walls_single = []
    G_single = []
    for j, w in enumerate(walls):
        cell_ids = w.get('cell_ids', [])
        # Keep if wall bounds 2+ cells OR is the only cell
        if len(cell_ids) >= 2 or (len(cell_ids) == 1 and len(A_cells) == 1):
            walls_single.append({
                'length': w['length'],
                'thickness': w['thickness'],
                'cell_ids': [0]  # All reference the single merged cell
            })
            if j < len(G_vals):
                G_single.append(G_vals[j])
    
    return A_single, walls_single, G_single

def split_polygon_by_xcuts(xy: np.ndarray, xcuts: List[float]) -> List[np.ndarray]:
    if xy is None or len(xy) < 4:
        return []
    poly = Polygon(xy)
    if not poly.is_valid:
        try:
            poly = ops.make_valid(poly)
        except Exception:
            try:
                poly = poly.buffer(0)
            except Exception:
                pass

    cutters = [LineString([(float(xc), -1e6), (float(xc), 1e6)]) for xc in sorted(xcuts)]
    multi = MultiLineString(cutters)
    try:
        pieces = ops.split(poly, multi)
        geoms = list(getattr(pieces, "geoms", [pieces]))
    except Exception:
        geoms = [poly]

    out = []
    for g in geoms:
        if isinstance(g, Polygon) and g.area > 1e-15:
            coords = np.array(g.exterior.coords)
            if coords.shape[0] >= 4:
                if not np.allclose(coords[0], coords[-1], atol=1e-12):
                    coords = np.vstack([coords, coords[0]])
                out.append(coords.astype(float))

    if not out:
        try:
            if isinstance(poly, Polygon) and poly.is_valid and poly.area > 1e-15:
                coords = np.array(poly.exterior.coords)
                if coords.shape[0] >= 4:
                    out = [coords.astype(float)]
        except Exception:
            out = []

    out.sort(key=lambda arr: float(np.min(arr[:, 0])) if arr.size else 0.0)
    return out

# ---------------------- DIAGNOSTICS ----------------------

def report_low_fos(name: str, arr: np.ndarray, span_arr: np.ndarray, threshold: float = 3.0, n_show: int = 8):
    arr = np.asarray(arr, dtype=float)
    finite = np.isfinite(arr)
    if not finite.any():
        print(f"{name}: no finite values")
        return
    mn = float(np.nanmin(arr))
    mean = float(np.nanmean(arr))
    mx = float(np.nanmax(arr))
    print(f"{name}: min {mn:.3g}, mean {mean:.3g}, max {mx:.3g}")
    low_idx = np.where((finite) & (arr <= threshold))[0]
    if low_idx.size == 0:
        print(f"  no values <= {threshold}")
        return
    worst = low_idx[np.argsort(arr[low_idx])][:n_show]
    print(f"  worst indices (<= {threshold}):")
    for i in worst:
        print(f"    i={int(i)}, x={span_arr[int(i)]:.4f} m, {name}={arr[int(i)]:.3g}")

# ---------------------- MAIN ANALYSIS ----------------------

def run_analysis_professional(inp: Inputs, mat: Materials) -> Dict[str, pd.DataFrame]:
    # Validate critical inputs
    if inp.chord_m is None or inp.chord_m <= 0:
        raise ValueError("Chord (inp.chord_m) must be provided and > 0. Do not rely on AUW-based chord derivation.")
    if not (0.0 <= inp.web1_p <= 1.0 and 0.0 <= inp.web2_p <= 1.0):
        raise ValueError("web1_p and web2_p must be between 0 and 1.")
    if inp.web2_p <= inp.web1_p:
        raise ValueError("web2_p must be greater than web1_p (webs must be ordered root->tip in x/c).")

    # Load and process airfoil
    pts = load_airfoil_points(inp.airfoil_file)
    loop = normalize_and_order_loop(pts)
    xg, yu_g, yl_g = build_smooth_envelopes(loop, ngrid=601, smooth=0.0)

    span_m = inp.span_ft * 0.3048
    chord_m = inp.chord_m

    raw_outline = make_closed_outline_from_envelopes(
        xg, yu_g, yl_g, 
        te_flat=2e-4, 
        scale_chord=chord_m, 
        close_le=True, 
        le_radius_pts=40
    )

    # define local thickness and web cut positions
    t_face = mat.t_skin_total
    t_core = mat.t_core
    t_web = mat.t_web_total
    p1 = float(np.clip(inp.web1_p, 0.0, 1.0))
    p2 = float(np.clip(inp.web2_p, 0.0, 1.0))
    x_cut1 = p1 * chord_m
    x_cut2 = p2 * chord_m

    # print diagnostics
    print(f"DEBUG: chord={chord_m:.6f} m, span={span_m:.4f} m, sweep={inp.sweep_deg:.3f} deg")
    print(f"DEBUG: web positions x1={x_cut1:.6f} m, x2={x_cut2:.6f} m")
    print(f"DEBUG: t_face={t_face:.6g} m, t_core={t_core:.6g} m, t_web={t_web:.6g} m")
    # compute skin CLT G estimate placeholder (we'll compute properly after sb)
    try:
        xy_poly = sanitize_polygon_xy(raw_outline)
    except Exception as e:
        raise RuntimeError(f"sanitize_polygon_xy failed: {e}")

    upper = np.column_stack([xg * chord_m, yu_g * chord_m])
    lower = np.column_stack([xg * chord_m, yl_g * chord_m])
    
    # --- Geometry-derived scalars and CLT sandwich properties (must come before any use of sb/EI_total) ---
    # web cut positions (already computed earlier as x_cut1/x_cut2, but recompute to be safe)
    p1 = float(np.clip(inp.web1_p, 0.0, 1.0))
    p2 = float(np.clip(inp.web2_p, 0.0, 1.0))
    x_cut1 = p1 * chord_m
    x_cut2 = p2 * chord_m
    
    # local thicknesses at web cuts (airfoil geometry)
    h1 = thickness_at(x_cut1, upper, lower)
    h2 = thickness_at(x_cut2, upper, lower)
    h1_web = max(0.0, h1 - 2.0 * t_face)
    h2_web = max(0.0, h2 - 2.0 * t_face)
    
    # Validate web heights
    if h1_web <= 0 or h2_web <= 0:
        raise ValueError(
            f"Web height invalid: h1_web={h1_web*1000:.3f}mm, h2_web={h2_web*1000:.3f}mm. "
            f"Skins (t_face={t_face*1000:.3f}mm) may be too thick for airfoil."
        )
    
    # neutral axis / aerodynamic center fallback for moment arm calculation
    if inp.x_sc_c is not None:
        x_sc_c = float(inp.x_sc_c)
        print(f"DEBUG: using user x_sc_c = {x_sc_c:.4f} c")
    else:
        # Default to quarter-chord aerodynamic center (typical for subsonic airfoils)
        x_sc_c = 0.25
        print(f"DEBUG: x_sc_c not provided, defaulting to quarter-chord x_sc_c = {x_sc_c:.4f} c")
    
    # Compute sandwich/CLT section properties once and for all (used by bending and transforms)
    sb = sandwich_bending_properties(
        chord_m, mat,
        inp.rod_width_m, inp.rod_height_m, inp.rod_gap_mm * 1e-3,
        n_rods_top=2, n_rods_bot=2
    )
    EI_total = sb.get('EI_total', 1e-12)
    
    # Define consistent reference modulus for ALL transformed-section calculations
    # Use CLT-derived effective bending modulus from D matrix
    D_face = sb.get('D_face')
    t_skin_total = sb.get('t_skin_total', mat.t_skin_total)
    if D_face is not None and t_skin_total > 0:
        E_ref = 12.0 * float(D_face[0,0]) / max(t_skin_total**3, 1e-18)
    else:
        E_ref = float(mat.E1_fiber)  # fallback
    print(f"DEBUG: Using E_ref = {E_ref/1e9:.2f} GPa for transformed-section calculations")
    
    # Define explicit rod chordwise positions (fraction of chord)
    # 2 rods forward (near web1), 2 rods aft (near web2)
    rod_positions_xc = [0.2, 0.2, 0.65, 0.65]  # 2 @ 20%, 2 @ 65% chord
    rod_x_positions = [xc * chord_m for xc in rod_positions_xc]
    
    # Count rods forward of each web (for Q calculations)
    n_rods_fwd_web1 = sum(1 for x in rod_x_positions if x <= x_cut1)
    n_rods_fwd_web2 = sum(1 for x in rod_x_positions if x <= x_cut2)
    
    # Split into top/bottom (assume equal distribution)
    n_top_fwd_1 = n_rods_fwd_web1 // 2
    n_bot_fwd_1 = n_rods_fwd_web1 - n_top_fwd_1
    n_top_fwd_2 = n_rods_fwd_web2 // 2
    n_bot_fwd_2 = n_rods_fwd_web2 - n_top_fwd_2
    
    print(f"DEBUG: Rod positions: {rod_positions_xc} (x/c)")
    print(f"DEBUG: Rods forward of web1 ({inp.web1_p:.2f}c): {n_rods_fwd_web1} (top:{n_top_fwd_1}, bot:{n_bot_fwd_1})")
    print(f"DEBUG: Rods forward of web2 ({inp.web2_p:.2f}c): {n_rods_fwd_web2} (top:{n_top_fwd_2}, bot:{n_bot_fwd_2})")

    # === Weight summary (full wing) ===
    df_weight = compute_wing_weight(
        span_m=span_m,
        chord_m=chord_m,
        mat=mat,
        inp=inp,
        h1_web=h1_web,
        h2_web=h2_web,
        sb=sb,
        skin_gsm_per_ply=200.0,     # FCIM255 total fabric weight (±45° biax)
        web_dry_gsm=200.0,          # FCIM255 fabric weight (±45° biax)
        web_aw_is_biax_total=True,  # True: 200 gsm is total for ±45° pair
        Vf=0.55                     # fiber volume fraction
    )
    # make it available to the Excel writer scope
    inp._df_weight = df_weight

    # Keep physical span fixed regardless of sweep
    rad = math.radians(inp.sweep_deg)
    cos_sweep = max(math.cos(rad), 1e-6)

    n = max(2, int(inp.n_stations))
    structural_half = span_m / 2.0
    xspan = np.linspace(0.0, structural_half, n)
    dx = np.zeros_like(xspan)
    dx[1:] = np.diff(xspan)

    # --- aerodynamic input handling (V_ms, rho -> q_pa) and lift selection ---
    
    # 1) Derive q if needed
    q_from_vrho = (
        inp.V_ms is not None and np.isfinite(inp.V_ms) and
        inp.rho  is not None and np.isfinite(inp.rho)
    )
    if (inp.q_pa is None or not np.isfinite(inp.q_pa)) and q_from_vrho:
        inp.q_pa = 0.5 * float(inp.rho) * float(inp.V_ms)**2
        print(f"DEBUG: computed q_pa from V_ms/rho = {inp.q_pa:.6g} Pa (V_ms={inp.V_ms}, rho={inp.rho})")
    
    # Treat non-positive q as “not provided”
    q_pa_local = float(inp.q_pa) if (inp.q_pa is not None and np.isfinite(inp.q_pa) and inp.q_pa > 0) else None
    
    # 2) Spanwise lift SHAPE (dimensionless, half-wing 0 → b/2)
    b_half = span_m / 2.0
    eta = xspan / max(b_half, 1e-12)  # 0 → 1
    
    # Tuned analytic shape to give:
    # - root ~0.9 of peak
    # - broad plateau inboard
    # - decay to zero at geometric tip
    n_shape = 0.20   # very gentle elliptic falloff
    a_shape = 0.80   # mid-span bump → "shoulders"
    
    base = 1.0 - eta**2
    base = np.clip(base, 0.0, None)
    
    shape_raw = base**n_shape * (1.0 + a_shape * eta**2 * base)
    
    # Normalize so max(shape) = 1
    shape_max = np.max(shape_raw)
    shape = shape_raw / shape_max if shape_max > 0.0 else shape_raw
    
    # 3) Build initial lift_per_unit_span; area is normalized later to total_lift_N
    if (q_pa_local is not None) and (inp.cl is not None) and np.isfinite(inp.cl):
        # Aero inputs available → convenient starting magnitude q*c*CL
        lift_per_unit_span = q_pa_local * chord_m * float(inp.cl) * shape   # N/m
        source = (
            "analytic spanwise shape: w(x) ~ q*c*CL*shape(eta); "
            "magnitude normalized to match half-wing lift at g_limit"
        )
    else:
        # No usable q,CL → just use shape, scaled later from AUW/g-limit
        lift_per_unit_span = shape
        source = (
            "analytic spanwise shape (no aero inputs); "
            "magnitude normalized to match half-wing lift at g_limit"
        )
# ====================== END: Half-wing load normalization =====================

    # 4) IMPORTANT: Do NOT multiply lift magnitude by 1/cos(sweep) here.
    # If you want sweep effects, project later in the internal-loads step (e.g., use cos(Λ) on bending component),
    # but keep the TOTAL lift equal to AUW * g_limit.


    # diagnostic print showing what source is used and key values
    print(f"DEBUG lift source: {source}")
    try:
        print(f"DEBUG lift_per_unit_span[0]={float(lift_per_unit_span[0]):.6g} N/m, chord_m={chord_m:.6g} m")
    except Exception:
        pass

    # Use lift offset where available; if it cancels x_sc_c, fall back to using cm_ac to create a pitching moment per unit span
    moment_arm = (inp.lift_offset_c - x_sc_c) * chord_m
    if abs(moment_arm) < 1e-12:
        # produce a moment per unit span from an aerodynamic moment coefficient about aerodynamic center
        # cm_ac is nondimensional moment coefficient per unit span referenced to chord (user-provided)
        moment_arm = (inp.cm_ac) * chord_m  # units: m (this produces a small nonzero lever proportional to chord)
    m_t_lift = lift_per_unit_span * moment_arm
    m_t_cm = np.zeros_like(xspan)
    m_t = m_t_lift + m_t_cm

    # ===================== BEGIN: Half-wing load normalization =====================
    # Goal: ensure ∫_0^{span_m} lift_per_unit_span dx == Target half-wing lift at g_limit
    # This fixes the "half of half" bug and makes root V/M correct.
    
    import numpy as _np
    
    # 1) Figure out target half-wing lift @ g_limit
    # Prefer aerodynamic inputs if available (q, cl); otherwise use wing loading.
    g_limit = getattr(inp, "g_limit", 5.0)  # fallback if not present
    
    # Dynamic pressure q (Pa): use existing inp.q_pa if already computed upstream
    q_pa = None
    try:
        q_pa = float(inp.q_pa) if _np.isfinite(inp.q_pa) else None
    except Exception:
        q_pa = None
    
    have_aero = (q_pa is not None and _np.isfinite(q_pa)
                 and getattr(inp, "cl", None) is not None
                 and _np.isfinite(inp.cl))
    
    if have_aero:
        # Aerodynamic path: full-wing lift at 1 g
        L_full_1g = q_pa * float(inp.cl) * float(chord_m) * float(span_m)
        # Target half-wing lift at g_limit
        L_half_target = 0.5 * L_full_1g * float(g_limit)
    else:
        # Wing-loading path (Inputs has wing_loading_lbft2)
        WL_lb_ft2 = float(getattr(inp, "wing_loading_lbft2", 0.0))
        # 1 lb/ft^2 = 47.88025898 N/m^2
        WL_N_m2   = WL_lb_ft2 * 47.88025898
        A_full    = float(span_m) * float(chord_m)     # full-wing area
        L_full_1g = WL_N_m2 * A_full                           # N @ 1g for full wing
        L_half_target = 0.5 * L_full_1g * float(g_limit)       # target N for one semi-span
    
    # 2) Integrate your presently-built distribution to see what it totals
    # Use trapz for numpy < 1.21.0 compatibility
    try:
        L_half_integrated = _np.trapezoid(lift_per_unit_span, xspan)
    except AttributeError:
        L_half_integrated = _np.trapz(lift_per_unit_span, xspan)
    
    # 3) Scale the distribution so it matches the target exactly
    _eps = 1e-12
    if not _np.isfinite(L_half_integrated) or abs(L_half_integrated) < _eps:
        raise ValueError("lift_per_unit_span integration is zero/invalid; cannot normalize loads.")
    
    scale = L_half_target / L_half_integrated
    lift_per_unit_span = lift_per_unit_span * scale
    
    # 4) Optional sanity prints
    print(f"[LOAD-NORM] Target half-wing lift @ g_limit: {L_half_target:.3f} N")
    try:
        integrated_check = _np.trapezoid(lift_per_unit_span, xspan)
    except AttributeError:
        integrated_check = _np.trapz(lift_per_unit_span, xspan)
    print(f"[LOAD-NORM] Integrated half-wing lift (after scale): {integrated_check:.3f} N")
    print(f"[LOAD-NORM] Applied scale factor: {scale:.6f}")
    # ====================== END: Half-wing load normalization ======================
    # -------- Diagnostic plot: spanwise lift distribution --------
    # -------- Diagnostic plot: NON-normalized spanwise lift distribution --------
    # -------- Plot lift distribution over the full wingspan --------
    if inp.make_plots:
        import matplotlib.pyplot as plt

        # Half-span data (already in meters, 0 → b/2)
        y_half = xspan                  # [m]
        w_half = lift_per_unit_span     # [N/m]

        # Build symmetric full-span arrays: -b/2 → +b/2
        # (skip the root point once to avoid duplicating y=0)
        y_full = np.concatenate((-y_half[::-1], y_half[1:]))
        w_full = np.concatenate(( w_half[::-1], w_half[1:]))

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(y_full, w_full, "-o", linewidth=2)
        ax.set_xlabel("y [m] Along Span (−b/2 to +b/2)")
        ax.set_ylabel("Lift Per Unit Span [N/m]")
        ax.set_title("Spanwise Lift Distribution (Full Wingspan)")
        ax.grid(True)
        fig.tight_layout()
        fig.savefig("lift_distribution_full_span.png", dpi=300)
        #plt.close(fig)
        # for interactive use instead of saving:
        plt.show()

    # --- Sweep-resolved internal loads (geometry unchanged) ---
    # Geometry/ABD/section properties stay fixed. Only the load vector is resolved through sweep.
    # Flapwise bending and shear come from the component L*cos(Λ).
    # Torsion from the chordwise lever of lift also scales with cos(Λ).
    
    rad = math.radians(inp.sweep_deg)
    cos_sweep = max(math.cos(rad), 1e-6)
    sin_sweep = math.sin(rad)  # kept for future in-plane effects if you add them
    
    # Decompose lift per unit span
    w_flap = lift_per_unit_span * cos_sweep    # N/m, drives flapwise V and M
    # Optional: in-plane component if you later model axial/membrane effects
    w_inpl = lift_per_unit_span * sin_sweep    # N/m (currently unused)
    
    # Torsional per-span moment from lift at chordwise lever arm, projected by cos(Λ)
    m_t_proj = m_t * cos_sweep                 # N·m per m
    
    # Initialize resultants (cantilever integration tip -> root)
    V = np.zeros_like(xspan)   # flapwise shear
    M = np.zeros_like(xspan)   # flapwise bending moment
    T = np.zeros_like(xspan)   # torsional moment
    
    print("\n-- sweep-projected loads (first 12 stations) --")
    for ii in range(0, min(12, len(xspan))):
        print(f"  i={ii:02d}, x={xspan[ii]:.4f} m, w_flap={w_flap[ii]:.3g} N/m, m_t_proj={m_t_proj[ii]:.3g} N·m/m")
    print("-- end sample --\n")
    
    # Integrate (trapezoidal) from tip to root
    for i in range(n-2, -1, -1):
        dxn = xspan[i+1] - xspan[i]
        # Shear from flapwise load
        V[i] = V[i+1] + 0.5 * (w_flap[i] + w_flap[i+1]) * dxn
        # Bending from shear
        M[i] = M[i+1] + 0.5 * (V[i] + V[i+1]) * dxn
        # Torsion from projected per-span torque
        T[i] = T[i+1] + 0.5 * (m_t_proj[i] + m_t_proj[i+1]) * dxn
    
    # Sign tidy (same convention as before)
    if M[0] < 0:
        V = -V; M = -M; T = -T


    # Bending stresses - use CLT-derived EI scaling for skins
    # Compute nominal linear bending stress as sigma = M * z / EI_total, and scale by effective material modulus
    # E_ref is defined earlier (right after sandwich_bending_properties)

    # Remove provisional rod stress calculation here. We'll compute final skin and rod stresses after curvature (kappa) is available.
    sigma_skin_top = np.zeros_like(M)  # placeholder array to be filled after kappa
    sigma_skin_bot = np.zeros_like(M)
    sigma_rod_top = np.zeros_like(M)
    sigma_rod_bot = np.zeros_like(M)

    # Multicell setup: split polygon by xcuts
    cell_polys = split_polygon_by_xcuts(xy_poly, [x_cut1, x_cut2])

    def cell_areas_from_upper_lower(upper_arr, lower_arr, c, p1_local, p2_local):
        xvals = upper_arr[:, 0]
        t = upper_arr[:, 1] - lower_arr[:, 1]
        try:
            A_total = np.trapezoid(t, xvals)
            A1 = np.trapezoid(t[xvals <= p1_local * c], xvals[xvals <= p1_local * c]) if np.any(xvals <= p1_local * c) else 0.0
            sel2 = (xvals > p1_local * c) & (xvals <= p2_local * c)
            A2 = np.trapezoid(t[sel2], xvals[sel2]) if np.any(sel2) else 0.0
        except AttributeError:
            A_total = np.trapz(t, xvals)
            A1 = np.trapz(t[xvals <= p1_local * c], xvals[xvals <= p1_local * c]) if np.any(xvals <= p1_local * c) else 0.0
            sel2 = (xvals > p1_local * c) & (xvals <= p2_local * c)
            A2 = np.trapz(t[sel2], xvals[sel2]) if np.any(sel2) else 0.0
        A3 = A_total - A1 - A2
        return [max(A1, 1e-12), max(A2, 1e-12), max(A3, 1e-12)]

    if len(cell_polys) < 3:
        A_cells = cell_areas_from_upper_lower(upper, lower, chord_m, p1, p2)
        perim = np.sum(np.sqrt(np.sum(np.diff(raw_outline, axis=0) ** 2, axis=1)))
        L_top = perim * 0.33
        L_mid = perim * 0.34
        L_bot = perim * 0.33
        # web vertical separation should equal rod centroid vertical separation (top-bottom)
        rod_vert_sep = max(1e-6, abs(sb.get('z_rod_top', 0.0) - sb.get('z_rod_bot', 0.0)))
        walls = [
            {'length': L_top, 'thickness': t_face, 'cell_ids': [0, 1]},
            {'length': rod_vert_sep, 'thickness': t_web, 'cell_ids': [1]},
            {'length': L_mid, 'thickness': t_face, 'cell_ids': [1, 2]},
            {'length': rod_vert_sep, 'thickness': t_web, 'cell_ids': [2]},
            {'length': L_bot, 'thickness': t_face, 'cell_ids': [0, 2]},
        ]
    else:
        cell_polys_sh = [Polygon(poly) for poly in cell_polys]
        A_cells = [max(float(p.area), 1e-12) for p in cell_polys_sh]
        perims = [max(float(p.length), 1e-12) for p in cell_polys_sh]
        vertical_lines = [
            LineString([(x_cut1, -1e3), (x_cut1, 1e3)]),
            LineString([(x_cut2, -1e3), (x_cut2, 1e3)])
        ]
        web_lengths = []
        web_cell_pairs = []
        poly_all = Polygon(xy_poly)
        for xl in vertical_lines:
            try:
                inter = xl.intersection(poly_all)
                length = float(inter.length) if not inter.is_empty else 0.0
            except Exception:
                length = 0.0
            touching = [i for i, p in enumerate(cell_polys_sh) if p.intersects(xl)]
            if not touching:
                xs_centers = [float(p.centroid.x) for p in cell_polys_sh]
                nearest = int(np.argmin([abs(xc - xl.centroid.x) for xc in xs_centers]))
                touching = [nearest]
            web_lengths.append(max(length, 1e-12))
            web_cell_pairs.append(touching)

        walls_temp = []
        Ncells = len(cell_polys_sh)
        if Ncells >= 2:
            L01 = 0.25 * (perims[0] + perims[1])
            walls_temp.append({'length': max(1e-6, L01), 'thickness': t_face, 'cell_ids': [0, 1]})
        else:
            walls_temp.append({'length': max(1e-6, perims[0]), 'thickness': t_face, 'cell_ids': [0]})

        # use rod vertical separation for web lengths (rods run spanwise and connect top/bottom)
        rod_vert_sep = max(1e-6, abs(sb.get('z_rod_top', 0.0) - sb.get('z_rod_bot', 0.0)))
        walls_temp.append({'length': rod_vert_sep, 'thickness': t_web, 'cell_ids': web_cell_pairs[0]})

        if Ncells >= 3:
            L12 = 0.25 * (perims[1] + perims[2])
            walls_temp.append({'length': max(1e-6, L12), 'thickness': t_face, 'cell_ids': [1, 2]})
        else:
            idx = 1 if Ncells > 1 else 0
            walls_temp.append({'length': max(1e-6, perims[idx]), 'thickness': t_face, 'cell_ids': [idx]})

        if len(web_lengths) > 1:
            walls_temp.append({'length': rod_vert_sep, 'thickness': t_web, 'cell_ids': web_cell_pairs[1]})
        else:
            walls_temp.append({'length': rod_vert_sep, 'thickness': t_web, 'cell_ids': [min(Ncells-1, 1)]})

        if Ncells >= 3:
            L0N = 0.25 * (perims[0] + perims[-1])
            walls_temp.append({'length': max(1e-6, L0N), 'thickness': t_face, 'cell_ids': [0, Ncells-1]})
        else:
            walls_temp.append({'length': max(1e-6, perims[0]), 'thickness': t_face, 'cell_ids': [0]})

        walls = walls_temp

    def _validate_and_fix_walls(A_cells_list, walls_list):
        nc_local = len(A_cells_list)
        fixed = []
        for w in walls_list:
            raw_ids = [int(ci) for ci in w.get('cell_ids', [])]
            ids_valid = [i for i in raw_ids if 0 <= i < nc_local]
            if not ids_valid:
                ids_valid = [max(0, min(nc_local-1, nc_local // 2))]
            seen = set()
            ids_clean = []
            for ii in ids_valid:
                if ii not in seen:
                    ids_clean.append(ii)
                    seen.add(ii)
            fixed.append({
                'length': float(w.get('length', 1e-6)),
                'thickness': float(w.get('thickness', 1e-6)),
                'cell_ids': ids_clean
            })
        return fixed

    walls = _validate_and_fix_walls(A_cells, walls)

    print(f"\n=== TORSION DIMENSIONAL CHECK ===")
    print(f"chord_m = {chord_m:.6f} m")
    print(f"A_cells = {A_cells}")
    print(f"A_cells[0] = {A_cells[0]:.6e} m² (physical)")
    print(f"A_cells[0]/(chord_m**2) = {A_cells[0]/(chord_m**2):.4f} (normalized check)")
    print(f"Expected normalized value: ~0.01-0.15 for typical airfoil")
    print("====================================\n")

    # Assign G per wall using CLT where available: skin uses sb['G_skin_eff'], web uses web CLT
    # skin membrane shear modulus from sandwich computation
    G_skin_eff = sb.get('G_skin_eff', max(mat.E1_skin/(2*(1+mat.nu12_skin)), 1e6))
    
    # build web CLT using explicit web laminate properties and web ply thickness
    ply_props_web = {'E1': mat.E1_web, 'E2': mat.E2_web, 'G12': mat.G12_web, 'nu12': mat.nu12_web}
    A_web, B_web, D_web, t_web_total = laminate_ABD_from_stack(mat.web_stack, ply_props_web, mat.t_ply_web)
    G_web_eff = laminate_membrane_G(A_web, t_web_total)  # Use A_web and membrane function!
    G_vals = []
    for w in walls:
        if abs(w['thickness'] - mat.t_web_total) < 1e-9:
            G_vals.append(max(G_web_eff, 1e6))
        else:
            G_vals.append(max(G_skin_eff, 1e6))

    # Check if we need to convert to single-cell topology
    if getattr(inp, '_no_webs', False):
        print("\n=== NO-WEBS MODE: Converting to single-cell topology ===")
        print(f"Original: {len(A_cells)} cells, {len(walls)} walls")
        A_cells_torsion, walls_torsion, G_vals_torsion = convert_to_single_cell_topology(A_cells, walls, G_vals)
        print(f"Single-cell: {len(A_cells_torsion)} cell, {len(walls_torsion)} walls")
        print(f"A_single = {A_cells_torsion[0]:.6e} m² (merged)")
    else:
        A_cells_torsion = A_cells
        walls_torsion = walls
        G_vals_torsion = G_vals

    # Torsion analysis
    tau_skin_tor = np.zeros_like(T)
    tau_web1 = np.zeros_like(T)
    tau_web2 = np.zeros_like(T)

    print("DEBUG: moment_arm (m) =", moment_arm)
    print("DEBUG: sample m_t[0..3] =", m_t[:4])
    print("DEBUG: sample T[0..3] =", T[:4])

    for i in range(n):
        mc = multicell_bredt_solver(A_cells_torsion, walls_torsion, G_vals_torsion, T[i])
        tw = mc.get('tau_wall', np.zeros(len(walls_torsion)))
        tau_skin_tor[i] = max(
            tw[0] if len(tw) > 0 else 0.0,
            tw[2] if len(tw) > 2 else 0.0,
            tw[-1] if len(tw) > 0 else 0.0
        )
        # In no-webs mode, set torsional web stresses to zero
        if getattr(inp, '_no_webs', False):
            tau_web1[i] = 0.0
            tau_web2[i] = 0.0
        else:
            tau_web1[i] = tw[1] if len(tw) > 1 else 0.0
            tau_web2[i] = tw[3] if len(tw) > 3 else 0.0

    # Web shear stress calculation
    # NOTE: Q and I_trans computed for optional section property diagnostics only.
    # They are NOT used for web transverse shear stress (which uses GA-based direct stress).
    # VQ/It would be needed only if implementing full thin-wall shear flow distribution
    # under transverse shear V (not currently implemented - GA partition is the chosen model).
    tau_web1_arr = np.zeros_like(V)
    tau_web2_arr = np.zeros_like(V)
    I_trans_arr = np.zeros_like(V)
    Q1_arr = np.zeros_like(V)
    Q2_arr = np.zeros_like(V)
    
    # Web geometry - compute vertical span between internal attachment points (rod centroids / inner walls)
    # Prefer the distance between spar-rod centroids if rods are inside the foam as modeled in sandwich_bending_properties.
    try:
        zrod_top = float(sb.get('z_rod_top', 0.0))
        zrod_bot = float(sb.get('z_rod_bot', 0.0))
        # distance between rod centroids (attachment points)
        h_web_effective = max(1e-6, abs(zrod_top - zrod_bot))
        # sanity: if this is unphysically small or greater than core, clamp to core thickness
        if h_web_effective < 1e-6 or h_web_effective > max(0.0, t_core * 1.5):
            h_web_effective = max(1e-6, t_core)
            print(f"DEBUG: computed rod-to-rod web height out of bounds, clamped to t_core = {t_core:.6g} m")
        else:
            print(f"DEBUG: web effective height from rod centroids = {h_web_effective:.6g} m (zrod_top={zrod_top:.6g}, zrod_bot={zrod_bot:.6g})")
    except Exception:
        h_web_effective = max(1e-6, t_core)
        print("DEBUG: failed to compute rod-based web height, using t_core as fallback")

    t_web_actual = max(mat.t_web_total, 1e-6)
    
    # Transformed section properties (E_ref defined earlier)
    A_skin = chord_m * t_face
    A_rod = inp.rod_width_m * inp.rod_height_m
    A_skin_trans = A_skin * (mat.E1_fiber / E_ref)
    A_rod_trans = A_rod * (mat.E_rod / E_ref)
    I_rod_self = (inp.rod_width_m * inp.rod_height_m ** 3) / 12.0
    
    # Rod positions (from sandwich computation). Keep gap for diagnostics only.
    z_skin_top = sb['z_skin_top']
    z_skin_bot = sb['z_skin_bot']
    gap = inp.rod_gap_mm * 1e-3   # m, used by sandwich_bending_properties; kept here for diagnostics if needed
    
    # Use the correct distances already computed in sandwich_bending_properties()
    zrod_top = sb['z_rod_top']
    zrod_bot = sb['z_rod_bot']

    
    # Web neutral axis positions (at foam midplane since they span full foam height)
    z_web1 = 0.0
    z_web2 = 0.0
    
    # Rod distribution (use explicitly defined positions from geometry setup)
    # These were computed earlier based on actual rod x-positions
    n_rods_top = sb.get('n_rods_top', 2)
    n_rods_bot = sb.get('n_rods_bot', 2)
    # Use geometrically-derived rod counts (already computed above)
    n_top_p1 = n_top_fwd_1  # from explicit rod positions
    n_bot_p1 = n_bot_fwd_1
    n_top_p2 = n_top_fwd_2
    n_bot_p2 = n_bot_fwd_2
    
    # Check if tau_web_shear was calculated correctly
    tau_web_shear_check = 0.5 * (tau_web1_arr + tau_web2_arr)
    print(f"\n=== WEB SHEAR GEOMETRY ===")
    print(f"Web effective height: {h_web_effective*1000:.3f} mm")
    print(f"Web thickness: {t_web_actual*1000:.3f} mm")
    print(f"t_web_actual: {t_web_actual:.6g} m")
    print(f"mat.t_web_total: {mat.t_web_total:.6g} m")
    print(f"Web 1 position: {p1*100:.1f}% chord = {p1*chord_m*1000:.2f} mm")
    print(f"Web 2 position: {p2*100:.1f}% chord = {p2*chord_m*1000:.2f} mm")
    print(f"Web positions: z_web1={z_web1:.6g} m, z_web2={z_web2:.6g} m")
    print(f"Skin positions: z_skin_top={sb['z_skin_top']:.6g} m, z_skin_bot={sb['z_skin_bot']:.6g} m")
    print(f"Rod positions (top): {zrod_top:.6g} m, (bot): {zrod_bot:.6g} m")
    print(f"Rod distribution: top({n_top_p1}|{n_top_p2}), bot({n_bot_p1}|{n_bot_p2})")
    print("===========================\n")
    
    #START PATCH: precise Q and I_trans using existing polygon
    pts_poly = np.asarray(xy_poly)  # closed polygon as used earlier, exact geometry preserved
    # ensure array shape (N,2) and closed
    if pts_poly.shape[0] > 1 and not np.allclose(pts_poly[0], pts_poly[-1], atol=1e-12):
        pts_poly = np.vstack([pts_poly, pts_poly[0]])
    
    def forward_skin_segment_stats(pts, x_web, top=True):
        # Walk polygon in order and collect contiguous segments that are:
        #   - on the requested surface sign (top: y>=0, bottom: y<=0)
        #   - with x <= x_web (forward of web) using polygon x-coordinates
        # This preserves geometry and uses exact points only (no geometry modification).
        L = 0.0
        z_times_len = 0.0
        I_self = 0.0
        # iterate edges in polygon order
        for k in range(len(pts)-1):
            p0 = pts[k]; p1 = pts[k+1]
            # robust forward selection: include the segment if any endpoint is forward of x_web
            # or if the segment spans the web x (covers curved segments reliably)
            x0, x1 = p0[0], p1[0]
            y0, y1 = p0[1], p1[1]
            # decide if the segment belongs to requested surface by sign of midpoint y (keeps top/bottom test)
            ym = 0.5*(y0 + y1)
            on_top = (ym >= 0.0)
            on_bottom = (ym <= 0.0)
            eps = 1e-12
            is_forward = (x0 <= x_web + eps) or (x1 <= x_web + eps) or (min(x0, x1) <= x_web + eps and max(x0, x1) >= x_web - eps)
            pick = (is_forward and ((top and on_top) or (not top and on_bottom)))
            if pick:
                seg_len = math.hypot(x1 - x0, y1 - y0)
                if seg_len <= 0:
                    continue
                L += seg_len
                z_times_len += 0.5*(p0[1] + p1[1]) * seg_len
                # thin-strip self I about centroid approx: (length * t^3)/12
                I_self += (seg_len * (t_face**3)) / 12.0
        if L <= 0:
            return 0.0, 0.0, 0.0
        centroid_z = z_times_len / L
        return L, centroid_z, I_self
    
    # precompute rod area and self I (use variables already defined)
    A_rod = inp.rod_width_m * inp.rod_height_m
    I_rod_self = (inp.rod_width_m * inp.rod_height_m**3) / 12.0
    # E_ref defined earlier (after sandwich_bending_properties)
    
    # Rod counts forward of each web already computed correctly at lines 1447-1450
    # (n_top_p1, n_bot_p1, n_top_p2, n_bot_p2 contain absolute counts, not incremental)
    # Use those values directly - DO NOT redefine to avoid double-counting
    n_top_fwd_1 = n_top_p1
    n_bot_fwd_1 = n_bot_p1
    n_top_fwd_2 = n_top_p2  # FIXED: was n_top_p1 + n_top_p2 (double-count bug)
    n_bot_fwd_2 = n_bot_p2  # FIXED: was n_bot_p1 + n_bot_p2 (double-count bug)
    
    # E_ref is defined once at the top (after sandwich_bending_properties)
    # and used consistently for all transformed-section calculations
    E_face_eff = E_ref  # consistent reference

    # Material stiffness ratios (use these for transformed-area calculations)
    n_skin = E_face_eff / max(E_ref, 1e-18)  # = 1.0
    n_rod  = float(mat.E_rod) / max(E_ref, 1e-18)

    # Skin total area per face (geometric)
    A_skin_total = chord_m * t_face
    
    # Transformed areas (use the canonical E_ref through n_skin/n_rod)
    A_skin_trans = n_skin * A_skin_total
    A_rod_geom   = inp.rod_width_m * inp.rod_height_m
    A_rod_trans  = n_rod * A_rod_geom
    I_rod_self   = (inp.rod_width_m * inp.rod_height_m**3) / 12.0  # geometric self-I
    
    # Precompute forward skin segment stats once (they don't vary with i)
    Lt1, zt1, Iself_t1 = forward_skin_segment_stats(pts_poly, x_cut1, top=True)
    Lb1, zb1, Iself_b1 = forward_skin_segment_stats(pts_poly, x_cut1, top=False)
    Lt2, zt2, Iself_t2 = forward_skin_segment_stats(pts_poly, x_cut2, top=True)
    Lb2, zb2, Iself_b2 = forward_skin_segment_stats(pts_poly, x_cut2, top=False)
    
    # Transform the thin-strip self-I from skins by n_skin as well
    Iself_t1_trans = n_skin * Iself_t1
    Iself_b1_trans = n_skin * Iself_b1
    Iself_t2_trans = n_skin * Iself_t2
    Iself_b2_trans = n_skin * Iself_b2
    
    # Counts forward of webs (already computed above)
    # n_top_fwd_1, n_bot_fwd_1, n_top_fwd_2, n_bot_fwd_2
    
    # ========================================================================
    # TRANSFORMED SECTION PROPERTIES (NEUTRAL AXIS AND I_TOTAL)
    # Compute once for the section (geometry-independent of station loading)
    # ========================================================================
    
    # --- Step 1: Find true neutral axis z_NA using transformed areas ---
    # Reference: z_NA = Σ(n_i * A_i * z_i) / Σ(n_i * A_i)
    
    sum_EA_z = 0.0  # numerator: sum of E*A*z products
    sum_EA = 0.0    # denominator: sum of E*A
    
    # Material stiffness ratios (already defined earlier)
    # n_skin = E_face_eff / E_ref  (should be ~1.0)
    # n_rod = mat.E_rod / E_ref
    
    # Integrate skins over chord using fine grid
    nint_NA = 401
    x_grid_NA = np.linspace(0.0, chord_m, nint_NA)
    yu_NA = np.interp(x_grid_NA, upper[:, 0], upper[:, 1])
    yl_NA = np.interp(x_grid_NA, lower[:, 0], lower[:, 1])
    
    # Skin contribution (top and bottom faces)
    for j in range(len(x_grid_NA) - 1):
        dx_j = x_grid_NA[j+1] - x_grid_NA[j]
        # Top skin strip (centroid is offset outward from surface by t_face/2)
        z_surface_top = (yu_NA[j] + yu_NA[j+1]) / 2.0
        z_top = z_surface_top + t_face / 2.0
        dA_top = t_face * dx_j
        sum_EA_z += n_skin * dA_top * z_top
        sum_EA += n_skin * dA_top
        # Bottom skin strip (centroid is offset outward from surface by t_face/2)
        z_surface_bot = (yl_NA[j] + yl_NA[j+1]) / 2.0
        z_bot = z_surface_bot - t_face / 2.0
        dA_bot = t_face * dx_j
        sum_EA_z += n_skin * dA_bot * z_bot
        sum_EA += n_skin * dA_bot
    
    # Rod contribution (all 4 rods: 2 top, 2 bottom)
    sum_EA_z += n_rods_top * (n_rod * A_rod_geom * zrod_top)
    sum_EA += n_rods_top * (n_rod * A_rod_geom)
    sum_EA_z += n_rods_bot * (n_rod * A_rod_geom * zrod_bot)
    sum_EA += n_rods_bot * (n_rod * A_rod_geom)
    
    # Core contribution (usually negligible, but include for completeness)
    A_core_total = chord_m * t_core
    n_core = float(mat.E_core) / max(E_ref, 1e-18)
    sum_EA_z += n_core * A_core_total * 0.0  # core at z=0
    sum_EA += n_core * A_core_total
    
    # Neutral axis position
    z_NA = sum_EA_z / max(sum_EA, 1e-18)
    
    print(f"\n=== TRANSFORMED SECTION PROPERTIES ===")
    print(f"Neutral axis position: z_NA = {z_NA*1000:.4f} mm (from geometric midplane)")
    print(f"  (z=0 is foam midplane; positive = above midplane)")
    
    # --- Step 2: Compute total section I about z_NA (parallel axis theorem) ---
    I_total_trans = 0.0
    
    # Skin contribution (integrate I about z_NA)
    for j in range(len(x_grid_NA) - 1):
        dx_j = x_grid_NA[j+1] - x_grid_NA[j]
        # Top skin (centroid is offset outward from surface by t_face/2)
        z_surface_top = (yu_NA[j] + yu_NA[j+1]) / 2.0
        z_top = z_surface_top + t_face / 2.0
        dA_top = t_face * dx_j
        # Self-inertia of thin strip about its own centroid (perpendicular-axis theorem)
        dI_self_top = (dx_j * t_face**3) / 12.0
        I_total_trans += n_skin * (dI_self_top + dA_top * (z_top - z_NA)**2)
        
        # Bottom skin (centroid is offset outward from surface by t_face/2)
        z_surface_bot = (yl_NA[j] + yl_NA[j+1]) / 2.0
        z_bot = z_surface_bot - t_face / 2.0
        dA_bot = t_face * dx_j
        dI_self_bot = (dx_j * t_face**3) / 12.0
        I_total_trans += n_skin * (dI_self_bot + dA_bot * (z_bot - z_NA)**2)
    
    # Rod contribution (parallel axis: I_self + A*d^2)
    I_rod_self_geom = (inp.rod_width_m * inp.rod_height_m**3) / 12.0
    I_total_trans += n_rods_top * n_rod * (I_rod_self_geom + A_rod_geom * (zrod_top - z_NA)**2)
    I_total_trans += n_rods_bot * n_rod * (I_rod_self_geom + A_rod_geom * (zrod_bot - z_NA)**2)
    
    # Core contribution (usually small)
    I_core_self = (chord_m * t_core**3) / 12.0
    I_total_trans += n_core * (I_core_self + A_core_total * (0.0 - z_NA)**2)
    
    # Safety floor
    I_total_trans = max(I_total_trans, 1e-10)
    
    print(f"Total transformed I: I_total_trans = {I_total_trans:.6e} m^4")
    print(f"  (For comparison: EI_total from sandwich = {EI_total:.6e} N·m^2)")
    print(f"  Implied E_avg = {EI_total / I_total_trans / 1e9:.2f} GPa")
    print("=====================================\n")
    
    for i in range(n):
        # ---- Continuous integration for Q and physically based I_trans for both webs ----
        # Fine x-grid over chord for integration (meters)
        nint_local = 401
        x_grid_local = np.linspace(0.0, chord_m, nint_local)
        yu_i = np.interp(x_grid_local, upper[:, 0], upper[:, 1])
        yl_i = np.interp(x_grid_local, lower[:, 0], lower[:, 1])

        # transformed skin strip area per unit x (m^2 per m)
        # n_skin should already be set earlier; if not, fall back to a safe ratio
        n_skin = n_skin if 'n_skin' in locals() else (E_face_eff / max(E_ref, 1e-18))
        A_skin_strip = n_skin * t_face

        # web x positions in meters (use your existing x_cut1/x_cut2)
        xw1 = float(x_cut1)
        xw2 = float(x_cut2)
        z_web = 0.0  # neutral axis reference (midplane)

        # Masks for forward-of-web samples
        mask_fwd1 = x_grid_local <= (xw1 + 1e-12)
        mask_fwd2 = x_grid_local <= (xw2 + 1e-12)

        # ---- Q calculation: first moment of area FORWARD of each web, about z_NA ----
        # Use the consistent neutral axis z_NA (computed once above)
        
        # --- Web 1: Q = first moment of area forward of xw1, about z_NA ---
        if mask_fwd1.sum() >= 2:
            x_f1 = x_grid_local[mask_fwd1]
            yu_f1 = yu_i[mask_fwd1]
            yl_f1 = yl_i[mask_fwd1]
            # Q = ∫ A_strip * (z_centroid - z_NA) dx for top and bottom
            # Skin centroid is offset from surface by t_face/2
            z_top_centroid_1 = yu_f1 + t_face / 2.0
            z_bot_centroid_1 = yl_f1 - t_face / 2.0
            try:
                Q_skin_top_1 = np.trapezoid(A_skin_strip * (z_top_centroid_1 - z_NA), x_f1)
                Q_skin_bot_1 = np.trapezoid(A_skin_strip * (z_bot_centroid_1 - z_NA), x_f1)
            except AttributeError:
                Q_skin_top_1 = np.trapz(A_skin_strip * (z_top_centroid_1 - z_NA), x_f1)
                Q_skin_bot_1 = np.trapz(A_skin_strip * (z_bot_centroid_1 - z_NA), x_f1)
        else:
            Q_skin_top_1 = 0.0
            Q_skin_bot_1 = 0.0
        
        # --- Web 2: Q = first moment of area forward of xw2, about z_NA ---
        if mask_fwd2.sum() >= 2:
            x_f2 = x_grid_local[mask_fwd2]
            yu_f2 = yu_i[mask_fwd2]
            yl_f2 = yl_i[mask_fwd2]
            # Skin centroid is offset from surface by t_face/2
            z_top_centroid_2 = yu_f2 + t_face / 2.0
            z_bot_centroid_2 = yl_f2 - t_face / 2.0
            try:
                Q_skin_top_2 = np.trapezoid(A_skin_strip * (z_top_centroid_2 - z_NA), x_f2)
                Q_skin_bot_2 = np.trapezoid(A_skin_strip * (z_bot_centroid_2 - z_NA), x_f2)
            except AttributeError:
                Q_skin_top_2 = np.trapz(A_skin_strip * (z_top_centroid_2 - z_NA), x_f2)
                Q_skin_bot_2 = np.trapz(A_skin_strip * (z_bot_centroid_2 - z_NA), x_f2)
        else:
            Q_skin_top_2 = 0.0
            Q_skin_bot_2 = 0.0
        
        # Rods forward of each web (transformed area * distance from NA * count)
        Q_rods_1 = (n_top_fwd_1 * (A_rod_trans * (zrod_top - z_NA))
                   + n_bot_fwd_1 * (A_rod_trans * (zrod_bot - z_NA)))
        Q_rods_2 = (n_top_fwd_2 * (A_rod_trans * (zrod_top - z_NA))
                   + n_bot_fwd_2 * (A_rod_trans * (zrod_bot - z_NA)))
        
        # Sum contributions (signed, then take absolute value)
        Q1_val = (Q_skin_top_1 + Q_skin_bot_1) + Q_rods_1
        Q2_val = (Q_skin_top_2 + Q_skin_bot_2) + Q_rods_2
        
        # Store absolute Q values (shear stress magnitude)
        Q1_arr[i] = abs(float(Q1_val))
        Q2_arr[i] = abs(float(Q2_val))
        
        # Store I_total for diagnostics (same at all stations)
        I_trans_arr[i] = I_total_trans

    # Web shear stress computation moved to AFTER GA partition for consistency
    # (will use V_web instead of full V)
    
    # --- Core shear: GA stiffness partition (physically consistent) ---
    # Partition shear between webs and core using parallel spring model.
    # This avoids circular logic of extracting forces from VQ/It stress fields
    # computed with full V. VQ/It kept only for diagnostic comparison.
    
    # Core shear effective width (85% of chord for LE/TE ineffectiveness)
    b_eff = 0.85 * chord_m
    t_core = float(mat.t_core)
    G_core_val = float(mat.G_core)
    webs_active = not getattr(inp, '_no_webs', False)
    
    # Web geometry
    n_webs = 2  # web1 and web2
    h_web_avg = h_web_effective  # rod-to-rod distance (already computed)
    
    # Use CLT-derived web shear modulus (already computed in torsion section)
    # G_web_eff is in scope from line ~1403
    
    print(f"DEBUG: Core shear b_eff = {b_eff*1000:.2f} mm")
    print(f"DEBUG: Shear stiffness - GA_web = {n_webs} × {G_web_eff/1e9:.1f} GPa × {t_web_actual*1e3:.3f} mm × {h_web_avg*1e3:.2f} mm")
    print(f"DEBUG: Shear stiffness - GA_core = {G_core_val/1e6:.1f} MPa × {t_core*1e3:.2f} mm × {b_eff*1e3:.2f} mm")
    
    # Partition shear at each station using stiffness ratios
    V_web_total = np.zeros_like(V)
    V_core_arr = np.zeros_like(V)
    
    for i in range(n):
        V_web_total[i], V_core_arr[i] = shear_partition_GA(
            V=V[i],
            n_webs=n_webs,
            h_web=h_web_avg,
            t_web=t_web_actual,
            G_web=G_web_eff,
            b_eff=b_eff,
            t_core=t_core,
            G_core=G_core_val,
            webs_active=webs_active
        )
    
    # ---- Compute web shear stress directly from partitioned force (consistent with GA) ----
    tau_web1_arr = np.zeros_like(V)
    tau_web2_arr = np.zeros_like(V)
    
    # Use already-computed effective web height (h_web_effective from lines 1502-1517)
    # This is the rod-to-rod vertical distance (physical attachment points)
    h_web1_eff = h_web_effective
    h_web2_eff = h_web_effective
    
    # Web shear stiffnesses for stiffness-weighted split
    GA_web1 = G_web_eff * t_web_actual * h_web1_eff
    GA_web2 = G_web_eff * t_web_actual * h_web2_eff
    GA_web_total = GA_web1 + GA_web2
    
    for i in range(n):
        if not webs_active or V_web_total[i] <= 1e-12:
            tau_web1_arr[i] = 0.0
            tau_web2_arr[i] = 0.0
        else:
            # Stiffness-weighted split (robust for asymmetric webs)
            V_web1 = V_web_total[i] * (GA_web1 / GA_web_total) if GA_web_total > 0 else V_web_total[i] / 2.0
            V_web2 = V_web_total[i] * (GA_web2 / GA_web_total) if GA_web_total > 0 else V_web_total[i] / 2.0
            
            # Direct stress: tau = V_web / A_web (parallel spring model)
            tau_web1_arr[i] = abs(V_web1) / (h_web1_eff * t_web_actual)
            tau_web2_arr[i] = abs(V_web2) / (h_web2_eff * t_web_actual)
    
    # Update tau_web_shear (max of web1/web2 for conservative check)
    if webs_active:
        tau_web_shear = np.maximum(tau_web1_arr, tau_web2_arr)
    else:
        tau_web_shear = np.zeros_like(tau_web1_arr)
    
    # Core shear stress from partitioned load
    A_core_shear = b_eff * t_core
    if A_core_shear <= 0:
        tau_core = np.zeros_like(V)
    else:
        tau_core = V_core_arr / A_core_shear  # Pa
    
    # Skin membrane shear: NOT physical in sandwich theory
    tau_skin_membrane = np.zeros_like(V)  # Zero for sandwich model
    
    print(f"DEBUG: Load partition - V_max={np.nanmax(np.abs(V)):.3g} N")
    print(f"DEBUG:   V_web_total={np.nanmax(V_web_total):.3g} N ({100*np.nanmax(V_web_total)/max(np.nanmax(np.abs(V)),1e-12):.1f}%)")
    print(f"DEBUG:   V_core={np.nanmax(V_core_arr):.3g} N ({100*np.nanmax(V_core_arr)/max(np.nanmax(np.abs(V)),1e-12):.1f}%)")
    print(f"DEBUG: tau_core_max={np.nanmax(tau_core):.3g} Pa (from GA partition)")
    print(f"DEBUG: tau_web_max={np.nanmax(tau_web_shear):.3g} Pa (direct stress from GA partition)")
    print(f"DEBUG: b_eff={b_eff*1000:.2f} mm, t_core={t_core*1000:.2f} mm, A_core_shear={A_core_shear*1e6:.3f} mm²")
    
    """
    ===== DIAGNOSTIC ONLY: VQ/It with V_web (INCONSISTENT, DO NOT USE) =====
    This block shows what the old mixed model produced. It is WRONG because:
    1. Mixes GA partition (V_web) with beam theory (VQ/It)  
    2. Uses V_web from parallel springs in a formula derived for different physics
    3. Result is neither pure GA nor pure beam shear flow distribution
    Kept commented for educational/comparison purposes only.
    
    # Diagnostic: compute what full-V would give (for comparison only)
    tau_web1_arr_full_V = np.zeros_like(V)
    tau_web2_arr_full_V = np.zeros_like(V)
    
    for i in range(n):
        t_use = max(t_web_actual, 1e-6)
        
        # Diagnostic with full V (old inconsistent approach, for comparison)
        if webs_active and abs(V[i]) > 1e-12:
            tau_web1_arr_full_V[i] = abs(V[i]) * Q1_arr[i] / (I_total_trans * t_use)
            tau_web2_arr_full_V[i] = abs(V[i]) * Q2_arr[i] / (I_total_trans * t_use)
    
    if webs_active and np.nanmax(tau_web1_arr_full_V) > 0:
        print(f"DEBUG: [Diagnostic] tau_web with full V would be: {np.nanmax(np.maximum(tau_web1_arr_full_V, tau_web2_arr_full_V)):.3g} Pa (inconsistent, for comparison only)")
    """

    # DEBUG: Check web shear calculation
    print("\n=== WEB SHEAR DEBUG (GA-BASED DIRECT STRESS) ===")
    print(f"Web stiffness: GA_web1 = {GA_web1:.3g} N, GA_web2 = {GA_web2:.3g} N")
    print(f"Web geometry: h_web_eff = {h_web1_eff*1000:.3f} mm, t_web = {t_web_actual*1000:.3f} mm")
    print(f"tau_web1_arr: min={np.nanmin(tau_web1_arr):.3g}, max={np.nanmax(tau_web1_arr):.3g}, mean={np.nanmean(tau_web1_arr):.3g}")
    print(f"tau_web2_arr: min={np.nanmin(tau_web2_arr):.3g}, max={np.nanmax(tau_web2_arr):.3g}, mean={np.nanmean(tau_web2_arr):.3g}")
    print(f"V_web_total: min={np.nanmin(V_web_total):.3g}, max={np.nanmax(V_web_total):.3g}")
    print(f"V (shear force): min={np.nanmin(np.abs(V)):.3g}, max={np.nanmax(np.abs(V)):.3g}")
    print("==========================================\n")
    
    # ---- Combine transverse + torsion shear at component level ----
    # Webs experience both transverse shear (tau_web1/2_arr) and torsion (tau_web1/2 from Bredt)
    # Both are arrays over span (length n)
    # Combine using von Mises-like RSS (conservative for shear):
    tau_web1_total = np.sqrt(tau_web1_arr**2 + tau_web1**2)  # RSS combination, per station
    tau_web2_total = np.sqrt(tau_web2_arr**2 + tau_web2**2)  # RSS combination, per station
    tau_web_shear_total = np.maximum(tau_web1_total, tau_web2_total)
    
    # Skins: torsion only (transverse shear carried by webs in GA model)
    # tau_skin_tor already computed (line ~1471)
    
    # Core: transverse shear only (foam doesn't resist torsion significantly)
    # tau_core already computed (line ~1910)
    
    # FoS calculations
    # Bending stresses from curvature (physically consistent)
    # Compute curvature (vector-safe for M as array over span)
    kappa = np.asarray(M) / np.maximum(EI_total, 1e-18)   # 1/m

    # Extract in-plane axial modulus from CLT (A11/t, not D-derived E)
    A_face = sb.get('A_face')
    t_skin = float(sb.get('t_skin_total', mat.t_skin_total))

    if A_face is not None and t_skin > 0:
        E_face_axial = float(laminate_inplane_axial_modulus(A_face, t_skin))
    else:
        E_face_axial = float(mat.E1_skin)  # Fallback

    # Face stress from shared curvature strain field: eps = -kappa*z
    # Sign chosen to match existing M convention (top compression when M>0)
    z_top = float(sb['z_skin_top'])
    z_bot = float(sb['z_skin_bot'])

    sigma_skin_top = -E_face_axial * kappa * z_top
    sigma_skin_bot = -E_face_axial * kappa * z_bot

    # Rods: stress = E_rod * curvature * z_rod (physical)
    # SIGN CONVENTION: Negative sign for consistency with skin stress
    sigma_rod_top = -float(mat.E_rod) * kappa * sb['z_rod_top']
    sigma_rod_bot = -float(mat.E_rod) * kappa * sb['z_rod_bot']

    # Quick diagnostics (insert immediately after sigma_rod_top/sigma_rod_bot are computed via kappa)
    print("DIAG: kappa stats:", np.nanmin(kappa), np.nanmean(kappa), np.nanmax(kappa))
    print("DIAG: EI_total:", EI_total)
    print("DIAG: sb z_skin_top/bot:", sb.get('z_skin_top'), sb.get('z_skin_bot'))
    print("DIAG: sb z_rod_top/bot:", sb.get('z_rod_top'), sb.get('z_rod_bot'))
    print("DIAG: sigma_rod_top stats BEFORE cap:", np.nanmin(sigma_rod_top), np.nanmean(sigma_rod_top), np.nanmax(sigma_rod_top))
    print("DIAG: sigma_skin_top stats BEFORE cap:", np.nanmin(sigma_skin_top), np.nanmean(sigma_skin_top), np.nanmax(sigma_skin_top))

    # Diagnostic: prove rod stiffness effect via curvature reduction
    M_root = float(np.atleast_1d(M)[0])
    kappa_root = float(np.atleast_1d(kappa)[0])
    EI_skins_only = sb.get('EI_skins', 0.0) + sb.get('EI_core', 0.0)
    kappa_no_rods = M_root / max(EI_skins_only, 1e-18) if EI_skins_only > 0 else np.inf

    print(f"DIAG: Curvature with rods: kappa = {kappa_root:.3e} 1/m")
    print(f"DIAG: Curvature without rods (skins+core only): kappa = {kappa_no_rods:.3e} 1/m")
    if np.isfinite(kappa_no_rods) and kappa_no_rods > 0:
        reduction_factor = kappa_no_rods / max(kappa_root, 1e-18)
        print(f"DIAG: Rod stiffness reduces curvature by {reduction_factor:.1f}x (confirms rod dominance)")

    # FoS calculations (robust)
    FoS_skin_top = fos_normal(sigma_skin_top, mat.Xt_skin, mat.Xc_skin)
    FoS_skin_bot = fos_normal(sigma_skin_bot, mat.Xt_skin, mat.Xc_skin)
    FoS_rod_top = fos_normal(sigma_rod_top, mat.Xt_rod, mat.Xc_rod)
    FoS_rod_bot = fos_normal(sigma_rod_bot, mat.Xt_rod, mat.Xc_rod)

    # ---- Sandwich face-sheet wrinkling checks ----
    # Critical wrinkling stress using Allen/Zenkert continuous foundation model
    # E_face_axial already computed above for stress calculation
    
    # Base wrinkling stress (symmetric sandwich, c_coeff per Zenkert §6.4)
    # c=0.5 is anti-symmetric (one face only), c=0.7-0.9 is symmetric (both faces)
    sigma_wr_base = sandwich_wrinkling_stress(
        E_face_axial=E_face_axial,
        E_core=float(mat.E_core),
        G_core=float(mat.G_core),
        nu_face=mat.nu12_skin,
        c_coeff=0.9  # Symmetric sandwich (both faces active, upper bound)
    )
    
    # ---- Discrete stiffener correction for wrinkling ----
    # Pultruded UD carbon rods (138 GPa, 3×3mm) are bonded into sandwich structure
    # They provide discrete nodal restraints that prevent face wrinkling
    # 
    # Physical mechanism:
    # 1. Rods are orders of magnitude stiffer than thin face (0.26mm biax)
    # 2. Face must deform around bonded rods → cannot wrinkle freely
    # 3. Effect scales with rod spacing and bond quality
    #
    # Literature range for bonded discrete stiffeners at ~40mm spacing: 1.5-2.5×
    # Conservative mid-range estimate: 2.0×
    # 
    # Validation required: Panel test or FEA with explicit rod geometry
    
    rod_stiffening_factor = 2.0  # Empirical for bonded pultruded rods at ~40mm spacing
    
    # Design wrinkling stress (includes rod benefit)
    sigma_wr_stiffened = sigma_wr_base * rod_stiffening_factor
    
    # Rod spacing for documentation
    rod_positions_xc = [0.15, 0.15, 0.55, 0.55]  # x/c positions from geometry
    rod_x = [xc * chord_m for xc in rod_positions_xc]
    rod_spacings = [abs(rod_x[i+1] - rod_x[i]) for i in range(len(rod_x)-1) if abs(rod_x[i+1] - rod_x[i]) > 1e-9]
    s_min = min(rod_spacings) if rod_spacings else chord_m  # m
    
    print(f"\n=== WRINKLING ANALYSIS (REVISED) ===")
    print(f"Base wrinkling stress (symmetric sandwich, c=0.9): {sigma_wr_base/1e6:.1f} MPa")
    print(f"Rod spacing (chordwise): {s_min*1000:.1f} mm")
    print(f"Rod stiffening factor: {rod_stiffening_factor:.1f}× (bonded pultruded UD rods)")
    print(f"Design wrinkling stress (with rods): {sigma_wr_stiffened/1e6:.1f} MPa")
    print(f"Physical basis: Bonded stiffeners constrain wrinkling wavelength")
    print(f"                and provide local reinforcement")
    print(f"Validation: Requires panel test or FEA for final certification")
    print(f"========================================")
    
    # Wrinkling FoS: sigma_wr_stiffened / |sigma_face_compressive|
    # Check both top and bottom (both can be in compression depending on load case)
    
    eps_sigma = 1e-6  # Consistent threshold for compression check and denominator floor
    
    # Top face wrinkling (only relevant when in compression)
    FoS_skin_wrinkling_top = np.where(
        sigma_skin_top < -eps_sigma,  # Only check compressive stress
        sigma_wr_stiffened / np.maximum(np.abs(sigma_skin_top), eps_sigma),
        np.inf  # Tension → no wrinkling risk
    )
    
    # Bottom face wrinkling (only relevant when in compression)
    FoS_skin_wrinkling_bot = np.where(
        sigma_skin_bot < -eps_sigma,  # Only check compressive stress
        sigma_wr_stiffened / np.maximum(np.abs(sigma_skin_bot), eps_sigma),
        np.inf  # Tension → no wrinkling risk
    )
    
    print(f"\nINFO: Wrinkling check ENABLED (face stress from curvature strain, E_face_axial=A11/t)")
    print(f"      Includes discrete stiffener effects from rods")
    print(f"DEBUG: E_face_axial = {E_face_axial/1e9:.2f} GPa")
    print(f"DEBUG: sigma_skin_top[0] = {float(np.atleast_1d(sigma_skin_top)[0])/1e6:.1f} MPa")
    print(f"      FoS_wrinkling_top: min={np.nanmin(FoS_skin_wrinkling_top):.2f}, mean={np.nanmean(FoS_skin_wrinkling_top):.2f}")
    print(f"      FoS_wrinkling_bot: min={np.nanmin(FoS_skin_wrinkling_bot):.2f}, mean={np.nanmean(FoS_skin_wrinkling_bot):.2f}")
    
    # Web shear buckling (infinite-strip model, no discrete ribs)
    tau_cr_web = web_shear_buckling_stress_infinite_strip(
        E_web=float(mat.E1_web),
        h_web=h_web_effective,  # Vertical clear height (rod-to-rod)
        t_web=t_web_actual,
        nu_web=mat.nu12_web
    )
    
    print(f"DEBUG: Web shear buckling (infinite strip): tau_cr_web = {tau_cr_web/1e6:.1f} MPa")
    print(f"       (Conservative for continuous sandwich, no discrete panels)")
    
    # Web buckling FoS (use combined transverse + torsion stress)
    FoS_web_buckling = tau_cr_web / np.maximum(tau_web_shear_total, 1e-6)

    # Safe denominator floors (small but nonzero) to avoid infinite FoS when stresses are zero,
    # while still allowing sensitivity to small changes. Use physically sensible eps based on allowables.
    eps_web = max(1e-6, 1e-6 * float(mat.tau_web_allow))
    eps_skin_tor = max(1e-6, 1e-6 * float(mat.tau_skin_allow))
    eps_core = max(1e-9, 1e-6 * float(mat.tau_core_allow))

    # Compute simple shear FoS values using those floors
    # In no-webs mode, set web FoS to infinity so it doesn't govern
    no_webs_flag = getattr(inp, '_no_webs', False)
    print(f"\n*** DEBUG: _no_webs flag = {no_webs_flag} ***")
    print(f"*** DEBUG: tau_web_shear_total stats: min={np.nanmin(tau_web_shear_total):.3g}, max={np.nanmax(tau_web_shear_total):.3g} ***")
    
    if no_webs_flag:
        FoS_web_shear = np.full_like(tau_web_shear_total, 1e10)  # effectively infinite (web doesn't exist)
        print(f"*** NO-WEBS MODE ACTIVE: FoS_web_shear set to 1e10 ***")
        print(f"*** DEBUG: FoS_web_shear after setting: min={np.nanmin(FoS_web_shear):.3g}, max={np.nanmax(FoS_web_shear):.3g} ***")
    else:
        FoS_web_shear = float(mat.tau_web_allow) / np.maximum(np.abs(tau_web_shear_total), eps_web)
        print(f"*** NORMAL MODE: FoS_web_shear calculated from tau_web_shear_total (RSS of transverse + torsion) ***")
        print(f"*** DEBUG: FoS_web_shear after calc: min={np.nanmin(FoS_web_shear):.3g}, max={np.nanmax(FoS_web_shear):.3g} ***")
    
    FoS_skin_tor_simple = float(mat.tau_skin_allow) / np.maximum(np.abs(tau_skin_tor), eps_skin_tor)
    FoS_core_shear = float(mat.tau_core_allow) / np.maximum(np.abs(tau_core), eps_core)
    
    # Skin membrane shear: NOT physical in sandwich theory (removed from model)
    # Set to non-limiting value for backward compatibility in data structures
    FoS_skin_membrane = np.full_like(V, 1e10, dtype=float)  # Non-limiting
    
    # DEBUG: Print FoS statistics
    print(f"\n=== FOS STATISTICS ===")
    print(f"FoS_skin_top: min={np.nanmin(FoS_skin_top):.3f}, mean={np.nanmean(FoS_skin_top):.3f}")
    print(f"FoS_skin_bot: min={np.nanmin(FoS_skin_bot):.3f}, mean={np.nanmean(FoS_skin_bot):.3f}")
    print(f"FoS_rod_top: min={np.nanmin(FoS_rod_top):.3f}, mean={np.nanmean(FoS_rod_top):.3f}")
    print(f"FoS_rod_bot: min={np.nanmin(FoS_rod_bot):.3f}, mean={np.nanmean(FoS_rod_bot):.3f}")
    if not getattr(inp, '_no_webs', False):
        print(f"FoS_web_shear: min={np.nanmin(FoS_web_shear):.3f}, mean={np.nanmean(FoS_web_shear):.3f}")
    print(f"FoS_skin_tor: min={np.nanmin(FoS_skin_tor_simple):.3f}, mean={np.nanmean(FoS_skin_tor_simple):.3f}")
    print(f"FoS_core_shear: min={np.nanmin(FoS_core_shear):.3f}, mean={np.nanmean(FoS_core_shear):.3f}")
    # FoS_skin_membrane removed from model (not physical in sandwich theory)

    # Tsai-Wu composite check for skin (vectorized)
    # FIXED: Use conservative txy=0 (tau_skin_tor is Bredt wall stress, not laminate in-plane shear)
    tsai_idx_top = np.array([
        tsai_wu_index(
            sx=float(sigma_skin_top[i]) if np.isfinite(sigma_skin_top[i]) else 0.0,
            sy=0.0,
            txy=0.0,  # Conservative: ignore shear in Tsai-Wu (Bredt stress is not laminate stress)
            Xt=mat.Xt_skin,
            Xc=mat.Xc_skin,
            S=mat.tau_skin_allow
        )
        for i in range(len(M))
    ], dtype=float)

    # Invert only where tsai_idx_top > 0; treat non-positive index as non-limiting
    pos = tsai_idx_top > 1e-12
    FoS_skin_tor_tsai = np.full_like(tsai_idx_top, np.inf, dtype=float)
    FoS_skin_tor_tsai[pos] = 1.0 / tsai_idx_top[pos]

    # Merge shear-only FoS and Tsai–Wu: take the more conservative (minimum) while treating
    # "no Tsai limit" as non-limiting (inf)
    FoS_skin_tor = np.minimum(FoS_skin_tor_simple, FoS_skin_tor_tsai)

    # Clip to reasonable upper bound to avoid numerical infinities dominating plots/tables
    FoS_skin_tor = np.clip(FoS_skin_tor, 0.0, 1e4)

    # Diagnostics: print summary stats so you can see whether distributions are varying with sweep
    try:
        print("\n=== FoS DIAGNOSTICS ===")
        print(f"FoS_skin_tor: min={np.nanmin(FoS_skin_tor):.3g}, mean={np.nanmean(FoS_skin_tor):.3g}, max={np.nanmax(FoS_skin_tor):.3g}")
        print(f"FoS_core_shear: min={np.nanmin(FoS_core_shear):.3g}, mean={np.nanmean(FoS_core_shear):.3g}, max={np.nanmax(FoS_core_shear):.3g}")
        print(f"tau_skin_tor: min={np.nanmin(tau_skin_tor):.3g}, mean={np.nanmean(tau_skin_tor):.3g}, max={np.nanmax(tau_skin_tor):.3g}")
        print(f"tau_core: min={np.nanmin(tau_core):.3g}, mean={np.nanmean(tau_core):.3g}, max={np.nanmax(tau_core):.3g}")
        print("========================\n")
    except Exception:
        pass

    # Merge with simple shear-only FoS (take the more critical)
    FoS_skin_tor = np.minimum(FoS_skin_tor, FoS_skin_tor_tsai)
    
    # (optional) quick sanity
    # print("Tsai-Wu idx min/max:", np.nanmin(tsai_idx_top), np.nanmax(tsai_idx_top))
    # print("FoS_skin_tor_tsai  :", np.nanmin(FoS_skin_tor_tsai), np.nanmax(FoS_skin_tor_tsai))

    # Cap FoS values to prevent unrealistic infinities
    def _cap(a):
        a = np.asarray(a, dtype=float)
        a[~np.isfinite(a)] = np.nan
        a = np.clip(a, 0.0, 1e5)
        return a

    FoS_skin_top = _cap(FoS_skin_top)
    FoS_skin_bot = _cap(FoS_skin_bot)
    FoS_rod_top = _cap(FoS_rod_top)
    FoS_rod_bot = _cap(FoS_rod_bot)
    # Don't cap FoS_web_shear in no-webs mode (keep it at 1e10 so it never governs)
    no_webs_at_cap = getattr(inp, '_no_webs', False)
    print(f"\n*** DEBUG AT CAPPING: _no_webs = {no_webs_at_cap} ***")
    print(f"*** FoS_web_shear BEFORE capping: min={np.nanmin(FoS_web_shear):.3g}, max={np.nanmax(FoS_web_shear):.3g} ***")
    
    if not no_webs_at_cap:
        FoS_web_shear = _cap(FoS_web_shear)
        print(f"*** FoS_web_shear AFTER capping: min={np.nanmin(FoS_web_shear):.3g}, max={np.nanmax(FoS_web_shear):.3g} ***")
    else:
        print(f"*** NO-WEBS MODE: Skipping cap for FoS_web_shear (keeping at 1e10) ***")
    
    FoS_skin_tor = _cap(FoS_skin_tor)
    # (patch) do not cap FoS_core_shear so trends show clearly
    def minpos(a):
        s = np.asarray(a)[np.isfinite(a) & (np.asarray(a) > 0)]
        return float(s.min()) if s.size else float('nan')

    # Build overall FoS dictionary (exclude web shear in no-webs mode)
    # Note: Skin membrane shear removed - not physical in sandwich theory
    overall = {
        "Skin top": minpos(FoS_skin_top),
        "Skin bottom": minpos(FoS_skin_bot),
        "Rod top": minpos(FoS_rod_top),
        "Rod bottom": minpos(FoS_rod_bot),
        "Skin torsion": minpos(FoS_skin_tor),
        "Core shear": minpos(FoS_core_shear),
        "Skin wrinkling top": minpos(FoS_skin_wrinkling_top),
        "Skin wrinkling bot": minpos(FoS_skin_wrinkling_bot),
        "Web shear buckling": minpos(FoS_web_buckling),
    }
    
    # Only include web shear if not in no-webs mode
    no_webs_at_overall = getattr(inp, '_no_webs', False)
    print(f"\n*** DEBUG AT OVERALL DICT: _no_webs = {no_webs_at_overall} ***")
    print(f"*** FoS_web_shear value: min={np.nanmin(FoS_web_shear):.3g}, minpos={minpos(FoS_web_shear):.3g} ***")
    
    if not no_webs_at_overall:
        web_shear_fos = minpos(FoS_web_shear)
        overall["Web shear"] = web_shear_fos
        print(f"*** ADDING Web shear to overall dict with FoS = {web_shear_fos:.3g} ***")
    else:
        print(f"*** NO-WEBS MODE: EXCLUDING Web shear from overall dict ***")
    
    min_fos = min([v for v in overall.values() if np.isfinite(v) and v > 0] or [np.nan])
    max_g = min_fos * inp.g_limit if np.isfinite(min_fos) else np.nan
    
    print(f"\n=== GOVERNING FAILURE ANALYSIS ===")
    print(f"No-webs mode: {getattr(inp, '_no_webs', False)}")
    if getattr(inp, '_no_webs', False):
        print("  (Web shear excluded from analysis - webs neutralized)")
    for mode, fos_val in sorted(overall.items(), key=lambda x: x[1]):
        print(f"  {mode:20s}: FoS = {fos_val:.3f}")
    print(f"Minimum FoS: {min_fos:.3f} (mode: {min(overall, key=overall.get) if overall else 'N/A'})")

    print("\n===== FACTOR OF SAFETY DIAGNOSTICS =====")
    x_stations = xspan.copy()
    report_low_fos("FoS_skin_top", FoS_skin_top, x_stations, threshold=3.0)
    report_low_fos("FoS_skin_bot", FoS_skin_bot, x_stations, threshold=3.0)
    report_low_fos("FoS_rod_top", FoS_rod_top, x_stations, threshold=3.0)
    report_low_fos("FoS_rod_bot", FoS_rod_bot, x_stations, threshold=3.0)
    if not getattr(inp, '_no_webs', False):
        report_low_fos("FoS_web_shear", FoS_web_shear, x_stations, threshold=3.0)
    report_low_fos("FoS_skin_tor", FoS_skin_tor, x_stations, threshold=3.0)
    report_low_fos("FoS_core_shear", FoS_core_shear, x_stations, threshold=3.0)

    print("\n===== DRIVER STRESS STATISTICS =====")
    drivers = {
        "tau_web1_arr": tau_web1_arr,
        "tau_web2_arr": tau_web2_arr,
        "tau_web_shear": tau_web_shear,
        "tau_skin_tor": tau_skin_tor,
        "tau_core": tau_core,
        "tau_skin_membrane": tau_skin_membrane,
        "I_trans_arr": I_trans_arr,
        "sigma_skin_top": sigma_skin_top,
        "sigma_rod_top": sigma_rod_top
    }
    for dname, arr in drivers.items():
        a = np.asarray(arr, dtype=float)
        valid = np.isfinite(a)
        if not valid.any():
            print(f"  {dname}: no finite values")
        else:
            print(f"  {dname}: {float(np.nanmin(a)):.3g} / {float(np.nanmean(a)):.3g} / {float(np.nanmax(a)):.3g}")

    print("\n===== CELL AND WALL INFO =====")
    print(f"A_cells (original): {A_cells}")
    print(f"A_cells (torsion): {A_cells_torsion}")
    print("Walls (torsion topology):")
    for j, w in enumerate(walls_torsion):
        print(f"  wall {j}: length={w.get('length'):.6g} m, thickness={w.get('thickness'):.6g} m, cell_ids={w.get('cell_ids')}")

    # Diagnostic plot (optional)
    if inp.make_plots:
        try:
            fig, ax = plt.subplots(2, 1, figsize=(10, 8))
            ax[0].plot(x_stations, tau_web1_arr, label='tau_web1', linewidth=2)
            ax[0].plot(x_stations, tau_web2_arr, label='tau_web2', linewidth=2)
            ax[0].plot(x_stations, tau_web_shear, label='tau_web_shear (avg)', linestyle='--', linewidth=2)
            ax[0].axhline(mat.tau_web_allow, color='r', linestyle=':', label='tau_web_allow')
            ax[0].legend(); ax[0].set_ylabel('Shear Stress (Pa)'); ax[0].set_title('Web Shear Stress Drivers'); ax[0].grid(True)
            ax[1].plot(x_stations, tau_skin_tor, label='tau_skin_tor', linewidth=2)
            ax[1].plot(x_stations, tau_core, label='tau_core', alpha=0.7, linewidth=2)
            ax[1].axhline(mat.tau_skin_allow, color='r', linestyle=':', label='tau_skin_allow')
            ax[1].axhline(mat.tau_core_allow, color='orange', linestyle=':', label='tau_core_allow')
            ax[1].legend(); ax[1].set_xlabel('Span position (m)'); ax[1].set_ylabel('Shear Stress (Pa)'); ax[1].set_title('Torsion and Core Shear Stress Drivers'); ax[1].grid(True)
            plt.tight_layout()
            diag_path = os.path.join(".", "diagnostic_drivers.png")
            fig.savefig(diag_path, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"\nDiagnostic plot saved to: {diag_path}")
        except Exception as e:
            print(f"Failed to save diagnostic plot: {e}")
        # Diagnostics: report nontrivial counts and percentiles for torsion/core/lift distributions
    try:
        tau_skin_tor_arr = np.asarray(tau_skin_tor) if ('tau_skin_tor' in locals() or 'tau_skin_tor' in globals()) else None
        tau_core_arr     = np.asarray(tau_core)     if ('tau_core' in locals() or 'tau_core' in globals()) else None
        lift_arr         = np.asarray(lift_per_unit_span) if ('lift_per_unit_span' in locals() or 'lift_per_unit_span' in globals()) else None

        if tau_skin_tor_arr is None or tau_core_arr is None or lift_arr is None:
            print("DIAG: one or more diagnostic arrays not yet defined:",
                  "tau_skin_tor defined?", tau_skin_tor_arr is not None,
                  "tau_core defined?", tau_core_arr is not None,
                  "lift_per_unit_span defined?", lift_arr is not None)
        else:
            n_nontriv_skin = int(np.count_nonzero(np.abs(tau_skin_tor_arr) > 1e-6))
            n_nontriv_core = int(np.count_nonzero(np.abs(tau_core_arr) > 1e-6))
            print(f"DIAG: tau_skin_tor nontrivial: {n_nontriv_skin} of {len(tau_skin_tor_arr)}")
            print(f"DIAG: tau_core nontrivial: {n_nontriv_core} of {len(tau_core_arr)}")
            print("DIAG percentiles (tau_skin_tor) p10,p50,p90:", np.nanpercentile(tau_skin_tor_arr, [10,50,90]))
            print("DIAG percentiles (tau_core)     p10,p50,p90:", np.nanpercentile(tau_core_arr, [10,50,90]))
            print(f"DIAG: lift_per_unit_span[0] = {float(lift_arr[0])} N/m")
    except Exception as e:
        print("DIAG: diagnostics failed:", e)
        
    print("="*50 + "\n")

    # Build DataFrames
    def _to_array(x, n):
        if isinstance(x, np.ndarray):
            return x.astype(float)
        elif isinstance(x, (list, tuple)):
            return np.array(x, dtype=float)
        else:
            try:
                return np.full(n, float(x))
            except:
                return np.zeros(n, dtype=float)

    xspan = _to_array(xspan, n)
    dx = _to_array(dx, n)
    lift_per_unit_span = _to_array(lift_per_unit_span, n)
    V = _to_array(V, n)
    M = _to_array(M, n)
    T = _to_array(T, n)
    sigma_skin_top = _to_array(sigma_skin_top, n)
    sigma_skin_bot = _to_array(sigma_skin_bot, n)
    sigma_rod_top = _to_array(sigma_rod_top, n)
    sigma_rod_bot = _to_array(sigma_rod_bot, n)
    tau_web_shear = _to_array(tau_web_shear, n)
    tau_skin_tor = _to_array(tau_skin_tor, n)
    tau_web1 = _to_array(tau_web1_arr, n)
    tau_web2 = _to_array(tau_web2_arr, n)
    tau_core = _to_array(tau_core, n)
    tau_skin_membrane = _to_array(tau_skin_membrane, n)
    FoS_skin_top = _to_array(FoS_skin_top, n)
    FoS_skin_bot = _to_array(FoS_skin_bot, n)
    FoS_rod_top = _to_array(FoS_rod_top, n)
    FoS_rod_bot = _to_array(FoS_rod_bot, n)
    FoS_web_shear = _to_array(FoS_web_shear, n)
    FoS_skin_tor = _to_array(FoS_skin_tor, n)
    FoS_core_shear = _to_array(FoS_core_shear, n)
    FoS_skin_membrane = _to_array(FoS_skin_membrane, n)
    FoS_skin_wrinkling_top = _to_array(FoS_skin_wrinkling_top, n)
    FoS_skin_wrinkling_bot = _to_array(FoS_skin_wrinkling_bot, n)
    FoS_web_buckling = _to_array(FoS_web_buckling, n)

    stations = pd.DataFrame({"i": np.arange(n), "x (m)": xspan, "dx (m)": dx})
    loads = pd.DataFrame({"i": np.arange(n), "x (m)": xspan, "lift_per_unit_span (N/m)": lift_per_unit_span, "V (N)": V, "M (N·m)": M, "T (N·m)": T})
    stresses = pd.DataFrame({
        "i": np.arange(n), "x (m)": xspan,
        "sigma_skin_top (Pa)": sigma_skin_top, "sigma_skin_bot (Pa)": sigma_skin_bot,
        "sigma_rod_top (Pa)": sigma_rod_top, "sigma_rod_bot (Pa)": sigma_rod_bot,
        "tau_web_shear (Pa)": tau_web_shear, "tau_skin_tor (Pa)": tau_skin_tor,
        "tau_web1_tor (Pa)": tau_web1, "tau_web2_tor (Pa)": tau_web2, "tau_core (Pa)": tau_core,
        "tau_skin_membrane (Pa)": tau_skin_membrane
    })
    fos_df = pd.DataFrame({
        "i": np.arange(n), "x (m)": xspan,
        "FoS_skin_top": FoS_skin_top, "FoS_skin_bot": FoS_skin_bot,
        "FoS_rod_top": FoS_rod_top, "FoS_rod_bot": FoS_rod_bot,
        "FoS_web_shear": FoS_web_shear, "FoS_skin_tor": FoS_skin_tor, "FoS_core_shear": FoS_core_shear,
        "FoS_skin_membrane": FoS_skin_membrane,
        "FoS_skin_wrinkling_top": FoS_skin_wrinkling_top,
        "FoS_skin_wrinkling_bot": FoS_skin_wrinkling_bot,
        "FoS_web_buckling": FoS_web_buckling,
    })
    summary = pd.DataFrame({
        "Parameter": [
            "Chord (m)", "Span (m)", "Sweep (deg)",
            "Web1 position (x/c)", "Web2 position (x/c)",
            "Web1 height (m)", "Web2 height (m)",
            "Skin thickness per face (mm)", "Core thickness (mm)",
            "EI_total (N·m^2)", "Min FoS @g_limit", "Max G capability", "Governing failure mode"
        ],
        "Value": [
            chord_m, span_m, inp.sweep_deg,
            inp.web1_p, inp.web2_p,
            h1, h2,
            sb.get('t_skin_total', t_face) * 1000, t_core * 1000,
            EI_total, min_fos, max_g,
            (min(overall, key=overall.get) if overall else "N/A")
        ]
    })

    return {
        "Inputs": pd.DataFrame([asdict(inp)]),
        "Materials": pd.DataFrame([asdict(mat)]),
        "Stations": stations,
        "Loads": loads,
        "Stresses": stresses,
        "FoS_per_station": fos_df,
        "Summary": summary,
        "Airfoil_xy": pd.DataFrame(raw_outline, columns=["x (m)", "y (m)"])
    }

# ---------------------- PLOTTING / IO ----------------------

def create_equations_appendix():
    eq = {
        "Category": [
            "Sandwich Bending", "Skin parallel-axis", "Rod contribution",
            "Core bending", "Bending stress", "Shear partition",
            "Torsion equilibrium", "Bredt shear flow", "Shear stress"
        ],
        "Equation": [
            "EI_total = Σ E_i * (I_i + A_i * z_i^2)",
            "EI_skin = laminate D contribution * z^2",
            "EI_rod = n_rods * E_rod * (I_self + A_rod * z^2)",
            "I_core = b * t^3 / 12",
            "σ = M * z * E_material / EI_total (approx)",
            "k_i = (G_i * t_i) / Σ(G_j * t_j) ; V_i = k_i * V_total",
            "T = 2 * Σ (q_j * A_cell_contrib_j)",
            "Solve linear system for wall q_j using twist compatibility and torque equilibrium",
            "τ = q / t ; τ_web = V*Q/(I_trans * t_web) for vertical webs"
        ],
        "Description": [
            "Parallel-axis sum of skins, rods and core bending terms",
            "Face sheets treated with CLT D matrix",
            "Discrete spar caps as rectangular rods",
            "Core self-inertia (small)",
            "Beam bending law with composite laminate scaling",
            "Stiffness-based shear partition: skins, webs, core act as parallel springs",
            "Torque balance across cells",
            "Bredt-Batho linear solve (multi-cell)",
            "Convert shear flow to shear stress by dividing by thickness"
        ]
    }
    return pd.DataFrame(eq)

def plot_all_results(results: Dict, inp: Inputs, output_dir: str = "."):
    all_sweeps = results.get("All_Sweeps", pd.DataFrame())
    af_df = results.get("Airfoil_xy", pd.DataFrame())

    os.makedirs(output_dir, exist_ok=True)

    # ---------------- Airfoil cross-section with webs ----------------
    fig, ax = plt.subplots(figsize=(5, 4))
    if af_df is None or af_df.empty:
        ax.text(0.5, 0.5, "No airfoil data", ha="center", va="center")
    else:
        af = af_df.copy()
        xs = af["x (m)"].to_numpy(dtype=float)
        ys = af["y (m)"].to_numpy(dtype=float)

        # Close if not already closed
        if not np.allclose([xs[0], ys[0]], [xs[-1], ys[-1]], atol=1e-12):
            xs = np.append(xs, xs[0])
            ys = np.append(ys, ys[0])

        ax.plot(xs, ys, "k-", linewidth=2)

        chord = inp.chord_m if inp.chord_m is not None else (xs.max() - xs.min())
        x_web1 = xs.min() + inp.web1_p * chord
        x_web2 = xs.min() + inp.web2_p * chord

        ax.axvline(x_web1, linestyle="--", color="r", label=f"Web1 @ {inp.web1_p:.0%}c")
        ax.axvline(x_web2, linestyle="--", color="b", label=f"Web2 @ {inp.web2_p:.0%}c")
        ax.set_xlabel("Chordwise Position (m)")
        ax.set_ylabel("Height (m)")
        ax.set_title("Airfoil Cross-Section with Web Locations")
        ax.axis("equal")
        ax.grid(True)
        ax.legend(fontsize=8)

    airfoil_path = os.path.join(output_dir, "airfoil_with_webs.png")
    fig.savefig(airfoil_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {airfoil_path}")

    # If no sweep data, we are done
    if all_sweeps is None or all_sweeps.empty:
        print("No All_Sweeps data – skipping spanwise plots.")
        return

    # ---------------- Governing FoS column (min of available FoS) ----------------
    fos_cols = [
        c for c in [
            "FoS_skin_top", "FoS_skin_bot",
            "FoS_rod_top", "FoS_rod_bot",
            "FoS_web_shear", "FoS_skin_tor",
            "FoS_core_shear", "FoS_skin_membrane",
            "FoS_skin_buckling",
        ]
        if c in all_sweeps.columns
    ]
    if fos_cols:
        all_sweeps["FoS_governing"] = all_sweeps[fos_cols].min(axis=1)

    # ---------------- Helper to save one curve family per file ----------------
    def save_group(col: str,
                   scale: float,
                   ylabel: str,
                   filename: str,
                   ylog: bool = False,
                   ylim: Optional[Tuple[float, float]] = None):
        if col not in all_sweeps.columns:
            print(f"Skipping {col}: column not found in All_Sweeps.")
            return

        fig, ax = plt.subplots(figsize=(5, 4))
        for s, grp in all_sweeps.groupby("sweep_deg"):
            ax.plot(
                grp["x (m)"],
                grp[col] * scale,
                label=f"Λ={int(float(s))}°",
                linewidth=1.2,
            )
        ax.set_xlabel("Spanwise Position from Root (m)")
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.grid(True)
        ax.legend(fontsize=7)

        if ylog:
            ax.set_yscale("log")
        if ylim is not None:
            ax.set_ylim(*ylim)

        outpath = os.path.join(output_dir, filename)
        fig.savefig(outpath, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {outpath}")

    # ---------------- One file per graph ----------------
    # Loads
    save_group("M (N·m)",            1.0,  "Bending Moment (N·m)",          "bending_moment_vs_span.png")
    save_group("V (N)",              1.0,  "Shear Force (N)",               "shear_force_vs_span.png")
    save_group("T (N·m)",            1.0,  "Torque (N·m)",                  "torque_vs_span.png")

    # Normal stresses (Pa -> MPa)
    save_group("sigma_rod_top (Pa)", 1e-6, "Rod Stress Top (MPa)",          "rod_stress_top_vs_span.png")
    save_group("sigma_rod_bot (Pa)", 1e-6, "Rod Stress Bottom (MPa)",       "rod_stress_bottom_vs_span.png")
    save_group("sigma_skin_top (Pa)",1e-6, "Skin Stress Top (MPa)",         "skin_stress_top_vs_span.png")
    save_group("sigma_skin_bot (Pa)",1e-6, "Skin Stress Bottom (MPa)",      "skin_stress_bottom_vs_span.png")

    # Shear / torsion stresses (Pa -> MPa)
    # Skip web plots if --no_webs mode
    if not getattr(inp, '_no_webs', False):
        save_group("tau_web_shear (Pa)", 1e-6, "Web Shear Stress (MPa)",        "web_shear_stress_vs_span.png")
        save_group("tau_web1_tor (Pa)",  1e-6, "Web1 Torsion Stress (MPa)",     "web1_torsion_stress_vs_span.png")
        save_group("tau_web2_tor (Pa)",  1e-6, "Web2 Torsion Stress (MPa)",     "web2_torsion_stress_vs_span.png")
    save_group("tau_skin_tor",       1e-6, "Skin Torsion Stress (MPa)",     "skin_torsion_stress_vs_span.png")
    save_group("tau_core (Pa)",      1e-6, "Core Shear Stress (MPa)",       "core_shear_stress_vs_span.png")
    save_group("tau_skin_membrane (Pa)", 1e-6, "Skin Membrane Shear Stress (MPa)", "skin_membrane_shear_vs_span.png")

    # Factors of Safety (log scale)
    save_group("FoS_rod_top",        1.0,  "FoS: Rod Top",                  "FoS_rod_top_vs_span.png",        ylog=True)
    save_group("FoS_rod_bot",        1.0,  "FoS: Rod Bottom",               "FoS_rod_bottom_vs_span.png",     ylog=True)
    save_group("FoS_skin_top",       1.0,  "FoS: Skin Top",                 "FoS_skin_top_vs_span.png",       ylog=True)
    save_group("FoS_skin_bot",       1.0,  "FoS: Skin Bottom",              "FoS_skin_bottom_vs_span.png",    ylog=True)
    if not getattr(inp, '_no_webs', False):
        save_group("FoS_web_shear",      1.0,  "FoS: Web Shear",                "FoS_web_shear_vs_span.png",      ylog=True)
    save_group("FoS_skin_tor",       1.0,  "FoS: Skin Torsion",             "FoS_skin_torsion_vs_span.png",   ylog=True)
    save_group("FoS_core_shear",     1.0,  "FoS: Core Shear",               "FoS_core_shear_vs_span.png",     ylog=True)
    save_group("FoS_skin_membrane",  1.0,  "FoS: Skin Membrane Shear",      "FoS_skin_membrane_vs_span.png",  ylog=True)
    save_group("FoS_skin_buckling",  1.0,  "FoS: Skin Buckling (Top)",      "FoS_skin_buckling_vs_span.png",  ylog=True)

    # Governing FoS (min over all FoS types)
    if "FoS_governing" in all_sweeps.columns:
        save_group("FoS_governing",  1.0,  "Governing Factor of Safety",    "FoS_governing_vs_span.png",      ylog=True)


def run_sweeps_professional(inp_base: Inputs, mat: Materials, sweeps: List[float]):
    combined_rows = []
    summaries_rows = []
    per_sweep = {}
    last_res = None

    for s in sweeps:
        print(f"\n{'='*60}")
        print(f"Running sweep angle: {s}°")
        print('='*60)
        inp = Inputs(**{**asdict(inp_base), "sweep_deg": float(s)})
        # Preserve the _no_webs flag (not part of dataclass fields)
        if hasattr(inp_base, '_no_webs'):
            inp._no_webs = inp_base._no_webs
        res = run_analysis_professional(inp, mat)
        print(f"Sweep {s}°: T root = {res['Loads']['T (N·m)'].iloc[0]:.3g}, M root = {res['Loads']['M (N·m)'].iloc[0]:.3g}")
        last_res = res
        merged = res["Loads"].merge(res["Stresses"], on=["i", "x (m)"]).merge(res["FoS_per_station"], on=["i", "x (m)"])
        # create a column that is explicitly the sweep value repeated for each station
        merged["sweep_deg"] = float(s)
        # ensure dtype numeric
        merged["sweep_deg"] = pd.to_numeric(merged["sweep_deg"], errors="coerce")
        # put sweep_deg first for readability
        cols = merged.columns.tolist()
        cols = ["sweep_deg"] + [c for c in cols if c != "sweep_deg"]
        merged = merged[cols]
        combined_rows.append(merged)
        val_map = dict(zip(res["Summary"]["Parameter"], res["Summary"]["Value"]))
        val_map["Sweep (deg)"] = float(s)
        summaries_rows.append(val_map)
        tag = f"{int(s):02d}"
        per_sweep[f"Loads_Λ{tag}"] = res["Loads"]
        per_sweep[f"Stresses_Λ{tag}"] = res["Stresses"]
        per_sweep[f"FoS_Λ{tag}"] = res["FoS_per_station"]
        per_sweep[f"Summary_Λ{tag}"] = res["Summary"]

    all_df = pd.concat(combined_rows, ignore_index=True) if combined_rows else pd.DataFrame()
    summ_df = pd.DataFrame(summaries_rows) if summaries_rows else pd.DataFrame()
    equations_df = create_equations_appendix()

    if not summ_df.empty and "Min FoS @g_limit" in summ_df.columns:
        worst_idx = summ_df["Min FoS @g_limit"].idxmin()
        worst_sweep = summ_df.loc[worst_idx, "Sweep (deg)"]
        worst_fos = summ_df.loc[worst_idx, "Min FoS @g_limit"]
    else:
        worst_sweep = np.nan
        worst_fos = np.nan

    final_summary_df = pd.DataFrame([{"Worst Sweep (deg)": worst_sweep, "Min FoS @Design G": worst_fos}])

    out = {
        "Final_Summary": final_summary_df,
        "Inputs": pd.DataFrame([asdict(inp_base)]),
        "Materials": pd.DataFrame([asdict(mat)]),
        "All_Sweeps": all_df,
        "Summaries_All_Sweeps": summ_df,
        "Equations_Appendix": equations_df,
        "Airfoil_xy": last_res["Airfoil_xy"] if last_res else pd.DataFrame()
    }
    out.update(per_sweep)
    return out

# ---------------------- SIMPLE TESTS / VALIDATION ----------------------

def validate_bredt_single_cell():
    H = 0.1; B = 0.2; t = 1e-3
    A_cell = H * B
    walls = [
        {'length': B, 'thickness': t, 'cell_ids': [0]},
        {'length': H, 'thickness': t, 'cell_ids': [0]},
        {'length': B, 'thickness': t, 'cell_ids': [0]},
        {'length': H, 'thickness': t, 'cell_ids': [0]},
    ]
    A_cells = [A_cell]
    G_vals = [1e9] * 4
    T = 10.0
    mc = multicell_bredt_solver(A_cells, walls, G_vals, T)
    q = mc['q']
    tau = mc['tau_wall']
    tau_analytic = T / (2.0 * A_cell * t)
    return {'mc': mc, 'tau_analytic': tau_analytic, 'tau': tau}

# ---------------------- MAIN / CLI ----------------------

def main():
    parser = argparse.ArgumentParser(description="Professional wing structural analysis - CLT-enabled")
    parser.add_argument("--airfoil", type=str, default="Modified_NACA_6_Series.txt")
    parser.add_argument("--span_ft", type=float, default=4.53)
    parser.add_argument("--sweeps", type=str, default=None, help="comma separated sweep angles")
    parser.add_argument("--chord_m", type=float, default=None)
    parser.add_argument("--out", type=str, default="wing_results_final.xlsx")
    parser.add_argument("--no_plots", action="store_true")
    parser.add_argument("--no_webs", action="store_true", 
                        help="Run analysis without structural webs (webs neutralized)")
    args = parser.parse_args()
    
    # ========================================================================
    # CONFIGURATION SECTION - GEOMETRY AND FLIGHT CONDITIONS
    # ========================================================================
    # Material properties are defined in the Materials dataclass above
    # (FCIM255 fabric with T700S fiber, ROHACELL 31 IG-F core)
    
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
    
    # ========================================================================
    
    # Create Materials instance (uses dataclass defaults for FCIM255/T700S)
    if args.no_webs:
        print("\n*** NO-WEBS MODE ACTIVATED ***")
        print("Webs neutralized: ultra-thin (1 micron) and ultra-weak (1 MPa)")
        # Override web properties to neutralize them structurally
        mat = Materials(
            ply_t_web_mm=0.001,  # 1 micron (negligible)
            E1_web=1e6, E2_web=1e6, G12_web=1e5, nu12_web=0.30,
            tau_web_allow=1e6    # 1 MPa (very weak)
        )
    else:
        # Use dataclass defaults (FCIM255 with T700S fiber)
        mat = Materials()
    
    # Create Inputs instance using configuration
    inp_base = Inputs(
        airfoil_file=args.airfoil,
        span_ft=args.span_ft if args.span_ft != 4.53 else SPAN_FT,
        chord_m=CHORD_M if args.chord_m is None else float(args.chord_m),
        web1_p=WEB1_POSITION,
        web2_p=WEB2_POSITION,
        rod_width_m=ROD_WIDTH_MM / 1000.0,
        rod_height_m=ROD_HEIGHT_MM / 1000.0,
        rod_gap_mm=ROD_GAP_MM,
        g_limit=G_LIMIT,
        n_stations=N_STATIONS,
        V_ms=VELOCITY_MS,
        rho=AIR_DENSITY,
        cl=LIFT_COEFFICIENT,
        out_excel=args.out,
        make_plots=(not args.no_plots)
    )
    
    # Store no_webs flag for plot filtering
    inp_base._no_webs = args.no_webs
    
    # Diagnostic output for no-webs mode
    if args.no_webs:
        print(f"\n{'='*70}")
        print("*** NO-WEBS MODE SUMMARY ***")
        print(f"Web ply thickness: {mat.ply_t_web_mm:.3f} mm (negligible)")
        print(f"Web modulus: {mat.E1_web/1e6:.1f} MPa (ultra-weak)")
        print("Web results will be hidden from outputs")
        print(f"{'='*70}\n")
    
    # Build sweep list from configuration or CLI override
    if args.sweeps:
        sweeps = [float(s.strip()) for s in args.sweeps.split(",") if s.strip()]
    else:
        sweeps = SWEEP_ANGLES
    
    test_bredt = validate_bredt_single_cell()
    print("Bredt single-cell test status:", test_bredt['mc']['status'], "cond:", test_bredt['mc'].get('cond'))
    print("  tau analytic:", test_bredt['tau_analytic'])
    print("  tau solver:", test_bredt['tau'])

    print(f"\n{'='*70}")
    print(f"WING STRUCTURAL ANALYSIS - CLT ENABLED")
    print(f"{'='*70}")
    print(f"Running analysis for sweep angles: {sweeps}")
    print(f"{'='*70}\n")

    results = run_sweeps_professional(inp_base, mat, sweeps)

    print(f"\n{'='*70}")
    print(f"Writing results to Excel: {inp_base.out_excel}")
    print(f"{'='*70}")
    
    # Filter web-related columns if --no_webs mode is active
    if args.no_webs:
        print("Filtering web-related columns from output...")
        # Remove web-related columns from All_Sweeps
        if "All_Sweeps" in results and not results["All_Sweeps"].empty:
            web_cols = [c for c in results["All_Sweeps"].columns 
                        if 'web' in c.lower()]
            results["All_Sweeps"] = results["All_Sweeps"].drop(columns=web_cols, errors='ignore')
        
        # Filter web columns from other sheets
        for sheet_name in list(results.keys()):
            if any(x in sheet_name for x in ['FoS', 'Stresses', 'Loads']):
                df = results[sheet_name]
                if not df.empty:
                    web_cols = [c for c in df.columns if 'web' in c.lower()]
                    results[sheet_name] = df.drop(columns=web_cols, errors='ignore')

    with pd.ExcelWriter(inp_base.out_excel, engine="xlsxwriter") as writer:
        results["Final_Summary"].to_excel(writer, sheet_name="0_Final_Summary", index=False)
        results["Equations_Appendix"].to_excel(writer, sheet_name="0_Equations", index=False)
        results["Inputs"].to_excel(writer, sheet_name="1_Inputs", index=False)
        results["Materials"].to_excel(writer, sheet_name="1_Materials", index=False)
        results["Airfoil_xy"].to_excel(writer, sheet_name="1_Airfoil", index=False)
        results["All_Sweeps"].to_excel(writer, sheet_name="2_All_Sweeps", index=False)
        results["Summaries_All_Sweeps"].to_excel(writer, sheet_name="2_Summaries_All", index=False)
        # ================== BEGIN: Weight_Summary (build from existing results) ==================
        import pandas as _pd
        import numpy as _np
        
        _inp  = results["Inputs"]
        _mat  = results["Materials"]
        _summ = results["Summaries_All_Sweeps"]  # has Chord (m), Span (m), Web heights, etc.
        
        # Use the FIRST ROW (geometry is constant across sweeps)
        _chord_m     = float(_summ.loc[0, "Chord (m)"])
        _span_full_m = float(_summ.loc[0, "Span (m)"])      # this is already FULL span
        _h1_web      = float(_summ.loc[0, "Web1 height (m)"])
        _h2_web      = float(_summ.loc[0, "Web2 height (m)"])
        _t_core      = float(_mat.loc[0, "t_core"])
        _rho_core    = float(_mat.loc[0, "rho_core"])
        
        # Rod geometry from Inputs; rod density from Materials (fallback = 1600 kg/m^3)
        _rod_w   = float(_inp.loc[0, "rod_width_m"])
        _rod_h   = float(_inp.loc[0, "rod_height_m"])
        _rho_rod = float(_mat["rho_rod"].iloc[0]) if "rho_rod" in _mat.columns else 1600.0
        _n_rods_top = 2
        _n_rods_bot = 2
        _n_rods_total = _n_rods_top + _n_rods_bot
        
        # Skins: FCIM255 fabric (200 gsm total for ±45° biax), 1 fabric layer per face
        _skin_gsm_per_ply = 200.0  # FCIM255 total fabric weight (±45° biax)
        _skin_plies_per_face = 1   # number of fabric layers per face
        _skin_areal_kg_m2_per_face = (_skin_gsm_per_ply * _skin_plies_per_face) / 1e3
        _skin_total_areal_kg_m2    = 2.0 * _skin_areal_kg_m2_per_face  # top + bottom faces
        
        # Planform area (FULL wing)
        _A_planform_full_m2 = _span_full_m * _chord_m
        
        # Skins mass
        _mass_skins_kg = _skin_total_areal_kg_m2 * _A_planform_full_m2
        
        # Core mass
        _vol_core_m3  = _A_planform_full_m2 * _t_core
        _mass_core_kg = _vol_core_m3 * _rho_core
        
        # Webs: FCIM255 (HiMax CGL4012) ±45° biax labeled 200 gsm total
        _WEB_DRY_GSM        = 200.0
        _AW_IS_TOTAL_FOR_BIAX = True   # True: 200 gsm is total for ±45° biax pair
        _rho_f, _rho_r, _Vf = 1780.0, 1150.0, 0.55
        _res_over_f = ((1.0 - _Vf) * _rho_r) / (_Vf * _rho_f)              # ~0.53 at Vf=0.55
        _web_dry_each = _WEB_DRY_GSM if _AW_IS_TOTAL_FOR_BIAX else 2.0 * _WEB_DRY_GSM
        _web_cured_gsm_per_layer = _web_dry_each * (1.0 + _res_over_f)     # g/m^2 per ±45 layer
        _n_web_layers = 1  # your stack is [(+45,1),(-45,1)] -> treat as one ±45 "layer"
        _A_webs_full_m2 = _span_full_m * (_h1_web + _h2_web)               # two webs combined
        _mass_webs_kg = _n_web_layers * (_web_cured_gsm_per_layer / 1e3) * _A_webs_full_m2
        
        # Rods mass
        _A_rod_cs     = _rod_w * _rod_h
        _vol_rods_m3  = _n_rods_total * _span_full_m * _A_rod_cs
        _mass_rods_kg = _vol_rods_m3 * _rho_rod
        
        # Totals
        _mass_total_kg = _mass_skins_kg + _mass_core_kg + _mass_webs_kg + _mass_rods_kg
        _g0 = 9.80665
        _W_total_N  = _mass_total_kg * _g0
        _W_total_lb = _W_total_N / 4.4482216153
        
        _rows = [
            ("Planform area (full) [m^2]", _A_planform_full_m2),
            ("Skins mass [kg]",            _mass_skins_kg),
            ("Core mass [kg]",             _mass_core_kg),
            ("Webs mass [kg]",             _mass_webs_kg),
            ("Rods mass [kg]",             _mass_rods_kg),
            ("TOTAL mass [kg]",            _mass_total_kg),
            ("TOTAL weight [N]",           _W_total_N),
            ("TOTAL weight [lb]",          _W_total_lb),
            ("Half-wing mass [kg]",        0.5 * _mass_total_kg),
            ("Skin gsm/ply (impregnated)", _skin_gsm_per_ply),
            ("Web dry gsm label",          _WEB_DRY_GSM),
            ("AW is total for biax?",      _AW_IS_TOTAL_FOR_BIAX),
            ("Vf (web)",                   _Vf),
            ("Web cured gsm/layer",        _web_cured_gsm_per_layer),
            ("Web layers (assumed)",       _n_web_layers),
            ("Rod count (top+bot)",        _n_rods_total),
            ("Rod CS [mm x mm]",           f"{_rod_w*1e3:.1f} x {_rod_h*1e3:.1f}"),
            ("Core t [mm], ρ [kg/m^3]",    f"{_t_core*1e3:.2f}, {_rho_core:.1f}"),
        ]
        _df_weight = _pd.DataFrame(_rows, columns=["Parameter", "Value"])
        _df_weight.to_excel(writer, sheet_name="Weight_Summary", index=False)
        # =================== END: Weight_Summary (build from existing results) ===================
        for name, df in results.items():
            if name in ["Final_Summary", "Equations_Appendix", "Inputs", "Materials", "All_Sweeps", "Summaries_All_Sweeps", "Airfoil_xy"]:
                continue
            sheet_name = name[:31]
            try:
                df.to_excel(writer, sheet_name=sheet_name, index=False)
            except Exception as e:
                print(f"Warning: failed to write sheet {sheet_name}: {e}")

    print(f"Excel written successfully: {inp_base.out_excel}")

    if inp_base.make_plots:
        print(f"\n{'='*70}")
        print(f"Generating plots...")
        print(f"{'='*70}")
        try:
            plot_all_results(results, inp_base)
        except Exception as e:
            print(f"Warning: plotting failed: {e}")

    print(f"\n{'='*70}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'='*70}")

    if not results["Final_Summary"].empty:
        worst_sweep = results["Final_Summary"]["Worst Sweep (deg)"].iloc[0]
        min_fos = results["Final_Summary"]["Min FoS @Design G"].iloc[0]
        print(f"Worst case sweep angle: {worst_sweep:.1f}°")
        print(f"Minimum Factor of Safety: {min_fos:.3f}")

    print(f"\n{'='*70}\n")

    # Use results dataframes instead of local arrays (avoid NameError)
    try:
        all_sweeps = results.get("All_Sweeps", pd.DataFrame())
        if not all_sweeps.empty:
            tau_skin_nontriv = int((all_sweeps["tau_skin_tor (Pa)"].abs() > 1e-6).sum()) if "tau_skin_tor (Pa)" in all_sweeps.columns else int((all_sweeps["tau_skin_tor"].abs() > 1e-6).sum())
            tau_core_nontriv = int((all_sweeps["tau_core (Pa)"].abs() > 1e-6).sum()) if "tau_core (Pa)" in all_sweeps.columns else int((all_sweeps["tau_core"].abs() > 1e-6).sum())
            print("diag: tau_skin_tor nontrivial:", tau_skin_nontriv, "of", len(all_sweeps))
            print("diag: tau_core nontrivial:", tau_core_nontriv, "of", len(all_sweeps))
            # print a representative lift value for sweep 0 if present
            if "lift_per_unit_span (N/m)" in all_sweeps.columns:
                print("diag: lift_per_unit_span[0] =", float(all_sweeps["lift_per_unit_span (N/m)"].iloc[0]))
        else:
            print("diag: All_Sweeps dataframe empty; no final diagnostics available")
    except Exception as _e:
        print("diag: final diagnostics failed:", _e)

if __name__ == "__main__":
    main()
