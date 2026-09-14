"""Validate the retained numerical results and regenerate the application figures."""
import csv
from pathlib import Path

import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE = Path(__file__).resolve().parent
SRC = HERE
FIGURES = HERE / "results"
FIGURES.mkdir(exist_ok=True)
FIELDS = SRC / "forward" / "exp" / "oblate_xz_fields.mat"
EXACT = SRC / "forward" / "exp" / "oblate_xz_exact.mat"
SWEEP = SRC / "spheroid-analytic-forward"
RESULTS = SRC / "inverse" / "results"
RECON = [("padding8", "direct"), ("padding8", "gp"), ("tv_strict", "tv")]

plt.rcParams.update({"font.size": 8, "savefig.dpi": 200, "axes.linewidth": 0.6})


def fig1_forward():
    d = loadmat(FIELDS, simplify_cells=True)
    x, z, theta = d["x"], d["z"], d["theta_deg"]
    shape = (len(x), len(z))
    inside = d["inside"].reshape(shape).astype(float)
    xs, zs = (x >= -6.4) & (x <= 6.4), (z >= -4) & (z <= 8)
    extent = [x[xs][0], x[xs][-1], z[zs][0], z[zs][-1]]

    angles = [int(np.argmin(abs(theta - 0))), int(np.argmin(abs(theta - 15)))]
    models = ("exact", "bornEx", "rytovEx")
    amp, pha = [], []
    for k in angles:
        u0 = d["u0"][:, k].reshape(shape)
        for m in models:
            field = d[m][:, k].reshape(shape)
            amp.append(np.abs(field)[np.ix_(xs, zs)].T)
            pha.append(np.angle(field / u0)[np.ix_(xs, zs)].T)
    amax = max(a.max() for a in amp)
    pmin, pmax = min(p.min() for p in pha), max(p.max() for p in pha)

    fig, axes = plt.subplots(2, 6, figsize=(7.0, 2.85), layout="constrained")
    for j in range(6):
        for i, (data, vmin, vmax, cmap) in enumerate(
                [(amp[j], 0, amax, "viridis"), (pha[j], pmin, pmax, "magma")]):
            ax = axes[i, j]
            im = ax.imshow(data, origin="lower", extent=extent, aspect="equal",
                           vmin=vmin, vmax=vmax, cmap=cmap)
            ax.contour(x[xs], z[zs], inside[np.ix_(xs, zs)].T, [0.5],
                       colors="white", linewidths=0.5)
            ax.set_xticks([-5, 0, 5]); ax.set_yticks([0, 5])
            if j: ax.set_yticklabels([])
            if i == 0: ax.set_xticklabels([])
            if i == 0 and j == 0: ax.set_ylabel("z (μm)")
            if i == 1 and j == 0: ax.set_ylabel("z (μm)")
            if i == 1: ax.set_xlabel("x (μm)")
            if i == 0:
                ax.set_title(f"{('Exact', 'Born', 'Rytov')[j % 3]}\nθ = {(0, 15)[j // 3]}°", fontsize=8)
            if j == 5:
                fig.colorbar(im, ax=axes[i, :], shrink=0.9, pad=0.01,
                             label="|E$_x$| / |E$_0$|" if i == 0 else "phase (rad)")
    fig.savefig(FIGURES / "fig1-forward.png")
    plt.close(fig)


def fig2_Q():
    """Plot the omitted logarithmic-gradient term at normal incidence."""
    d = loadmat(EXACT, simplify_cells=True)
    x, z = d["x"], d["z"]
    shape = (len(x), len(z))
    k = int(np.argmin(abs(d["theta_deg"])))
    inside = d["inside"].reshape(shape).astype(float)
    q = np.abs(d["Q"][:, k]).reshape(shape) / float(d["f0"])
    q[inside == 0] = np.nan
    xs, zs = (x >= -5.3) & (x <= 5.3), (z >= -3.3) & (z <= 3.3)
    extent = [x[xs][0], x[xs][-1], z[zs][0], z[zs][-1]]

    fig, ax = plt.subplots(figsize=(3.0, 1.95), layout="constrained")
    cmap = matplotlib.colormaps["magma"].copy(); cmap.set_bad("white")
    im = ax.imshow(q[np.ix_(xs, zs)].T, origin="lower", extent=extent, aspect="equal",
                   norm=LogNorm(vmin=0.01, vmax=20), cmap=cmap)
    ax.contour(x[xs], z[zs], inside[np.ix_(xs, zs)].T, [0.5], colors="black", linewidths=0.5)
    ax.set_xticks([-5, 0, 5]); ax.set_yticks([-3, 0, 3])
    ax.set_xlabel("x (μm)"); ax.set_ylabel("z (μm)")
    fig.colorbar(im, ax=ax, shrink=0.85, pad=0.03, label="|Q| / |f|")
    fig.savefig(FIGURES / "fig2-Q.png")
    plt.close(fig)


def fig3_size_contrast():
    """Mean Born/Rytov error over 81 directions vs. size factor (1/) and Δn factor (2/)."""
    fig, axes = plt.subplots(1, 2, figsize=(5.6, 2.4), layout="constrained")
    for ax, sub, xl in zip(axes, ("1", "2"), ("size factor", "Δn factor")):
        fac = {r["case_id"]: float(r["factor"])
               for r in csv.DictReader(open(SWEEP / sub / "parameters.csv"))}
        err = {}
        for r in csv.DictReader(open(SWEEP / sub / "summary.csv")):
            err.setdefault(r["metric"], {})[fac[r["case_id"]]] = 100 * float(r["mean_error"])
        for metric, style, lab in (("born_scalar", "o-", "Born"), ("rytov_scalar", "s-", "Rytov")):
            xs = sorted(err[metric])
            ax.loglog(xs, [err[metric][v] for v in xs], style, ms=3.5, lw=1, label=lab)
        ax.set_xlabel(xl); ax.set_ylabel("relative error (%)")
        ax.grid(True, which="both", lw=0.3, alpha=0.5)
    axes[1].legend(loc="upper left", frameon=False, fontsize=7)
    fig.savefig(FIGURES / "fig3-size-contrast.png")
    plt.close(fig)


def fig4_reconstruction():
    g = loadmat(RESULTS / "padding8_geometry.mat", simplify_cells=True)
    truth, nm = g["truth"], float(g["nm"])
    delta = float(g["np"]) - nm
    x, y, z = g["x"], g["y"], g["z"]
    ix, iy, iz = truth.shape[0] // 2, truth.shape[1] // 2, truth.shape[2] // 2
    dx = x[1] - x[0]
    ext_xy = [x[0] - dx/2, x[-1] + dx/2, y[0] - dx/2, y[-1] + dx/2]
    ext_xz = [x[0] - dx/2, x[-1] + dx/2, z[0] - dx/2, z[-1] + dx/2]
    vols = []
    for label, method in RECON:
        assert np.array_equal(loadmat(RESULTS / f"{label}_geometry.mat",
                                      simplify_cells=True)["truth"], truth)
        vols.append((loadmat(RESULTS / f"{label}_{method}.mat", simplify_cells=True)["R"]["n"] - nm) / delta)
    vmin = min(0.0, min(float(v.min()) for v in vols))
    vmax = max(float(v.max()) for v in vols)
    tb = (truth - nm) / delta

    fig, axes = plt.subplots(2, 3, figsize=(6.0, 3.55), layout="constrained",
                             gridspec_kw={"height_ratios": [len(y), len(z)]})
    for j, v in enumerate(vols):
        for i, (sl, tsl, ext, yl) in enumerate([
                (v[:, :, iz].T, tb[:, :, iz].T, ext_xy, "y (μm)"),
                (v[:, iy, :].T, tb[:, iy, :].T, ext_xz, "z (μm)")]):
            ax = axes[i, j]
            im = ax.imshow(sl, origin="lower", extent=ext, aspect="equal",
                           vmin=vmin, vmax=vmax, cmap="viridis")
            ax.contour(x, y if i == 0 else z, tsl, [0.5], colors="white", linewidths=0.6)
            if j: ax.set_yticklabels([])
            else: ax.set_ylabel(yl)
            if i == 0:
                ax.set_xticklabels([])
                ax.set_title(("Direct", "GP", "TV")[j])
            else: ax.set_xlabel("x (μm)")
    fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02, label="(n − n$_m$) / Δn")
    fig.savefig(FIGURES / "fig4-reconstruction.png")
    plt.close(fig)



def verify_results():
    """Catch wrong case selection, corrupted volumes, and unreported nonconvergence."""
    for label in ("padding8", "tv_strict", "linear_control"):
        g = loadmat(RESULTS / f"{label}_geometry.mat", simplify_cells=True)
        truth, nm = g["truth"], float(g["nm"])
        with (RESULTS / f"{label}_metrics.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            method = row["method"]
            r = loadmat(RESULTS / f"{label}_{method}.mat", simplify_cells=True)["R"]
            n = r["n"]
            assert n.shape == truth.shape and np.isfinite(n).all()
            error = np.linalg.norm(n - truth) / np.linalg.norm(truth - nm)
            assert np.isclose(error, float(row["relative_ri_error"]), rtol=1e-9)
            assert bool(r["info"]["converged"]) == bool(int(row["converged"]))
            if method != "direct":
                assert n.min() >= nm - 1e-13
            print(f"{label}/{method}: RI error {100*error:.4f}%, converged={bool(r['info']['converged'])}")
    d = loadmat(EXACT, simplify_cells=True)
    q = np.abs(d["Q"][:, 0]) / float(d["f0"])
    q[~d["inside"].astype(bool)] = -1
    ix, iz = np.unravel_index(q.argmax(), (len(d["x"]), len(d["z"])))
    assert np.isclose(q.max(), 14.32723618, rtol=1e-6)
    assert np.isclose(d["x"][ix], 0) and np.isclose(d["z"][iz], -0.2)
    print(f"Interior |Q|/|f| maximum: {q.max():.6f} at (x,z)=(0,-0.2) um")
    for sub in ("1", "2"):
        with (SWEEP / sub / "summary.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        values = {(int(r["case_id"]), r["metric"]): float(r["mean_error"]) for r in rows}
        assert np.isclose(100*values[1, "born_co"], 91.31903, atol=1e-4)
        assert np.isclose(100*values[1, "rytov_co"], 10.50928, atol=1e-4)
        assert values[3, "born_scalar"] < .01 and values[3, "rytov_scalar"] < .01
    with (SRC / "forward/tests/rytov_analysis.csv").open() as stream:
        row = next(csv.DictReader(stream))
    assert np.isclose(100*float(row["born_error"]), 197.55961, atol=1e-4)
    assert np.isclose(100*float(row["corrected_rytov_error"]), 25.73520, atol=1e-4)


if __name__ == "__main__":
    verify_results()
    fig1_forward()
    fig2_Q()
    fig3_size_contrast()
    fig4_reconstruction()
    print("RESULTS_AND_FIGURES_PASS")
