"""Recompute saved RI metrics independently and plot the inverse experiments."""
from pathlib import Path
import csv
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / "results"
METHODS = ("direct", "gp", "tv")
COLORS = ("#4477AA", "#EE6677", "#228833")


def main():
    plt.rcParams.update({"font.size": 10, "savefig.dpi": 180})
    checked = []
    tv_profiles = {}
    for path in sorted(OUT.glob("*_metrics.csv")):
        label = path.name.removesuffix("_metrics.csv")
        g = loadmat(OUT / f"{label}_geometry.mat", simplify_cells=True)
        truth, nm, npart = g["truth"], float(g["nm"]), float(g["np"])
        delta = npart - nm
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        results = {}
        for row in rows:
            method = row["method"]
            r = loadmat(OUT / f"{label}_{method}.mat", simplify_cells=True)["R"]
            n = r["n"]
            ix, iy = n.shape[0]//2, n.shape[1]//2
            assert n.shape == truth.shape and np.isfinite(n).all()
            inside = truth > nm
            denominator = np.linalg.norm(truth - nm)
            values = {
                "relative_ri_error": np.linalg.norm(n - truth) / denominator,
                "object_error": np.linalg.norm(n[inside] - npart) / denominator,
                "background_error": np.linalg.norm(n[~inside] - nm) / denominator,
                "mean_object_n": n[inside].mean(),
                "relative_index_bias": (n[inside].mean() - npart) / delta,
            }
            for key, value in values.items():
                assert np.isclose(value, float(row[key]), rtol=1e-8, atol=1e-10), (label, method, key)
            if "axial_span_censored" in row:
                profile = n[n.shape[0]//2, n.shape[1]//2, :] > nm + .5*delta
                assert bool(profile[0] or profile[-1]) == bool(int(row["axial_span_censored"]))
            if method != "direct":
                assert n.min() >= nm - 1e-13
            results[method] = r
            checked.append((label, method, values["relative_ri_error"]))
            if method == "tv" and row["source"] == "exact":
                tv_profiles[label] = (g["z"], (n[ix, iy, :]-nm)/delta,
                                      (truth[ix, iy, :]-nm)/delta, values["relative_ri_error"])
            if label == "tv_strict" and method == "tv":
                contrast = (n-nm)/delta
                x, z = g["x"], g["z"]
                dx = x[1]-x[0]
                fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
                image = axes[0].imshow(contrast[:, iy, :].T, origin="lower",
                    extent=[x[0]-dx/2, x[-1]+dx/2, z[0]-dx/2, z[-1]+dx/2],
                    vmin=min(0, contrast.min()), vmax=max(1.05, contrast.max()), cmap="viridis")
                axes[0].contour(x, z, (truth[:, iy, :].T-nm)/delta, [.5], colors="white", linewidths=.8)
                axes[0].set(xlabel="x (μm)", ylabel="z (μm)", title="White contour: true boundary")
                fig.colorbar(image, ax=axes[0], label="RI contrast / true Δn")
                axes[1].plot(z, (truth[ix, iy, :]-nm)/delta, "k--", label="Truth")
                axes[1].plot(z, contrast[ix, iy, :], label="Strict refined TV")
                axes[1].set(xlabel="z (μm)", ylabel="Central RI contrast / true Δn")
                axes[1].legend(); axes[1].grid(alpha=.2)
                fig.suptitle(f"Refined TV: RI error {100*values['relative_ri_error']:.2f}%")
                fig.savefig(OUT / "tv_strict_reconstruction.png")
                plt.close(fig)
        if not all(m in results for m in METHODS):
            continue
        x, y, z = g["x"], g["y"], g["z"]
        dx = x[1] - x[0]
        extent_xy = [x[0]-dx/2, x[-1]+dx/2, y[0]-dx/2, y[-1]+dx/2]
        extent_xz = [x[0]-dx/2, x[-1]+dx/2, z[0]-dx/2, z[-1]+dx/2]
        volumes = [truth] + [results[m]["n"] for m in METHODS]
        minimum = min(float((v.min()-nm)/delta) for v in volumes)
        maximum = max(float((v.max()-nm)/delta) for v in volumes)
        titles = ["Truth"] + [f"{m.upper()}: RI error {100*results[m]['metric']['relative_ri_error']:.2f}%" for m in METHODS]
        fig, axes = plt.subplots(2, 4, figsize=(14, 6.2), layout="constrained")
        for j, (volume, title) in enumerate(zip(volumes, titles)):
            contrast = (volume - nm) / delta
            for ax, plane, extent in (
                (axes[0, j], contrast[:, :, contrast.shape[2]//2].T, extent_xy),
                (axes[1, j], contrast[:, contrast.shape[1]//2, :].T, extent_xz),
            ):
                image = ax.imshow(plane, origin="lower", extent=extent, cmap="viridis", vmin=minimum, vmax=maximum)
                ax.set_xlabel("x (μm)")
            axes[0, j].set_title(title)
        axes[0, 0].set_ylabel("y (μm)")
        axes[1, 0].set_ylabel("z (μm)")
        fig.colorbar(image, ax=axes, shrink=.86, label="(n − nₘ) / true Δn")
        fig.suptitle(f"{label}: NA illumination 0.5, detector 1.0; Δn × 0.01")
        fig.savefig(OUT / f"{label}_reconstruction.png")
        plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
        ix, iy = truth.shape[0]//2, truth.shape[1]//2
        axes[0].plot(z, (truth[ix, iy, :] - nm)/delta, "k--", label="Truth")
        for method, color in zip(METHODS, COLORS):
            r = results[method]
            axes[0].plot(z, (r["n"][ix, iy, :] - nm)/delta, color=color, label=method.upper())
            history = np.atleast_2d(r["info"]["history"])
            if method == "direct":
                axes[1].axhline(history[-1, 3], color=color, ls="--", label=method.upper())
            else:
                axes[1].semilogy(np.arange(len(history))+1, history[:, 3], color=color, label=method.upper())
        axes[0].axhline(.5, color=".7", lw=.8)
        axes[0].set(xlabel="z (μm)", ylabel="Central-line RI contrast / true Δn")
        axes[1].set(xlabel="Iteration (TV: cumulative inner iterations)", ylabel="Weighted Ewald relative residual")
        for ax in axes:
            ax.legend(); ax.grid(alpha=.2)
        fig.suptitle(label)
        fig.savefig(OUT / f"{label}_profiles_convergence.png")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(6, 4), layout="constrained")
        ax.pcolormesh(g["axis_qx"], g["axis_qz"], g["coverage_xz"].T, cmap="Greys", shading="nearest")
        ax.set(xlabel="qₓ (μm⁻¹)", ylabel="q_z (μm⁻¹)", xlim=(-19, 19), ylim=(-8, 8),
               title="Interpolation support + inferred Hermitian symmetry")
        fig.savefig(OUT / f"{label}_missing_cone.png")
        plt.close(fig)
    assert checked, "No completed reconstruction files."
    groups = [("TV weight", ("tv_half", "tv_long", "tv_double")),
              ("Inner stopping and penalties", ("exact", "tv_long", "tv_penalty")),
              ("Padding and tolerance", ("tv_penalty", "tv_refined", "tv_strict"))]
    if all(label in tv_profiles for _, labels in groups for label in labels):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), layout="constrained")
        for ax, (title, labels) in zip(axes, groups):
            z, _, truth_profile, _ = tv_profiles[labels[0]]
            ax.plot(z, truth_profile, "k--", label="Truth")
            for label in labels:
                z, profile, _, error = tv_profiles[label]
                ax.plot(z, profile, label=f"{label}: {100*error:.2f}%")
            ax.set(title=title, xlabel="z (μm)", ylabel="RI contrast / true Δn")
            ax.legend(fontsize=8); ax.grid(alpha=.2)
        fig.suptitle("Same physical data; legends give whole-volume RI error")
        fig.savefig(OUT / "tv_sensitivity.png")
        plt.close(fig)
    print(f"INVERSE_ARTIFACT_CHECK_PASS: {len(checked)} reconstructions")
    for label, method, error in checked:
        print(f"{label:20s} {method:6s} RI error {100*error:.6f}%")


if __name__ == "__main__":
    main()
