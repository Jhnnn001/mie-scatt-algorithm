"""Render the computed total Ex fields with one physical scale and color scale."""
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
import numpy as np
from scipy.io import loadmat

out = Path(__file__).resolve().parent
d = loadmat(out / "oblate_xz_fields.mat", squeeze_me=True)
x, z, angles = d["x"], d["z"], d["theta_deg"]
fields = [d[name] for name in ("exact", "bornEx", "rytovEx")]
assert np.allclose(angles, [0, 5, 10, 15])
assert all(a.shape == (len(x) * len(z), 4) and np.isfinite(a).all() for a in fields)
vmax = np.ceil(max(np.abs(a).max() for a in fields) * 10) / 10
norm = Normalize(0, vmax)
plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
fig, axs = plt.subplots(4, 3, figsize=(11.2, 14), sharex=True, sharey=True,
                        layout="constrained")
extent = [x[0] - np.diff(x)[0] / 2, x[-1] + np.diff(x)[0] / 2,
          z[0] - np.diff(z)[0] / 2, z[-1] + np.diff(z)[0] / 2]
for row, theta in enumerate(angles):
    for col, (name, field) in enumerate(zip(("Exact", "Born", "Rytov"), fields)):
        ax = axs[row, col]
        amplitude = np.abs(field[:, row]).reshape((len(z), len(x)), order="F")
        im = ax.imshow(amplitude, origin="lower", extent=extent, cmap="viridis",
                       norm=norm, interpolation="nearest", aspect="equal")
        ax.add_patch(Ellipse((0, 0), 10, 6, fill=False, edgecolor="white", lw=0.8))
        ax.set(xlim=(-8, 8), ylim=(-5, 9), xticks=[-8, -4, 0, 4, 8], yticks=[-4, 0, 4, 8])
        ax.tick_params(length=3)
        if row == 0:
            ax.set_title(name, fontsize=15, pad=9)
        if col == 0:
            ax.set_ylabel(r"$z$ ($\mu$m)")
            ax.text(0.04, 0.94, rf"$\theta={theta:.0f}^\circ$", transform=ax.transAxes,
                    va="top", color="white", bbox={"facecolor": "black", "alpha": .45, "pad": 3})
        if row == 3:
            ax.set_xlabel(r"$x$ ($\mu$m)")
        rad = np.deg2rad(theta)
        start = np.array([-6., -3.8])
        end = start + 1.4 * np.array([np.sin(rad), np.cos(rad)])
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "color": "white", "lw": 1.3})
fig.colorbar(im, ax=axs, shrink=.87, pad=.025, label=r"Total $|E_x|/|E_0|$")
fig.suptitle("Total electric-field amplitude in the xz plane (y = 0)\n"
             r"$a_z=3\,\mu$m, $b_{xy}=5\,\mu$m; $n_p=1.365$, $n_m=1.335381534$, $\lambda_0=532$ nm"
             "\nRaw fields; white outline: spheroid boundary; arrows: incident direction",
             fontsize=12)
for suffix in ("png", "pdf"):
    fig.savefig(out / f"oblate_Ex_xz_4x3.{suffix}", dpi=220)
print(f"4 x 3 panels saved; shared color limits [0, {vmax:g}].")
with (out / "oblate_xz_plot_errors.csv").open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["theta_deg", "model", "region", "complex_relative_L2", "amplitude_relative_L2"])
    scalar_ex = d["scalar"] * np.cos(np.deg2rad(angles))
    for row, theta in enumerate(angles):
        for name, field in zip(("Born", "Rytov", "Scalar exact"), [*fields[1:], scalar_ex]):
            for region, mask in (("whole xz window", np.ones(len(field), dtype=bool)),
                                 ("spheroid interior", d["inside"].astype(bool))):
                target, estimate = fields[0][mask, row], field[mask, row]
                denom = np.linalg.norm(target)
                writer.writerow([theta, name, region, np.linalg.norm(estimate - target) / denom,
                                 np.linalg.norm(np.abs(estimate) - np.abs(target)) / denom])
