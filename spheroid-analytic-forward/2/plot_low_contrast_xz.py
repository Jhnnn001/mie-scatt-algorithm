"""Plot five xz incidences: Rytov/exact magnitude and incident-relative phase."""
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.io import loadmat

out = Path(__file__).resolve().parent
d = loadmat(out / "low_contrast_xz_fields.mat", squeeze_me=True, struct_as_record=False)
x, z, angles = d["x"], d["z"], d["theta_deg"]
p = d["p"]
assert np.allclose(angles, [0, 5, 10, 15, 20])
assert np.isclose(p.n_p, p.n_m + .01*(1.365-p.n_m), rtol=0, atol=1e-14)
assert p.a == 3 and p.b == 5 and getattr(p, "lambda") == .532
fields = [d["rytovEx"], d["exactEx"]]
incident = d["incidentEx"]
shape = (len(x)*len(z), len(angles))
assert all(f.shape == shape and np.isfinite(f).all() for f in [*fields, incident])
assert np.allclose(np.abs(incident), np.cos(np.deg2rad(angles))[None, :])
magnitudes = [np.abs(f) for f in fields]
phases = [np.angle(f/incident) for f in fields]
# No wrapping or unwrapping ambiguity in this weak-contrast phase range.
assert max(np.abs(v).max() for v in phases) < .1
phase_limit = np.ceil(max(np.abs(v).max() for v in phases)*1000)/1000
phase_norm = Normalize(-phase_limit, phase_limit)
plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42, "axes.titlesize": 12})
fig = plt.figure(figsize=(17.8, 18.4), layout="constrained")
gs = fig.add_gridspec(5, 6, width_ratios=[1, 1, .045, 1, 1, .045])
extent = [x[0]-(x[1]-x[0])/2, x[-1]+(x[1]-x[0])/2,
          z[0]-(z[1]-z[0])/2, z[-1]+(z[1]-z[0])/2]
titles = ["Rytov\nmagnitude |Ex|", "Maxwell exact\nmagnitude |Ex|",
          "Rytov\nphase relative to incidence", "Maxwell exact\nphase relative to incidence"]
records = []
for row, theta in enumerate(angles):
    values = [v[:, row] for v in magnitudes]
    low = min(v.min() for v in values)
    high = max(v.max() for v in values)
    pad = .02*(high-low)
    amp_norm = Normalize(np.floor((low-pad)*10000)/10000,
                         np.ceil((high+pad)*10000)/10000)
    for col, grid_col in enumerate([0, 1, 3, 4]):
        ax = fig.add_subplot(gs[row, grid_col])
        is_phase = col >= 2
        model = col % 2
        values = phases[model] if is_phase else magnitudes[model]
        arr = values[:, row].reshape((len(z), len(x)), order="F")
        im = ax.imshow(arr, origin="lower", extent=extent, aspect="equal", interpolation="nearest",
                       cmap="RdBu_r" if is_phase else "viridis",
                       norm=phase_norm if is_phase else amp_norm)
        color = "black" if is_phase else "white"
        ax.add_patch(Ellipse((0, 0), 2*p.b, 2*p.a, fill=False, edgecolor=color, lw=.8))
        rad = np.deg2rad(theta)
        start = np.array([-6., -3.8]); end = start + 1.4*np.array([np.sin(rad), np.cos(rad)])
        ax.annotate("", xy=end, xytext=start,
                    arrowprops={"arrowstyle": "->", "color": color, "lw": 1.2})
        ax.set(xlim=(-8, 8), ylim=(-5, 9), xticks=[-8, -4, 0, 4, 8], yticks=[-4, 0, 4, 8])
        if row == 0:
            ax.set_title(titles[col], pad=9)
        if col == 0:
            ax.set_ylabel(f"θ = {theta:.0f}°\nz (μm)", fontsize=12)
        else:
            ax.tick_params(labelleft=False)
        if row == 4:
            ax.set_xlabel("x (μm)")
        else:
            ax.tick_params(labelbottom=False)
        if col == 1:
            cb = fig.colorbar(im, cax=fig.add_subplot(gs[row, 2]))
            cb.locator = MaxNLocator(4); cb.update_ticks()
            cb.ax.tick_params(labelsize=8)
        if is_phase:
            phase_image = im
    for model, label in enumerate(["Rytov", "Maxwell exact"]):
        records.append({"theta_deg": theta, "model": label,
                        "incident_Ex_magnitude": float(np.cos(np.deg2rad(theta))),
                        "magnitude_min": magnitudes[model][:, row].min(),
                        "magnitude_max": magnitudes[model][:, row].max(),
                        "relative_phase_min_rad": phases[model][:, row].min(),
                        "relative_phase_max_rad": phases[model][:, row].max(),
                        "magnitude_color_min": amp_norm.vmin, "magnitude_color_max": amp_norm.vmax,
                        "phase_color_min_rad": -phase_limit, "phase_color_max_rad": phase_limit})
cb = fig.colorbar(phase_image, cax=fig.add_subplot(gs[:, 5]))
cb.set_label("arg(Ex / Ex,inc) (rad)")
fig.suptitle("Low-contrast oblate spheroid: Δn = 0.01 × current contrast\n"
             "nₚ = 1.335677718660; nₘ = 1.335381534; λ₀ = 532 nm; a(z) = 3 μm, b(x,y) = 5 μm\n"
             "y = 0; raw total fields; |E₀| = 1; magnitude scale shared within each row; one phase scale for all rows",
             fontsize=13)
for suffix in ["png", "pdf"]:
    fig.savefig(out / f"low_contrast_xz_5x4.{suffix}", dpi=200, bbox_inches="tight", pad_inches=.15)
plt.close(fig)
with (out / "low_contrast_xz_plot_stats.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(records[0]))
    writer.writeheader(); writer.writerows(records)
np.savez_compressed(out / "low_contrast_xz_plot_data.npz", x=x, z=z, theta_deg=angles,
                    magnitude_rytov=magnitudes[0], magnitude_exact=magnitudes[1],
                    phase_rytov=phases[0], phase_exact=phases[1],
                    absolute_phase_rytov=np.angle(fields[0]), absolute_phase_exact=np.angle(fields[1]))
print(f"LOW_CONTRAST_XZ_PLOT_PASS: 5 x 4 maps; phase scale ±{phase_limit:g} rad")
for row in records:
    print(row)
