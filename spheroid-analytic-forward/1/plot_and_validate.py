"""Render both experiments and independently check their saved CSV summaries.

Run with the existing NumPy/Matplotlib environment; no MATLAB execution here.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
COLORS = ["#174a7e", "#d16a17", "#29835d"]
METRICS = ["born_scalar", "rytov_scalar", "born_co", "rytov_co",
           "born_full", "rytov_full", "scalar_co", "scalar_full",
           "vector_born_co", "vector_born_full"]


def read(path):
    def number_or_text(value):
        try:
            return float(value)
        except ValueError:
            return value
    with path.open(newline="") as stream:
        return [{key: number_or_text(value)
                 for key, value in row.items()} for row in csv.DictReader(stream)]


def weights(rows):
    na = np.array([r["illumination_NA"] for r in rows])
    rings = np.unique(na)
    edges = np.r_[0, (rings[:-1] + rings[1:]) / 2, .5]
    w = np.zeros(len(rows))
    for q, area in zip(rings, np.diff(edges**2) / .5**2):
        mask = na == q
        w[mask] = area / mask.sum()
    assert abs(w.sum() - 1) < 1e-14
    return w


def main():
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180})
    for project in ("1", "2"):
        out = ROOT / project
        parameters = read(out / "parameters.csv")
        angles = read(out / "angles.csv")
        summary = read(out / "summary.csv")
        assert len(parameters) == 3 and len(angles) == 243 and len(summary) == 30
        lookup = {(int(r["case_id"]), r["metric"]): r for r in summary}
        signals = []
        for case in (1, 2, 3):
            rows = [r for r in angles if r["case_id"] == case]
            assert len(rows) == 81
            assert len({(r["illumination_NA"], r["phi_deg"]) for r in rows}) == 81
            assert all(np.isfinite(v) for r in rows for v in r.values())
            w = weights(rows)
            signals.append(np.sqrt(np.dot(w, [r["scalar_norm2"] for r in rows])))
            for metric in METRICS:
                values = np.array([r[metric] for r in rows])
                norm_key = ("scalar_norm2" if metric.endswith("scalar") else
                            "co_norm2" if metric.endswith("co") else "vector_norm2")
                norm2 = np.array([r[norm_key] for r in rows])
                actual = lookup[case, metric]
                expected = {"mean_error": np.dot(w, values),
                            "pooled_L2": np.sqrt(np.dot(w, values**2 * norm2) / np.dot(w, norm2)),
                            "minimum": values.min(), "maximum": values.max(),
                            "fraction_below_1pct": w[values < .01].sum(),
                            "fraction_below_5pct": w[values < .05].sum()}
                for key, value in expected.items():
                    assert np.isclose(actual[key], value, rtol=1e-11, atol=1e-13), (project, case, metric, key)
            checks = read(out / f"case_{case}_detector_convergence.csv")
            assert len(checks) == 7
            assert max(r["sampling_field_gap"] for r in checks) < 1e-4
            assert max(r["window_field_gap"] for r in checks) < .002
        ref = read(out / "reference_validation.csv")
        assert len(ref) == 3
        assert max(r["cutoff_vector"] for r in ref) < 1e-6
        assert max(r["cutoff_scalar"] for r in ref) < 1e-6
        assert max(r["public_far_gap"] for r in ref) < 1e-5
        assert "PARAMETER_EXPERIMENT_PASS" in (out / "run.log").read_text()

        factors = np.array([p["factor"] for p in parameters])
        xname = "Linear size / current size" if project == "1" else "Index contrast / current contrast"
        fig, axes = plt.subplots(1, 3, figsize=(14, 4.3), layout="constrained")
        panels = [(axes[0], ["born_scalar", "rytov_scalar", "vector_born_full"],
                   ["Scalar Born / scalar exact", "Scalar Rytov / scalar exact", "Vector Born / Maxwell"],
                   "Approximation error"),
                  (axes[1], ["born_co", "rytov_co", "scalar_co"],
                   ["Scalar Born / Maxwell co", "Scalar Rytov / Maxwell co", "Scalar exact / Maxwell co"],
                   "Error against Maxwell co-polarization")]
        for ax, metrics, labels, title in panels:
            for color, metric, label, style in zip(COLORS, metrics, labels, ["o-", "s--", "^:"]):
                y = 100*np.array([lookup[case, metric]["mean_error"] for case in (1, 2, 3)])
                ax.loglog(factors, y, style, color=color, label=label, markerfacecolor="white")
            for value in (1, 5):
                ax.axhline(value, color=".7", lw=.7, ls=":")
            ax.set(xlabel=xname, ylabel="Pupil-weighted mean relative L2 error (%)", title=title)
            ax.grid(True, which="both", alpha=.15)
            ax.legend(fontsize=8, loc="best")
        axes[2].loglog(factors, np.array(signals)/signals[0], "o-", color=COLORS[0])
        axes[2].set(xlabel=xname, ylabel="Detected scalar scattered-field norm / current",
                    title="Scattering strength (fixed illumination)")
        axes[2].grid(True, which="both", alpha=.15)
        fig.suptitle("Oblate size sweep" if project == "1" else "Oblate index-contrast sweep")
        for suffix in ("png", "pdf"):
            fig.savefig(out / f"error_vs_parameter.{suffix}", bbox_inches="tight", pad_inches=.15)
        plt.close(fig)

        fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
        for ax, metric, title in zip(axes.flat,
                ["born_scalar", "rytov_scalar", "rytov_co", "vector_born_full"],
                ["Scalar Born / scalar exact", "Scalar Rytov / scalar exact",
                 "Scalar Rytov / Maxwell co", "Vector Born / Maxwell full field"]):
            for case, color, factor in zip((1, 2, 3), COLORS, factors):
                rows = [r for r in angles if r["case_id"] == case]
                theta = sorted({r["theta_deg"] for r in rows})
                groups = [[100*r[metric] for r in rows if r["theta_deg"] == angle] for angle in theta]
                ax.semilogy(theta, [np.mean(v) for v in groups], "o-", color=color, label=f"factor = {factor:g}")
                ax.fill_between(theta, [min(v) for v in groups], [max(v) for v in groups], color=color, alpha=.15)
            ax.set(xlabel="Incidence angle from z (degrees)", ylabel="Relative L2 error (%)", title=title)
            ax.grid(True, which="both", alpha=.15)
            ax.legend(fontsize=9, loc="best" if metric == "rytov_co" else "center right")
        fig.suptitle("Lines: azimuth mean; shading: azimuth range. Illumination NA <= 0.5; detector NA = 1.0")
        for suffix in ("png", "pdf"):
            fig.savefig(out / f"error_vs_angle.{suffix}", bbox_inches="tight", pad_inches=.15)
        plt.close(fig)
        print(f"EXPERIMENT_{project}_CSV_AND_PLOTS_PASS: 243 directions, 30 summaries, 21 detector checks")
        print("relative scalar signal norms:", np.array(signals)/signals[0])
    current1 = read(ROOT / "1" / "case_1_angles.csv")
    current2 = read(ROOT / "2" / "case_1_angles.csv")
    assert current1 == current2, "The common current condition must agree between projects."
    old = read(ROOT.parent / "forward" / "exp" / "oblate_window_51.csv")
    for new, previous in zip(current1, old):
        for metric in ("born_scalar", "rytov_scalar", "born_co", "rytov_co"):
            assert np.isclose(new[metric], previous["continuous_"+metric], rtol=1e-10, atol=1e-12)
    print("COMMON_BASELINE_AND_PREVIOUS_EXPERIMENT_REPRODUCED")


if __name__ == "__main__":
    main()
