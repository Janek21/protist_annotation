import json, glob, os, sys, argparse, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def run(input_pattern, out=None):
    files = sorted(glob.glob(input_pattern))

    if not files:
        print(f"No files found matching: {input_pattern}")
        return

    folder = os.path.dirname(files[0])
    output = out if out else f"{folder}/busco_summary.png"

    results = {}
    for f in files:
        basename = os.path.basename(f)

        parts = basename.replace(".json", "").split("_")
        taxid_idx = next((i for i, p in enumerate(parts) if p.isdigit()), None)
        species = " ".join(parts[:taxid_idx]) if taxid_idx else basename
        assay = parts[taxid_idx + 1] if taxid_idx else ""
        label = f"{species} ({assay})"

        with open(f) as fh:
            data = json.load(fh)

        r = data["results"]
        results[label] = {
            "C_S": round(r["Complete percentage"] - r["Multi copy percentage"], 2),
            "C_D": round(r["Multi copy percentage"], 2),
            "F":   round(r["Fragmented percentage"], 2),
            "M":   round(r["Missing percentage"], 2),
            "total": round(r.get("n_markers", 0), 2),
        }

    if not results:
        print(f"No valid BUSCO JSON files found matching: {input_pattern}")
        return

    fig, ax = plt.subplots(figsize=(12, 0.5 * len(results) + 2))
    colors = ["#56b4e9", "#0072b2", "#f0e442", "#f04442"]
    legend_labels = [
        "Complete (C) and single-copy (S)",
        "Complete (C) and duplicated (D)",
        "Fragmented (F)",
        "Missing (M)"
    ]

    for i, (sp, v) in enumerate(results.items()):
        vals = [v["C_S"], v["C_D"], v["F"], v["M"]]
        left = 0
        for val, col in zip(vals, colors):
            ax.barh(i, val, left=left, color=col, height=0.6)
            left += val

        score_text = (
            f"C:{round(v['C_S'] + v['C_D'], 2)}[S:{v['C_S']},D:{v['C_D']}],"
            f"F:{v['F']},M:{v['M']},n:{v['total']}"
        )
        ax.text(0.5, i, score_text, va="center", ha="left", fontsize=7,
                color="black", fontweight="bold")

    for x in range(0, 101, 20):
        ax.axvline(x, color="gray", linewidth=0.8, linestyle="--", zorder=3)

    ax.set_yticks(range(len(results)))
    ax.set_yticklabels(list(results.keys()), fontsize=8)
    ax.set_xlabel("% BUSCOs")
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.5, len(results) - 0.5)
    ax.set_title("BUSCO Assessment Results", fontsize=12, pad=10)

    ax.legend(
        [mpatches.Patch(color=c) for c in colors], legend_labels,
        fontsize=9, bbox_to_anchor=(0.5, -0.02), loc="upper center",
        ncol=2, borderaxespad=0
    )

    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {output} with {len(results)} entries")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a combined BUSCO summary plot.")
    parser.add_argument("-i", "--input", required=True, help="Glob pattern to BUSCO JSON files (e.g. /path/to/folder/*.json)")
    parser.add_argument("-o", "--out", default=None, help="Output PNG file (default: busco_summary.png in the input folder)")
    args = parser.parse_args()
    run(args.input, args.out)
