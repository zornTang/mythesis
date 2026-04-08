from __future__ import annotations

from pathlib import Path
from string import ascii_uppercase

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu


TARGET_GENE = "ATP6V1B2"
MODULE_GENE = "ATP6V1A"

GROUP_COLORS = {
    "Myeloid": "#C44E52",
    "Lymphoid": "#4C72B0",
    "Non-immune": "#7A7A7A",
}
SUBTYPE_COLORS = {
    "AML": "#C44E52",
    "MPN/CML": "#DD8452",
    "Other myeloid": "#55A868",
}
TAIL_COLORS = {
    "Top 5% tail": "#4C72B0",
    "Top 15% tail": "#55A868",
}


def load_data(root: Path) -> pd.DataFrame:
    data_dir = root / "analysis" / "depmap_atp6v1b2" / "data"
    model = pd.read_csv(data_dir / "Model_25Q3.csv")

    effect = pd.read_csv(
        data_dir / "CRISPRGeneEffect_25Q3.csv",
        usecols=["Unnamed: 0", "ATP6V1B2 (526)", "ATP6V1A (523)"],
    ).rename(
        columns={
            "Unnamed: 0": "ModelID",
            "ATP6V1B2 (526)": TARGET_GENE,
            "ATP6V1A (523)": MODULE_GENE,
        }
    )

    keep = [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
    ]
    merged = effect.merge(model[keep], on="ModelID", how="left")
    merged["DependencyGroup"] = merged["OncotreeLineage"].map(assign_group)
    merged["DiseaseGroup"] = merged.apply(assign_myeloid_subtype, axis=1)
    return merged


def assign_group(lineage: str) -> str:
    if lineage == "Myeloid":
        return "Myeloid"
    if lineage == "Lymphoid":
        return "Lymphoid"
    return "Non-immune"


def assign_myeloid_subtype(row: pd.Series) -> str:
    disease = str(row.get("OncotreePrimaryDisease", "") or "")
    subtype = str(row.get("OncotreeSubtype", "") or "")
    text = f"{disease} | {subtype}"
    if "Acute Myeloid Leukemia" in text:
        return "AML"
    if "Myeloproliferative" in text or "Chronic Myeloid Leukemia" in text:
        return "MPN/CML"
    return "Other myeloid"


def summarize_groups(merged: pd.DataFrame) -> pd.DataFrame:
    order = ["Myeloid", "Lymphoid", "Non-immune"]
    rows = []
    for group in order:
        values = merged.loc[merged["DependencyGroup"] == group, TARGET_GENE].dropna()
        rows.append(
            {
                "group": group,
                "n": len(values),
                "mean": values.mean(),
                "median": values.median(),
                "q1": values.quantile(0.25),
                "q3": values.quantile(0.75),
                "fraction_le_neg_0_5": (values <= -0.5).mean(),
                "fraction_le_neg_1_0": (values <= -1.0).mean(),
                "fraction_le_neg_1_5": (values <= -1.5).mean(),
                "fraction_le_neg_2_0": (values <= -2.0).mean(),
            }
        )
    return pd.DataFrame(rows)


def summarize_tail_enrichment(merged: pd.DataFrame) -> pd.DataFrame:
    groups = {
        "Myeloid": merged["DependencyGroup"] == "Myeloid",
        "AML/MPN-like": (
            (merged["OncotreePrimaryDisease"].fillna("").str.contains("Acute Myeloid Leukemia", case=False))
            | (merged["OncotreePrimaryDisease"].fillna("").str.contains("Myeloproliferative", case=False))
        ),
        "Lymphoid": merged["DependencyGroup"] == "Lymphoid",
    }
    rows = []
    for q, label in [(0.05, "Top 5% tail"), (0.15, "Top 15% tail")]:
        cutoff = merged[TARGET_GENE].quantile(q)
        extreme = merged[TARGET_GENE] <= cutoff
        for group_name, mask in groups.items():
            in_tail = int((mask & extreme).sum())
            not_in_tail = int((mask & ~extreme).sum())
            bg_tail = int((~mask & extreme).sum())
            bg_not_tail = int((~mask & ~extreme).sum())
            odds_ratio, pvalue = fisher_exact(
                [[in_tail, not_in_tail], [bg_tail, bg_not_tail]],
                alternative="greater",
            )
            rows.append(
                {
                    "tail_label": label,
                    "quantile": q,
                    "cutoff": cutoff,
                    "group": group_name,
                    "in_tail": in_tail,
                    "group_total": int(mask.sum()),
                    "tail_fraction": in_tail / max(int(mask.sum()), 1),
                    "background_fraction": bg_tail / max(int((~mask).sum()), 1),
                    "odds_ratio": odds_ratio,
                    "pvalue": pvalue,
                }
            )
    return pd.DataFrame(rows)


def summarize_top_myeloid(merged: pd.DataFrame, n: int = 12) -> pd.DataFrame:
    myeloid = merged.loc[merged["DependencyGroup"] == "Myeloid"].copy()
    return myeloid.sort_values(TARGET_GENE).head(n)[
        ["CellLineName", "OncotreePrimaryDisease", "OncotreeSubtype", "DiseaseGroup", TARGET_GENE]
    ]


def compare_with_rest(merged: pd.DataFrame) -> pd.DataFrame:
    masks = {
        "Myeloid_vs_rest": merged["DependencyGroup"] == "Myeloid",
        "AML_MPN_like_vs_rest": (
            (merged["OncotreePrimaryDisease"].fillna("").str.contains("Acute Myeloid Leukemia", case=False))
            | (merged["OncotreePrimaryDisease"].fillna("").str.contains("Myeloproliferative", case=False))
        ),
        "Lymphoid_vs_rest": merged["DependencyGroup"] == "Lymphoid",
    }
    rows = []
    for label, mask in masks.items():
        left = merged.loc[mask, TARGET_GENE].dropna()
        right = merged.loc[~mask, TARGET_GENE].dropna()
        stat, pvalue = mannwhitneyu(left, right, alternative="two-sided")
        rows.append(
            {
                "comparison": label,
                "n_left": len(left),
                "n_right": len(right),
                "mean_left": left.mean(),
                "mean_right": right.mean(),
                "median_left": left.median(),
                "median_right": right.median(),
                "mannwhitney_u": stat,
                "pvalue": pvalue,
            }
        )
    return pd.DataFrame(rows)


def plot_distribution(ax: plt.Axes, merged: pd.DataFrame) -> None:
    order = ["Myeloid", "Lymphoid", "Non-immune"]
    plot_df = merged.dropna(subset=[TARGET_GENE]).copy()
    sns.violinplot(
        data=plot_df,
        x="DependencyGroup",
        y=TARGET_GENE,
        order=order,
        inner=None,
        linewidth=0,
        cut=0,
        palette=[GROUP_COLORS[g] for g in order],
        ax=ax,
    )
    sns.boxplot(
        data=plot_df,
        x="DependencyGroup",
        y=TARGET_GENE,
        order=order,
        width=0.22,
        showcaps=True,
        showfliers=False,
        boxprops={"facecolor": "white", "zorder": 3},
        whiskerprops={"linewidth": 1.2},
        medianprops={"color": "black", "linewidth": 1.4},
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="DependencyGroup",
        y=TARGET_GENE,
        order=order,
        color="black",
        alpha=0.28,
        size=2.5,
        jitter=0.18,
        ax=ax,
    )
    ax.axhline(-1.0, color="#8C1D18", linestyle="--", linewidth=1.0)
    ax.axhline(-0.5, color="#BEBEBE", linestyle=":", linewidth=1.0)
    ax.text(2.35, -0.98, "strong dependency = -1.0", fontsize=8, color="#8C1D18", ha="right", va="bottom")
    ax.text(2.35, -0.48, "dependency = -0.5", fontsize=8, color="#707070", ha="right", va="bottom")
    summary = summarize_groups(merged)
    for idx, group in enumerate(order):
        row = summary.loc[summary["group"] == group].iloc[0]
        ax.text(idx, row["q3"] + 0.13, f"n={int(row['n'])}\nmed={row['median']:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("ATP6V1B2 Chronos gene effect")
    ax.set_title("Lineage-level distribution")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_ecdf(ax: plt.Axes, merged: pd.DataFrame) -> None:
    order = ["Myeloid", "Lymphoid", "Non-immune"]
    for group in order:
        values = np.sort(merged.loc[merged["DependencyGroup"] == group, TARGET_GENE].dropna().to_numpy())
        y = np.arange(1, len(values) + 1) / len(values)
        ax.step(values, y, where="post", color=GROUP_COLORS[group], linewidth=2.0, label=group)
    ax.axvline(-1.0, color="#8C1D18", linestyle="--", linewidth=1.0)
    ax.axvline(-1.5, color="#BEBEBE", linestyle=":", linewidth=1.0)
    ax.set_xlabel("ATP6V1B2 Chronos gene effect")
    ax.set_ylabel("Cumulative fraction of cell lines")
    ax.set_title("Uniformity of strong dependency")
    ax.legend(frameon=False, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_tail_enrichment(ax: plt.Axes, tail_df: pd.DataFrame) -> None:
    groups = ["Myeloid", "AML/MPN-like", "Lymphoid"]
    x = np.arange(len(groups))
    width = 0.33

    for i, tail_label in enumerate(["Top 5% tail", "Top 15% tail"]):
        subset = tail_df[tail_df["tail_label"] == tail_label].set_index("group").loc[groups]
        positions = x + (i - 0.5) * width
        ax.bar(
            positions,
            subset["tail_fraction"],
            width=width,
            color=TAIL_COLORS[tail_label],
            edgecolor="white",
            linewidth=0.8,
            label=tail_label,
        )
        expected = 0.05 if tail_label == "Top 5% tail" else 0.15
        ax.axhline(expected, color=TAIL_COLORS[tail_label], linestyle=":", linewidth=1.0, alpha=0.7)
        for xpos, (_, row) in zip(positions, subset.iterrows()):
            p_label = f"p={row['pvalue']:.3f}" if row["pvalue"] >= 0.001 else "p<0.001"
            ax.text(
                xpos,
                row["tail_fraction"] + 0.018,
                f"{row['tail_fraction']*100:.1f}%\nOR={row['odds_ratio']:.2f}\n{p_label}",
                ha="center",
                va="bottom",
                fontsize=7.2,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylim(0, 0.38)
    ax.set_ylabel("Fraction of models in extreme-dependency tail")
    ax.set_title("Tail enrichment among hematopoietic contexts")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_top_myeloid(ax: plt.Axes, top_df: pd.DataFrame) -> None:
    plot_df = top_df.sort_values(TARGET_GENE, ascending=True).copy()
    colors = plot_df["DiseaseGroup"].map(SUBTYPE_COLORS)
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df[TARGET_GENE], color=colors, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["CellLineName"])
    ax.invert_yaxis()
    ax.set_xlabel("ATP6V1B2 Chronos gene effect")
    ax.set_title("Most dependent myeloid models")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = []
    labels = []
    for label, color in SUBTYPE_COLORS.items():
        handles.append(plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=color, markersize=8))
        labels.append(label)
    ax.legend(handles, labels, frameon=False, loc="lower right", fontsize=8)


def add_panel_labels(axes: list[plt.Axes]) -> None:
    for idx, ax in enumerate(axes):
        ax.text(-0.15, 1.06, ascii_uppercase[idx], transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def build_main_figure(merged: pd.DataFrame, tail_df: pd.DataFrame, top_df: pd.DataFrame, out_base: Path) -> None:
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )
    fig = plt.figure(figsize=(7.3, 6.2))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.34)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    plot_distribution(ax1, merged)
    plot_ecdf(ax2, merged)
    plot_tail_enrichment(ax3, tail_df)
    plot_top_myeloid(ax4, top_df)
    add_panel_labels([ax1, ax2, ax3, ax4])

    fig.suptitle("DepMap supports a stronger myeloid ATP6V1B2 dependency tail", y=0.995, fontsize=11.5)
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def build_module_figure(merged: pd.DataFrame, out_base: Path) -> None:
    plot_df = merged.melt(
        id_vars=["DependencyGroup"],
        value_vars=[TARGET_GENE, MODULE_GENE],
        var_name="Gene",
        value_name="Chronos",
    )
    order = ["Myeloid", "Lymphoid", "Non-immune"]

    sns.set_theme(style="whitegrid", context="paper")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=False)
    for ax, gene in zip(axes, [TARGET_GENE, MODULE_GENE]):
        subset = plot_df[plot_df["Gene"] == gene]
        sns.boxplot(
            data=subset,
            x="DependencyGroup",
            y="Chronos",
            order=order,
            width=0.48,
            showfliers=False,
            palette=[GROUP_COLORS[g] for g in order],
            ax=ax,
        )
        sns.stripplot(
            data=subset,
            x="DependencyGroup",
            y="Chronos",
            order=order,
            color="black",
            alpha=0.2,
            size=2,
            jitter=0.2,
            ax=ax,
        )
        ax.set_title(gene)
        ax.set_xlabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("Chronos gene effect")
    axes[1].set_ylabel("")
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def build_summary_note(
    root: Path,
    group_summary: pd.DataFrame,
    tail_summary: pd.DataFrame,
    comparisons: pd.DataFrame,
    top_myeloid: pd.DataFrame,
) -> None:
    def df_to_markdown(df: pd.DataFrame) -> str:
        header = "| " + " | ".join(df.columns.astype(str)) + " |"
        sep = "| " + " | ".join(["---"] * len(df.columns)) + " |"
        rows = []
        for _, row in df.iterrows():
            rows.append("| " + " | ".join(str(v) for v in row.tolist()) + " |")
        return "\n".join([header, sep] + rows)

    out_path = root / "analysis" / "depmap_atp6v1b2" / "depmap_atp6v1b2_thesis_note.md"
    myeloid = group_summary.set_index("group").loc["Myeloid"]
    nonimmune = group_summary.set_index("group").loc["Non-immune"]
    top5_myeloid = tail_summary.query("group == 'Myeloid' and tail_label == 'Top 5% tail'").iloc[0]
    top5_aml = tail_summary.query("group == 'AML/MPN-like' and tail_label == 'Top 5% tail'").iloc[0]
    lines = [
        "# ATP6V1B2 DepMap thesis note",
        "",
        "- Data source: DepMap Public 25Q3 cached local files.",
        f"- Myeloid median {myeloid['median']:.3f} vs non-immune median {nonimmune['median']:.3f}.",
        f"- Fraction with Chronos <= -1.0: myeloid {myeloid['fraction_le_neg_1_0']*100:.1f}% vs non-immune {nonimmune['fraction_le_neg_1_0']*100:.1f}%.",
        f"- Top 5% extreme-dependency tail: myeloid {top5_myeloid['tail_fraction']*100:.1f}% vs background {top5_myeloid['background_fraction']*100:.1f}% (OR {top5_myeloid['odds_ratio']:.2f}, p={top5_myeloid['pvalue']:.3f}).",
        f"- AML/MPN-like top 5% tail: {top5_aml['tail_fraction']*100:.1f}% vs background {top5_aml['background_fraction']*100:.1f}% (OR {top5_aml['odds_ratio']:.2f}, p={top5_aml['pvalue']:.3f}).",
        "",
        "## Top myeloid models",
        "",
        df_to_markdown(top_myeloid),
        "",
        "## Mann-Whitney comparisons",
        "",
        df_to_markdown(comparisons),
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "analysis" / "depmap_atp6v1b2"
    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    thesis_figures_dir = root / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    thesis_figures_dir.mkdir(parents=True, exist_ok=True)

    merged = load_data(root)
    group_summary = summarize_groups(merged)
    tail_summary = summarize_tail_enrichment(merged)
    top_myeloid = summarize_top_myeloid(merged)
    comparisons = compare_with_rest(merged)

    group_summary.to_csv(tables_dir / "atp6v1b2_group_summary_refined.csv", index=False)
    tail_summary.to_csv(tables_dir / "atp6v1b2_tail_enrichment.csv", index=False)
    top_myeloid.to_csv(tables_dir / "atp6v1b2_top_myeloid_models_refined.csv", index=False)
    comparisons.to_csv(tables_dir / "atp6v1b2_group_comparisons_refined.csv", index=False)

    main_base = figures_dir / "atp6v1b2_depmap_thesis_main"
    module_base = figures_dir / "vatpase_module_depmap_lineages"
    build_main_figure(merged, tail_summary, top_myeloid, main_base)
    build_module_figure(merged, module_base)

    for stem in ["atp6v1b2_depmap_thesis_main", "vatpase_module_depmap_lineages"]:
        for suffix in [".png", ".pdf"]:
            src = figures_dir / f"{stem}{suffix}"
            dst = thesis_figures_dir / f"{stem}{suffix}"
            dst.write_bytes(src.read_bytes())

    build_summary_note(root, group_summary, tail_summary, comparisons, top_myeloid)

    print(group_summary.to_string(index=False))
    print()
    print(tail_summary.to_string(index=False))
    print()
    print(top_myeloid.to_string(index=False))


if __name__ == "__main__":
    main()
