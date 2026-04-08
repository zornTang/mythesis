from __future__ import annotations

from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import mannwhitneyu


MODEL_URL = "https://ndownloader.figshare.com/files/51065297"
GENE_EFFECT_URL = "https://ndownloader.figshare.com/files/51064667"
TARGET_GENE = "ATP6V1B2"


def load_model() -> pd.DataFrame:
    model = pd.read_csv(MODEL_URL)
    keep = [
        "ModelID",
        "CellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
    ]
    return model[keep].copy()


def resolve_gene_column() -> str:
    header = pd.read_csv(GENE_EFFECT_URL, nrows=0)
    for column in header.columns:
        if column.startswith(f"{TARGET_GENE} "):
            return column
    raise KeyError(f"Could not find a CRISPRGeneEffect column for {TARGET_GENE}.")


def load_gene_effect(gene_column: str) -> pd.DataFrame:
    header = pd.read_csv(GENE_EFFECT_URL, nrows=0)
    first_col = header.columns[0]
    effect = pd.read_csv(GENE_EFFECT_URL, usecols=[first_col, gene_column])
    effect = effect.rename(columns={first_col: "ModelID", gene_column: TARGET_GENE})
    return effect


def assign_group(lineage: str) -> str:
    if lineage == "Myeloid":
        return "Myeloid"
    if lineage == "Lymphoid":
        return "Lymphoid"
    return "Non-immune"


def summarize(merged: pd.DataFrame) -> pd.DataFrame:
    summary = (
        merged.groupby("DependencyGroup")[TARGET_GENE]
        .agg(
            n="count",
            mean="mean",
            median="median",
            min="min",
            max="max",
            dependent_fraction=lambda s: (s <= -0.5).mean(),
        )
        .reset_index()
    )
    return summary


def compare_groups(merged: pd.DataFrame) -> pd.DataFrame:
    comparisons = [
        ("Myeloid", "Non-immune"),
        ("Lymphoid", "Non-immune"),
        ("Myeloid", "Lymphoid"),
    ]
    rows = []
    for left, right in comparisons:
        left_values = merged.loc[merged["DependencyGroup"] == left, TARGET_GENE].dropna()
        right_values = merged.loc[merged["DependencyGroup"] == right, TARGET_GENE].dropna()
        stat, pvalue = mannwhitneyu(left_values, right_values, alternative="two-sided")
        rows.append(
            {
                "group_a": left,
                "group_b": right,
                "n_a": len(left_values),
                "n_b": len(right_values),
                "mean_a": left_values.mean(),
                "mean_b": right_values.mean(),
                "median_a": left_values.median(),
                "median_b": right_values.median(),
                "mannwhitney_u": stat,
                "pvalue": pvalue,
            }
        )
    return pd.DataFrame(rows)


def build_plot(merged: pd.DataFrame, output_path: Path) -> None:
    order = ["Myeloid", "Lymphoid", "Non-immune"]
    plot_df = merged.dropna(subset=[TARGET_GENE]).copy()

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    sns.boxplot(
        data=plot_df,
        x="DependencyGroup",
        y=TARGET_GENE,
        hue="DependencyGroup",
        order=order,
        width=0.55,
        showfliers=False,
        palette=["#B33A3A", "#2C7FB8", "#7F8C8D"],
        legend=False,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="DependencyGroup",
        y=TARGET_GENE,
        order=order,
        color="black",
        alpha=0.35,
        size=3,
        jitter=0.22,
        ax=ax,
    )

    hl60 = plot_df[plot_df["CellLineName"] == "HL-60"]
    if not hl60.empty:
        hl60_y = float(hl60.iloc[0][TARGET_GENE])
        ax.scatter([0], [hl60_y], color="#F28E2B", s=85, zorder=5)
        ax.annotate(
            f"HL-60 ({hl60_y:.2f})",
            (0, hl60_y),
            xytext=(12, -8),
            textcoords="offset points",
            fontsize=10,
            color="#8A4F00",
        )

    ax.axhline(-0.5, linestyle="--", linewidth=1.2, color="#8C1D18")
    ax.text(
        2.35,
        -0.48,
        "dependency threshold = -0.5",
        fontsize=9,
        color="#8C1D18",
        va="bottom",
        ha="right",
    )
    ax.set_xlabel("")
    ax.set_ylabel(f"{TARGET_GENE} CRISPR gene effect")
    ax.set_title(f"DepMap stratified dependency for {TARGET_GENE}")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "analysis" / "depmap_atp6v1b2"
    figures_dir = out_dir / "figures"
    tables_dir = out_dir / "tables"
    data_dir = out_dir / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    global MODEL_URL, GENE_EFFECT_URL
    local_model = data_dir / "Model_25Q3.csv"
    local_effect = data_dir / "CRISPRGeneEffect_25Q3.csv"
    if not local_model.exists():
        local_model = data_dir / "Model_24Q4.csv"
    if not local_effect.exists():
        local_effect = data_dir / "CRISPRGeneEffect_24Q4.csv"
    if local_model.exists():
        MODEL_URL = str(local_model)
    if local_effect.exists():
        GENE_EFFECT_URL = str(local_effect)

    model = load_model()
    gene_column = resolve_gene_column()
    effect = load_gene_effect(gene_column)

    merged = effect.merge(model, on="ModelID", how="left")
    merged["DependencyGroup"] = merged["OncotreeLineage"].map(assign_group)

    summary = summarize(merged)
    comparisons = compare_groups(merged)
    myeloid_rank = (
        merged[merged["DependencyGroup"] == "Myeloid"]
        .sort_values(TARGET_GENE)
        [["ModelID", "CellLineName", "OncotreePrimaryDisease", TARGET_GENE]]
    )

    summary.to_csv(tables_dir / "atp6v1b2_group_summary.csv", index=False)
    comparisons.to_csv(tables_dir / "atp6v1b2_group_comparisons.csv", index=False)
    myeloid_rank.head(20).to_csv(
        tables_dir / "atp6v1b2_top_myeloid_dependencies.csv", index=False
    )
    merged.to_csv(tables_dir / "atp6v1b2_merged_scores.csv", index=False)

    build_plot(merged, figures_dir / "atp6v1b2_depmap_grouped.png")

    hl60 = merged[merged["CellLineName"] == "HL-60"]
    print(f"Resolved gene column: {gene_column}")
    print(f"Model source: {MODEL_URL}")
    print(f"Gene effect source: {GENE_EFFECT_URL}")
    print(summary.to_string(index=False))
    print()
    print("Group comparisons:")
    print(comparisons.to_string(index=False))
    print()
    if hl60.empty:
        print("HL-60 is not present in the aggregated CRISPRGeneEffect matrix for this release.")
    else:
        row = hl60.iloc[0]
        print(
            "HL-60:",
            {
                "ModelID": row["ModelID"],
                "DependencyGroup": row["DependencyGroup"],
                "OncotreePrimaryDisease": row["OncotreePrimaryDisease"],
                TARGET_GENE: round(float(row[TARGET_GENE]), 4),
            },
        )
    print()
    print("Top 10 myeloid dependencies:")
    print(myeloid_rank.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
