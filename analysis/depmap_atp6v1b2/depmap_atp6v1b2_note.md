# ATP6V1B2 DepMap stratified validation note

## Data source

- DepMap Public 25Q3, release date: 2025-09-25
- Files used:
  - `Model.csv`
  - `CRISPRGeneEffect.csv`
- Local cached paths:
  - `data/Model_25Q3.csv`
  - `data/CRISPRGeneEffect_25Q3.csv`

## Key outputs

- Grouped figure: `figures/atp6v1b2_depmap_grouped.png`
- Group summary: `tables/atp6v1b2_group_summary.csv`
- Group comparison statistics: `tables/atp6v1b2_group_comparisons.csv`
- Top myeloid models: `tables/atp6v1b2_top_myeloid_dependencies.csv`

## Main result

ATP6V1B2 shows strong dependency across multiple lineages in DepMap rather than a sharply myeloid-restricted pattern. In the 25Q3 aggregated `CRISPRGeneEffect` matrix, the mean ATP6V1B2 dependency score was `-2.375` in myeloid models (`n=43`), `-2.297` in lymphoid models (`n=93`), and `-2.157` in non-immune models (`n=1050`). All three groups were far below the conventional dependency threshold of `-0.5`, indicating that ATP6V1B2 behaves more like a broadly required factor than a lineage-exclusive vulnerability. The myeloid group trended toward stronger dependency than the non-immune group, but this shift did not reach conventional significance in the present comparison (`Mann-Whitney U p=0.073`).

## Important caveat

`HL-60 (ACH-000002)` is not present in the aggregated `CRISPRGeneEffect` matrix of the releases checked locally (`25Q3` and `24Q4`). Therefore, this DepMap analysis can support a lineage-level comparison, but it cannot be used to claim direct concordance between the DepMap matrix and the HL-60 experiment.

## Conservative paragraph for Chapter 4

为对ATP6V1B2的条件性依赖特征进行外部公共数据复核，本文进一步分析了DepMap Public 25Q3发布的`CRISPRGeneEffect.csv`与`Model.csv`。按细胞谱系将模型划分为髓系、淋系与非免疫来源后可见，ATP6V1B2在三组细胞中均表现出较强依赖，其基因效应分数中位数分别为-2.397、-2.391和-2.268，均明显低于常用依赖阈值-0.5。这说明ATP6V1B2并非仅在某一单一谱系中发挥作用，而更接近一个广泛参与细胞存活或稳态维持的关键节点。与此同时，髓系模型的平均依赖程度略强于非免疫模型（-2.375 vs -2.157），提示其在髓系背景下可能存在进一步增强的功能需求，但该差异在当前比较中尚未达到显著性水平（Mann-Whitney U检验，p=0.073）。因此，DepMap结果更适合被理解为对ATP6V1B2“广泛必需并在髓系中可能进一步增强”的补充支持，而非直接证明其具有严格的髓系专属性。

## Slightly stronger but still defensible paragraph

DepMap公共CRISPR依赖数据进一步表明，ATP6V1B2在多种肿瘤细胞背景中均维持较强依赖，但髓系模型整体分布仍较非免疫模型更偏向负值，提示其在髓系背景下可能存在额外增强的功能需求。结合本文在HL-60替代体系中观察到的PMA刺激后敲低效应进一步加重的现象，可以推测ATP6V1B2所体现的并非“从无到有”的谱系特异依赖，而更可能是一个在广泛基础必需性之上、于高应激免疫效应背景中进一步被放大的酸化动力学限制节点。这一结果与本文提出的“情境依赖性”概念并不矛盾，而是提示ATP6V1B2更接近应激放大型依赖蛋白，而非绝对谱系专有蛋白。
