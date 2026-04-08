# ATP6V1B2 DepMap thesis note

- Data source: DepMap Public 25Q3 cached local files.
- Myeloid median -2.397 vs non-immune median -2.268.
- Fraction with Chronos <= -1.0: myeloid 100.0% vs non-immune 94.4%.
- Top 5% extreme-dependency tail: myeloid 14.0% vs background 4.7% (OR 3.27, p=0.018).
- AML/MPN-like top 5% tail: 12.5% vs background 4.8% (OR 2.83, p=0.047).

## Top myeloid models

| CellLineName | OncotreePrimaryDisease | OncotreeSubtype | DiseaseGroup | ATP6V1B2 |
| --- | --- | --- | --- | --- |
| JURL-MK1 | Myeloproliferative Neoplasms | Chronic Myeloid Leukemia, BCR-ABL1+ | MPN/CML | -3.233030634826214 |
| EOL-1 | Acute Myeloid Leukemia | Acute Myeloid Leukemia | AML | -3.1531445824794577 |
| TF-1 | Acute Myeloid Leukemia | AML, NOS | AML | -3.126008625022296 |
| ARH-77 | Non-Cancerous | Immortalized Blood | Other myeloid | -3.1258424191524874 |
| HEL | Acute Myeloid Leukemia | AML, NOS | AML | -3.1088717320116306 |
| OCI-M2 | Acute Myeloid Leukemia | Acute Myeloid Leukemia | AML | -3.098213923644433 |
| KU812 | Myeloproliferative Neoplasms | Chronic Myeloid Leukemia, BCR-ABL1+ | MPN/CML | -2.9608996595572243 |
| AML-193 | Acute Myeloid Leukemia | AML with Myelodysplasia-Related Changes | AML | -2.885013551816702 |
| KYO-1 | Myeloproliferative Neoplasms | Chronic Myeloid Leukemia, BCR-ABL1+ | MPN/CML | -2.8309426275356597 |
| HEL 92.1.7 | Acute Myeloid Leukemia | Acute Myeloid Leukemia | AML | -2.821905867551924 |
| MUTZ-8 | Acute Myeloid Leukemia | Acute Myeloid Leukemia | AML | -2.8148950385260365 |
| HD-MY-Z | Acute Myeloid Leukemia | Acute Myeloid Leukemia | AML | -2.747559119451295 |

## Mann-Whitney comparisons

| comparison | n_left | n_right | mean_left | mean_right | median_left | median_right | mannwhitney_u | pvalue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Myeloid_vs_rest | 43 | 1143 | -2.375369998979494 | -2.168195586405128 | -2.396677534225488 | -2.2785642582105976 | 20851.0 | 0.0913153696463338 |
| AML_MPN_like_vs_rest | 40 | 1146 | -2.3796813290217544 | -2.1685874450753135 | -2.4014741172226906 | -2.277813329916249 | 19260.0 | 0.08569451225243367 |
| Lymphoid_vs_rest | 93 | 1093 | -2.297363359783653 | -2.165355601790759 | -2.3910230528078813 | -2.2724525360469414 | 45299.0 | 0.0814397890499071 |