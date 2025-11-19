# Key Findings

## Enrichment by PTM (<=10 Å)
| PTM | Glut yes | Glut no | Ctrl yes | Ctrl no | Odds Ratio | p-value |
|---|---:|---:|---:|---:|---:|---:|
| Acetylation | 285 | 574 | 890 | 3897 | 2.18 | 2.84e-20 |
| Malonylation | 111 | 748 | 308 | 4479 | 2.16 | 4.82e-10 |
| Methylation | 25 | 834 | 65 | 4722 | 2.20 | 1.75e-03 |
| N-linked Glycosylation | 3 | 856 | 64 | 4723 | 0.30 | 9.59e-03 |
| O-linked Glycosylation | 24 | 835 | 42 | 4745 | 3.27 | 1.79e-05 |
| Phosphorylation | 466 | 393 | 1615 | 3172 | 2.33 | 1.64e-29 |
| Succinylation | 105 | 754 | 240 | 4547 | 2.64 | 1.74e-13 |
| Ubiquitination | 158 | 701 | 515 | 4272 | 1.87 | 1.50e-09 |

## Modified vs. Unmodified Cysteines (<=10 Å, Any PTM)
- Odds Ratio (glut vs control): 2.26; p-value: 4.57e-27

| Group | Crosstalk ≤10 | Total | Rate |
|---|---:|---:|---:|
| glut | 553 | 859 | 64.38% |
| control | 2128 | 4787 | 44.45% |

## SASA (FreeSASA)
| Group | N | Mean SASA (Å²) | Median SASA (Å²) |
|---|---:|---:|---:|
| Glut + Crosstalk (≤10Å) | 550 | 17.70 | 5.16 |
| Glut (No Crosstalk) | 301 | 28.30 | 11.08 |
| Control (All Cys) | 4681 | 24.26 | 7.98 |

## Average SASA: Modified vs Unmodified (All Cysteines)
| Group | N | Mean SASA (Å²) | Median SASA (Å²) |
|---|---:|---:|---:|
| Modified (Glut Cys) | 851 | 21.45 | 6.42 |
| Unmodified (Control Cys) | 4681 | 24.26 | 7.98 |

## SASA vs Crosstalk Complexity (Glut Cys, FreeSASA)
More PTM types/sites within 10 Å associate with lower SASA.

### By PTM Types within 10 Å
| n_types | N | Mean SASA (Å²) | Median SASA (Å²) |
|---:|---:|---:|---:|
| 0 | 301 | 28.30 | 11.08 |
| 1 | 260 | 21.05 | 5.83 |
| 2 | 104 | 19.69 | 5.45 |
| 3+ | 186 | 11.90 | 3.50 |

### By PTM Sites within 10 Å
| n_sites | N | Mean SASA (Å²) | Median SASA (Å²) |
|---:|---:|---:|---:|
| 0 | 301 | 28.30 | 11.08 |
| 1 | 157 | 23.12 | 6.77 |
| 2 | 96 | 21.33 | 6.14 |
| 3-4 | 137 | 16.05 | 3.57 |
| 5+ | 160 | 11.60 | 2.58 |
