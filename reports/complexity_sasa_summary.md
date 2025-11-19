# PTM Crosstalk Complexity vs SASA

Definitions:
- n_sites: distinct PTM sites within 10 Å of the glutathionylated cysteine.
- n_types: distinct PTM types within 10 Å.

Correlations (glut cysteines with FreeSASA):
- Spearman(SASA, n_sites) = -0.176
- Spearman(SASA, n_types) = -0.148
- Pearson(SASA, n_sites) = -0.194
- Pearson(SASA, n_types) = -0.178

Mean SASA by n_types:
- 0: 28.30 Å² (n=301)
- 1: 21.05 Å² (n=260)
- 2: 19.69 Å² (n=104)
- 3+: 11.90 Å² (n=186)

Mean SASA by n_sites:
- 0: 28.30 Å² (n=301)
- 1: 23.12 Å² (n=157)
- 2: 21.33 Å² (n=96)
- 3-4: 16.05 Å² (n=137)
- 5+: 11.60 Å² (n=160)
