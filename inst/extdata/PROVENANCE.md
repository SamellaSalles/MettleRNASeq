# Example data and processing

[Expression matrix](bone_expression.csv) · [Sample metadata](bone_metadata.csv)

The example contains **eight mouse bone samples and 235 genes**, with four
FLASH and four STANDARD replicates.

The original RNA-seq data come from
[GEO GSE173944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173944),
published by [Velalopoulou et al. (2021)](https://doi.org/10.1158/0008-5472.CAN-21-1500).
**Samella Salles processed these data during her master's research**, which
underlies MettleRNASeq. The included subset contains normalized expression
values and genes preselected by PCA.

For inclusion in the package, these processed values were reorganized into
samples in rows and genes in columns, with a separate sample-metadata file.
Expression values were preserved during this conversion. Sample names F1B–F4B
denote FLASH and S1B–S4B denote STANDARD; B denotes bone.

The [MettleRNASeq preprint](https://doi.org/10.1101/2025.05.06.652387)
describes the method and reports the scientific analyses and results.
The bundled subset provides a small example for learning the package.

These are **processed expression values**, not raw counts. Skip normalization
for this example. For your own analysis, follow the manual with your expression
matrix and sample metadata, and choose parameters for your experimental design.
Attribute the original study and Samella Salles' processing when reusing these
example data.
