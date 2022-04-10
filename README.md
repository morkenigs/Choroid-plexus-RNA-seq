# Choroid plexus RNA-seq
Codes used for the analysis of RNA-seq data performed in the paper "Cholesterol 24-hydroxylase at the choroid plexus
restrains local inflammation and protects overall brain function".
## Recommendations
Before running any of the R scripts please install the following R libraries :
- **DESeq2** (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
- **ggplot2** (https://ggplot2.tidyverse.org)
- **pheatmap** (https://cran.r-project.org/web/packages/pheatmap/index.html)
- **FactoMineR** (https://cran.r-project.org/web/packages/FactoMineR/index.html)
- **apeglm** (https://bioconductor.org/packages/release/bioc/html/apeglm.html)
- **BiocParallel** (https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
- **Hmisc** (https://cran.r-project.org/web/packages/Hmisc/index.html)
- **stringr** (https://www.rdocumentation.org/packages/stringr/versions/1.4.0)
- **ggpubr** (https://rpkgs.datanovia.com/ggpubr/)
- **RColorBrewer** (https://cran.r-project.org/web/packages/RColorBrewer/index.html)
## List of scripts
- **CPE_cultures_24OH_analysis.R** - R script used to analyze the RNAseq data generated from choroid plexus epithelium cultures treated with 24-OH or DMSO as control.
- **CYP46A1_overexpression_analysis.R** - R script used to analyze the RNAseq data generated from choroid plexus 6 weeks post injection with AAV overexpressing CYP46A1 or GFP as control.
