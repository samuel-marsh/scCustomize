# Precompiled vignettes that depend on SeuratData Package or other data

library(knitr)

knitr::knit("vignettes/Read_and_Write_Functions.Rmd.orig", output = "vignettes/Read_and_Write_Functions.Rmd")
knitr::knit("vignettes/Sequencing_QC_Plots.Rmd.orig", output = "vignettes/Sequencing_QC_Plots.Rmd")
