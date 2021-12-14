# Code to pre-compile vignettes locally

all_files <- list.files("vignettes", full.names = T)
vignette_files <- grep(pattern = "*.orig", x = all_files, value = T)
rmd_names <- gsub(pattern = ".orig", replacement = "", x = vignette_files)


for (i in vignette_files) {
  knitr::knit(i, output = gsub(pattern = ".orig", replacement = "", x = i))
}
