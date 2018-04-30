.use.local <- FALSE; knitr::knit("de_analysis.Rmd"); rmarkdown::render("de_analysis.md")
source("make_pics.R")
source("make_pca.R")
knitr::knit("snr.Rmd"); rmarkdown::render("snr.md")
