# script to render vignettes for package.
# not included in package build.
library(rmarkdown)

render("package_developer_guide.Rmd",
       output_format = md_document("markdown_github"),
       output_file = "../inst/doc/package_develper_guide.md")

