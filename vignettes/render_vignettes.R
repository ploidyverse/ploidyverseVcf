# script to render vignettes for package.
# not included in package build.
library(rmarkdown)

render("package_developer_guide.Rmd",
       output_format = github_document(),
       output_file = "package_developer_guide.md")
file.rename("package_developer_guide.md", 
            "../inst/doc/package_developer_guide.md")
