# ploidyverseClasses
This R package extends the `VCF` S4 class from Bioconductor's 
[VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)
package, with new methods to facilitate import and export of VCF data using
ploidyverse software.

A draft specification for VCF within the ploidyverse can be found at
[inst/doc/ploidyverse_VCF_specification.md](https://github.com/ploidyverse/ploidyverseClasses/blob/master/inst/doc/ploidyverse_VCF_specification.md).

This package serves three purposes for R package developers in the ploidyverse:

1. Checking whether VCFs meet ploidyverse specifications.
2. Facilitating the creation of VCFs that meet ploidyverse specifications.
3. Providing easy access to genotype calls and other data in a ploidyverse VCF.

Overall, the goal of this package is to make it as painless as possible for 
ploidyverse packages to integrate with each other.  See the Package Developer 
Guide at [inst/doc/package_developer_guide.Rmd](https://github.com/ploidyverse/ploidyverseClasses/blob/master/inst/doc/package_developer_guide.md)
for more information.
