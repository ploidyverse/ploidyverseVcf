Using ploidyverseVcf in your R package
================
Lindsay V. Clark, University of Illinois, Urbana-Champaign
07 April 2019

## Purpose of this document

This document is intended for R package developers who wish to utilize
S4 classes from the ploidyverse in order to make their package
interoperable with other ploidyverse packages. It assumes basic
familiarity with the process of creating an R package.

## How to indicate dependency on `ploidyverseVcf`

The first question you should ask yourself is, will the ability to read
and write VCFs be mandatory for using your package? Or, will users
frequently want to import or export data in a different format, and not
care about incorporating other ploidyverse software into their workflow?

`ploidyverseVcf` depends on the Bioconductor package
[`VariantAnnotation`](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html).
Like many Bioconductor packages, `VariantAnnotation` can be unavoidably
cumbersome to install. We on the ploidyverse core team believe that it
is worthwhile for its flexibility in importing, storing, and exporting
genomic variants and metadata, but we want to minimize how much we force
it onto people who don’t want it.

### Option A: mandatory installation of `ploidyverseVcf`

In this case, you will indicate dependency on `ploidyverseVcf` in the
standard way that you would for any R package. In your `DESCRIPTION`
file, list `ploidyverseVcf` on the `Imports` line.

    Imports: ploidyverseVcf

In the `NAMESPACE` file, you should include the lines:

  - an `importFrom` line for any functions you are using from
    `ploidyverseVcf`
  - `importFrom` lines for any functions you are using from
    `VariantAnnotation` or other packages

For any ploidyverse or Bioconductor functions that you use, I recommend
looking at the help page to determine exactly which package it
originated from.

Remember that any package with an `importFrom` line will also need a
mention in the `Imports` line of your `DESCRIPTION` file.

### Option B: optional installation of `ploidyverseVcf`

If you don’t want to force someone to install `ploidyverseVcf` (and its
dependencies) in order to use your package, you should instead list
`ploidyverseVcf` and `methods` on the `Suggests` line of your
`DESCRIPTION` file. Also list any Bioconductor packages whose functions
you want to use in `Suggests`.

    Suggests: ploidyverseVcf, methods

Then within your function definitions, use the `::` operator to indicate
the package that a function comes from. For example, to use the
`seqinfo` function, type `GenomeInfoDb::seqinfo`.

``` r
# a trivial example
GetSeqInfo <- function(x){
  GenomeInfoDb::seqinfo(x)
}
```

If you want your function to behave differently depending on whether or
not `ploidyverseVcf` is installed, you can use the function
`requireNamespace`.

``` r
if(requireNamespace("ploidyverseVcf", quietly = TRUE)){
  cat("Let's make a VCF!")
} else {
  cat("Let's make something else.")
}
```

### Using Rcpp functions from `ploidyverseVcf`

The `ploidyverseVcf` package includes a set of utility functions
implemented in `Rcpp` to assist with calling, importing, and exporting
multiallelic polyploid genotypes. See `?nGen`.

To use any of those functions from within R, you can follow option A or
B above.

To use any of those functions from within C++:

  - Add the line `LinkingTo: ploidyverseVcf, Rcpp` to your `DESCRIPTION`
    file.
  - Add the lines `#include <ploidyverseVcf.h>` and `using namespace
    ploidyverseVCF;` to your C++ file(s) that use functions from
    `ploidyverseVcf`.

To be able to test your C++ files using `sourceCpp`, be sure to have
them in the `src` directory of a package with the `DESCRIPTION` file
modified as described above. If you put the C++ file in a different
directory, `sourceCpp` will not be able to find `ploidyverseVcf.h`.

## Importing or creating a VCF

VCF files can be imported using `VariantAnnotation`’s `readVcf`
function. The tutorial for your package might direct users to look at
the `readVcf` help page, or your package might include a custom import
function that uses `readVcf` internally. Note that the `param` argument
is very helpful for adjusting which samples, genomic regions, and fields
(GT, AD, etc.) get imported, saving valuable memory.

If you want to specify genomic regions to import, or if you want to loop
through a VCF in manageable chunks, you should compress the file with
`bgzip`, index it with `indexTabix`, and create a connection with
`TabixFile`.

Within R, the data from the VCF will be in an S4 object of class
`collapsedVCF` or `expandedVCF`. These two classes differ in terms of
how they handle sites with multiple alternative alleles. Both are
subclasses of the `VCF` virual class. An object of either of these
classes can be created with the `VCF` constructor function. If your
software will import data in a format other than VCF, then export to
VCF, you will need the constructor function.

## Helpful accessors from `VariantAnnotation`

Given a `VCF` object called `myvcf`:

  - `geno(myvcf)$AD` returns a two dimensional list (marker x sample) of
    integer vectors of allelic read depth. The first allele in each
    vector is the reference allele.
  - `geno(myvcf)$GP` returns genotype posterior probabilities in a
    similar format, in numeric vectors.
  - `samples(header(myvcf))` returns a character vector of sample names,
    taken from the column headers in the VCF file.
  - `rowRanges(myvcf)` returns the location of each marker, and other
    information such as reference and alternative alleles, in `GRanges`
    format.

## Adding genotype calls to a `VCF` object

The function `array3D_to_matrixList` can convert a 3D array, such as
those stored in the `genotypeLikelihood` and `posteriorProb` slots of a
`RADdata` object in `polyRAD`, into the matrix-list format required in
`VariantAnnotation`. So, if `myvcf` is a `VCF` object and `myrad` is a
`RADdata` object, one could do something like:

``` r
geno(myvcf)$GP <- 
  array3D_to_matrixList(myrad$posteriorProb[[1]][-OneAllelePerMarker(myrad)])
```

A matrix of posterior mean genotypes could to go into the `GN` slot.

*Some feedback from other package developers might be helpful here to
make it easy for them to get their genotype calls into the VCF.*

## Adding metadata suggested in the ploidyverse specifications

Note: `ploidyverseVcf` will include custom accessor functions for all of
these.

  - `sampleinfo` is used for adding metadata about samples.
  - `software` is used for adding metadata about the software used for
    genotype calling.

## Checking whether a `VCF` object meets ploidyverse specifications

To check whether a `VCF` object meets ploidyverse specifications, run

``` r
myvcf <- markValidity(myvcf)
```

There are three levels of validity that are checked:

  - “Precall” validity indicates that allelic read depths, in the AD
    field, are present, and hence it is possible to call genotypes.
  - “Postcall” validity indicates that genotype posterior probabilities,
    in the GP field, are present, and hence genotype calling has been
    performed and calls can be used in downstream analysis.
  - “Archival” validity indicates that information such as sample
    species and ploidy, reference genome version, and tag sequences in
    the absence of a reference, are present. Files meeting these
    standards are suitable for depositing at services such as Figshare
    and Zenodo.

If you want a `TRUE` or `FALSE` value to tell you whether the `VCF`
object meets these specifications (for example, for error checking
before running a function), after running `markValidity` you can run:

``` r
validPloidyverseVCF_Precall(myvcf)
validPloidyverseVCF_Postcall(myvcf)
validPloidyverseVCF_Archival(myvcf)
```

## Writing a VCF file

Before exporting a VCF, the `markValidity` function should be run on the
`VCF` object, as it will add header lines confirming whether or not the
file meets ploidyverse specifications. Any metadata that was added to
the object will also be exported to the file.

The object is then passed to the `writeVcf` function to be written to a
file.
