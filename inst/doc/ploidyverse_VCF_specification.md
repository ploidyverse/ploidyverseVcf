# Variant Call Format specification for the ploidyverse

The Variant Call Format, or VCF, is the standard in bioinformatics for storing
information about DNA sequence variation.  One VCF file can contain genotypes
across many loci and individuals.  Because VCF is highly flexible, widely used,
and self-documenting, we adopt it here as the primary format for transferring
data among ploidyverse packages, and between the ploidyverse and other
software.

VCF is stored as tab-delimited text, with a corresponding binary format called
BCF.  Here we describe guidelines for VCF files imported into or
exported from ploidyverse packages, and generally set standards for storing
genotype data from polyploid organisms.  These guidelines are meant to make it
easier to transfer data among software and researchers, rather than impose
restrictions.  If your needs deviate in some way from what is described here
and you would like to suggest some changes, you may file an issue or pull
request on GitHub, or start a discussion on the ploidyverse Google Group.

Additionally, the 
[`VariantAnnotation`](https://doi.org/doi:10.18129/B9.bioc.VariantAnnotation) 
package, from the Bioconductor collection of R packages, implements flexible 
S4 classes for storing any or all data from a VCF file within an R object, and
has functions for reading and writing VCF files.  The `ploidyverseVcf` R 
package extends `VariantAnnotation` and is designed to assist package 
developers and users with meeting the ploidyverse VCF guidelines.

## Guiding principles

VCF data within the ploidyverse should:

* Be fully compatible with the 
[current VCF specification](https://github.com/samtools/hts-specs).
* Be self-documenting.  Someone 50 years from now should be able to open up
the file and understand what is going on and how to use the data.
* Be adapted for use with non-model organisms.  There may be no reference
genome, or a reference genome for a different species from the one being
studied.
* Provide allelic read depths and posterior probability estimates for all
genotype calls.  Allele dosage will be uncertain without very high depth 
sequencing, and that uncertainty must be quantified in the file.

## Software information

For VCFs containing genotype calls, it is important to include information
on what software and model was used for generating those calls.  ploidyverse
VCFs should include a header formatted as follows:

```
##ploidyverse=<ID=GenotypeCalls,Software=EBG,Version=1,Model=hwe,Description="Genotype calls generated from EBG">
```

Additional lines can be added for other pieces of software used in the
pipeline for generating the VCF.  Lines with the following `ID` values are 
suggested:

* `SNPdiscovery`: How variant loci were identified.
* `ReadDepth`: How allelic read depth was calculated for each sample and locus.
* `GenotypeCalls`: How genotypes and genotype posterior probabilities were
inferred from allelic read depth.
* `Phasing`: How alleles were arranged into haplotypes.
* `Imputation`: How missing genotypes were imputed.

In `ploidyverseVcf`, the `software` function can be used to add these lines.

## Reference genome

If you aren't studying a model organism, it is likely that you are using a
reference genome that is likely to be updated in the next few years, and/or
a reference genome for a related species rather than the actual species that
you are studying.  Therefore it is important to indicate what reference was
used.

Lines beginning with `##contig` provide information about each contig (usually,
each chromosome) that was used for alignment.  Below are examples of contig 
lines that should be found in a ploidyverse VCF.

```
##contig=<ID=Chr01,length=150798994,assembly="Miscanthus sinensis v7.1 on Phytozome">
##contig=<ID=Chr02,length=146084461,assembly="Miscanthus sinensis v7.1 on Phytozome">
```

Optionally, you may wish to indicate the file name for the reference genome
with a line formatted as follows:

```
##reference=Msinensis_DH1_v7_0.fasta
```

If you are using a non-reference pipeline, include the following line:

```
##contig=<ID=NonRef>
```

And use "NonRef" in the CHROM column for all variants in the genotype table.

In `VariantAnnotation`, data from contig lines can be retrieved or assigned with 
the `seqinfo` function.

``` r
seqnames(seqinfo(vcf))   # The 'ID' field
seqlengths(seqinfo(vcf)) # The 'length' field
genome(seqinfo(vcf))     # The 'assembly' field
```

A URL and MD5 checksum can also be specified in the contig lines, although I am
not sure that these are currently supported by `VariantAnnotation`.

## SNP information

The `INFO` column of the genotypes table can be used to store
information about each variant in addition to what is already found in the
`CHROM`, `POS`, `REF`, `ALT`, `QUAL`, and `FILTER` columns.  Header lines 
beginning with `##INFO` describe all additional fields in the `INFO`
column.

When there is no reference genome, ploidyverse VCFs require an `INFO` field
called `TAG` that provides the sequence tag flanking the site.  When there
is a reference genome, we still recommend this field but do not require it.
For reduced representation technologies such as genotyping-by-sequencing (GBS)
or restriction site-associated DNA sequencing (RADseq), the entire sequence
tag beginning with the restriction cut site should be provided.  Either the
reference allele or most common allele should be represented here.  The
purpose is to facilitate alignment to future reference genomes and 
cross-referencing to other projects using the same reduced representation 
technology.

The following lines define these `INFO` fields:

```
##INFO=<ID=TAG,Number=1,Type=String,Description="DNA sequence of reference RAD tag">
##INFO=<ID=TAGREVCOMPL,Number=0,Type=Flag,Description="The tag aligns to the bottom strand of the reference genome">
##INFO=<ID=TAGPOS,Number=1,Type=Integer,Description="Position of variant with respect to the beginning of the RAD tag">
```

*Need to add details about how to infer the variant tags from this information*

## Sample information

In VCF, sample names are used for column headers in the genotype table.
It is also possible to store more data about samples in header lines that
begin with `##SAMPLE`.  (See "Sample field format" in the VCF specification.)
ploidyverse VCFs require the following information for each sample:

* **Species**.  Spell out the full genus and species name.  Optionally, add 
subspecies, variety, cultivar, and/or authority if necessary to unambiguously
identify the species.
* **Ploidy**.  If a reference genome was used, provide ploidy with respect
to the reference genome.  Otherwise, provide ploidy with respect to the
monoploid genome.  Use "2x" to represent a diploid, "4x" to represent an
autotetraploid, "2x+2x" to represent an allotetraploid, etc.  Note that
allohexaploid wheat, for example, would be listed as "2x" if hexaploid
wheat reference were used, but "2x+2x+2x" if the reference of a diploid
ancester were used.

The following header lines define the above two fields and are mandatory:

```
##META=<ID=Species,Type=String,Number=.,Description="Species name">
##META=<ID=Ploidy,Type=String,Number=.,Description="Ploidy with respect to reference genome">
```

Optional fields:

* **CollectionLocation**.  Provide decimal latitude and longitude, and/or
a description of the collection location of wild materials.
* **MaterialProvider**.  Give the name of the organization or person curating
the biological material who supplied it for this project.
* **PopulationDesign**.  If the samples in the VCF represent an artificial 
population (e.g. from a controlled cross), use this field to identify 
parents, progeny, controls, and/or family membership.  (See also the Pedigree
field in the VCF specification.)
* **BioSample**.  Provide the NCBI BioSample accession number in the format
`SAMN00000000`.
* **DOI**.  A digital object identifier uniquely identifying the sample, and
redirecting to a URL where more information about the sample is available.
Use either `doi:xxxxx` or `https://doi.org/xxxxx` format.
* Any custom fields may be added as long as `##META` lines are included that
contain descriptive definitions of those fields.  If you have pertinent 
information about your samples stored in a spreadsheet, it should go here.

Below are suggested header lines for these optional fields:

```
##META=<ID=CollectionLocation,Type=String,Number=.,Description="Collection location of wild samples">
##META=<ID=MaterialProvider,Type=String,Number=.,Description="Person or organization curating biological materials">
##META=<ID=PopulationDesign,Type=String,Number=.,Values=[Mother, Father, Progeny, Control],Description="Identity within artificial population">
##META=<ID=BioSample,Type=String,Number=.,Description="Accession number in NCBI BioSample database">
##META=<ID=DOI,Type=String,Number=.,Description="Digital Object Identifier">
```

And below are examples of what the sample lines might look like:

```
##SAMPLE=<ID=PMS-014,Species=Miscanthus sinensis,Ploidy=2x,BioSample=SAMN02213066,CollectionLocation=29.657 N 109.1195 E,Description="Miscanthus sinensis 'PMS-014'">
##SAMPLE=<ID=JM2014-S-4,Species=Miscanthus x giganteus,Ploidy=3x,BioSample=SAMN05752168,CollectionLocation=34.3014 N 134.1304 E,Description="Miscanthus x giganteus 'JM2014-S-4'">
```

In the above example, "PMS-014" and "JM2014-S-4" would also be found as column
headers for the genotype table.

Within the `VariantAnnotation` package, the table of sample information can be 
accessed by:

``` r
meta(header(vcf))$SAMPLE
```

or, using the `ploidyverseVcf` package:

``` r
sampleinfo(vcf)
```

where `vcf` is the name of a vcf object.

## Genotypes

Within the tab-delimited genotype table, for every sample by locus, multiple
pieces of information can be stored, separated by a colon (`:`).  The `FORMAT`
column indicates which fields are stored in the table and what order they are
in.  See the official VCF specification for more information.

In `VariantAnnotation`, data from the genotype table can be retrieved using the
`geno` function, which returns a named list of matrices.  Names of the list
correspond to the genotype fields found in the file.  Each matrix stores data
from the table for that field.  For fields such as `AD` or `GP` that can store
multiple comma-separated values, the matrix is a two-dimensional list of 
vectors.

### GT field

In keeping with the VCF specification and with existing software for calling
genotypes in polyploids, ploidyverse VCFs code unphased genotypes in the format
`0/0/0/1` and phased genotypes in the format `0|0|0|1`.  The example given is a
tetraploid with three copies of the reference allele and one copy of the first
alternative allele.  The genotype with the highest posterior probability should
be given in the GT field.

### AD field

This is the only genotype format field that is mandatory for VCFs that are to
be used as input to a ploidyverse genotype calling pipeline.  It is also
mandatory for VCFs output by the genotype calling pipelines, in order to enable
independent evaluation of genotyping quality.

This field contains integers separated by commas.  The first integer corresponds
to the reference allele, the second integer to the first alternative allele,
the third integer to the second alternative allele, and so on.  Every integer
indicates the read depth for the corresponding allele.  So, if the reference 
allele is A and there is only one alternative allele G, a value of `12,4`
would indicate 12 reads of A and four reads of G.  In a tetraploid that value 
in the AD field might correspond to the genotype `0/0/0/1`.

### GP field

This field contains genotype posterior probabilities, ranging from zero to one
and summing to one, for all possible genotypes.  These are comma-separated.
In order to conserve disk space, we recommend rounding to the nearest 0.001.

See the VCF specification for the ordering of genotypes.  For example, in a
tetraploid with one alternative allele, there are five possible genotypes,
and a posterior probability must be provided for each one.

If genotypes were called with a likelihood approach rather than a Bayesian
approach, uniform priors can be used to convert likelihoods to posterior
probabilities.

### PS field

This field is needed only if SNPs are phased.  If the data were produced by
RAD-seq or similar reduced-representation DNA sequencing technology, we
highly recommend phasing SNPs that were called within the span of one RAD
tag, *i.e.* alleles that can be unequivocally phased because they were on
the same sequencing read.

As suggested in the VCF specification, the phase set ID should be the
alignment position of the first SNP in the phase set.  For non-reference
pipelines, phase set IDs can simply be numbered sequentially from one.

### GN field

This is a custom field not described in the VCF specification, containing a
numeric genotype ranging from 0 to 1, indicating the estimated intra-individual
allele frequency for the alternative allele.  For example, if the genotype
`0/0/0/1` was called with high confidence, the value would be 0.25.

Generally this value should be the posterior mean genotype, *i.e.* the 
copy number of the alternative allele multiplied by the posterior probability
of that genotype, summed across all possible genotypes, divided by the ploidy.
If the value for this field is calculated differently, it should be indicated
in the file header.

This field should have one value for each alternative allele, separated by
commas when there are multiple alternative alleles.
