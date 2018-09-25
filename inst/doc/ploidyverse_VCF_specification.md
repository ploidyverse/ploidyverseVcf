# Variant Call Format specification for the ploidyverse

The Variant Call Format, or VCF, is the standard in bioinformatics for storing
information about DNA sequence variation.  One VCF file can contain genotypes
across many loci and individuals.  Because VCF is highly flexible, widely used,
and self-documenting, we adopt it here as the primary format for transferring
data among ploidyverse packages, and between the ploidyverse and other
software.

VCF is stored as tab-delimited text, with a corresponding binary format called
BCF.  Additionally, the 
[VariantAnnotation](https://doi.org/doi:10.18129/B9.bioc.VariantAnnotation) 
package from the Bioconductor 
collection of R packages implements flexible S4 classes for storing any or all
data from a VCF file within an R object, and has functions for reading and 
writing VCF files.  In the ploidyverse, we extend these S4 classes, adding
requirements for data that are mandatory for ploidyverse packages, and
generally setting standards for storing genotype data from polyploid organisms.

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
##ploidyverse=<ID=GenotypeTable,Software="EBG",Version=1,Model="hwe">
```

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
##META=<ID=DOI,Type=String,Number=.,Descrption="Digital Object Identifier">
```

And below are examples of what the sample lines might look like:

```
##SAMPLE=<ID=PMS-014,Species=Miscanthus sinensis,Ploidy=2x,BioSample=SAMN02213066,CollectionLocation=29.657 N 109.1195 E,Description="Miscanthus sinensis 'PMS-014'">
##SAMPLE=<ID=JM2014-S-4,Species=Miscanthus x giganteus,Ploidy=3x,BioSample=SAMN05752168,CollectionLocation=34.3014 N 134.1304 E,Description="Miscanthus x giganteus 'JM2014-S-4'">
```

In the above example, "PMS-014" and "JM2014-S-4" would also be found as column
headers for the genotype table.

Within the VariantAnnotation package, the table of sample information can be accessed by:

```{r}
meta(header(vcf))$SAMPLE
```

where `vcf` is the name of a vcf object.
