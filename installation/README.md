# Installing ARNAQ

## Installing dependencies

Apart from the R package requirements, ARNAQ also needs [pandoc](https://pandoc.org)
on the system. If you are using Linux, it is simple and convenient to use
[conda](https://anaconda.org/anaconda/conda) to 
create an environment with `pandoc` and the required R packages.

### Installing dependencies with conda

The file `ARNAQ_env.yml` can be used with `conda` to set up a suitable environment. The enviroment
file can be found [here](https://github.com/Queex/arnaq/blob/main/inst/extdata/ARNAQ_env_1.0.yml).
From the command line:

```{bash conda_1, eval=FALSE}
conda env create -f ARNAQ_env_1.0.yml
conda activate arnaq_1.0
```

When you are done, you can use the usual conda command to return to the standard environment.

```{bash conda_2, eval=FALSE}
conda deactivate
```

### Installing dependencies manually

Several of ARNAQ's required packages are not available for Windows through bioconda. For Windows
installation, or if you do not use conda on your system, a manual installation needs to be done.

First install pandoc, from [here](https://pandoc.org/installing.html).

Then make sure the `devtools` package is installed:

```{r pre_install, eval=FALSE}
install.packages("devtools")
```

#### Caveats

As of writing this guide, several dependencies of dependencies do not seem to be available for R 
version 4.4 via
this install method. You can try installing them from Bioconductor manually, or use a previous
version of R. If you follow the latter method, one package requires R 4.4 and up, which will raise
a warning, but this does not seem to prevent ARNAQ from running.

## Installing ARNAQ from github

```{r install, eval=FALSE}
devtools::install_github("Queex/arnaq", build_vignettes=TRUE)
```

Despite the "helpful" message, it is best not to update existing packages when prompted as conda
will have selected versions that are properly compatible with the other packages.

This command will also build the vignettes, which can take a few minutes to complete as they
include the generation of reports from exmaple data included in the package. If you would prefer
not to wait, set `build_vignettes` to `FALSE`.

## Support Files

ARNAQ uses two support files to process the data and create reports, and one file specifying 
metadata about the project. By default, ARNAQ will use internal versions of the first two. For the
third, you will need to create it for every project. You can create this file manually, copy over 
and edit one from a previous project, or create a blank template using the following
command:

```{r prepare, eval=FALSE}
arnaq.prepare("resources")
```

The same function can also be used to create copies of the other two files, if you need to alter
them.

These are stored as files outside R in order to make it possible for them to be automatically
generated (or at least generated in part) by an upstream pipeline.

## External Support Files (.gtf)

*Optional.*

ARNAQ can make use of an external data file, available from public repositories, for each reference
genome you intend to use ARNAQ with. They are not included in ARNAQ's distribution, as they must be
specific for your mapping pipeline and match the gene identifiers in your count files.

Specifically, ARNAQ presumes that Ensembl `.gtf` files for the organism are available. These should
include identifiers drawn from the same set as the gene identifiers uses in the count tables. There
should be a 'type' field, so ARNAQ can identify the gene-level entries. For the biotype portions of
ARNAQ output, there should also be a 'gene_biotype' field. You may wish to filter the Ensembl files
to just the entries with a `type` of `"gene"`, to reduce the size of the file and make ARNAQ runs
faster, but it is not required.

You will only need one such file for each reference genome you are working with.

### Directory Structure

ARNAQ presumes a particular directory structure for these files.

`resource_dir` and `genome_reference` in `resources.yml` define the directory ARNAQ looks in for
the `.gtf` file. The intention is that `resource_dir` gives the directory where a number of named
genomes are placed, and can be kept the same for every project. `genome_reference` names the 
directory in that location where the correct `.gtf` file can be found. ARNAQ will use the first
`.gtf` file it finds in that directory, so you don't have to rename files you download and risk
losing information as to their version.

If your directory looks like this:

```{text}
/some/directory/you/have/
├── GRCh38/
│   └── GRCh38.p14.20240712.gtf
└── GRCm39/
    └── GrCm39.20240317.gtf
```

The matching lines in `resources.yml` would be:

```{text}
resource_dir: /some/directory/you/have/
genome_reference: GRCh38
```

or

```{text}
resource_dir: /some/directory/you/have/
genome_reference: GRCm39
```

The directory path should be an absolute path; so you don't have to adjust it between projects in
different directories.

## ERCC Reference

*Optional.*

If you intend to use the ERCC spike-in functionality of ARNAQ, you will need to download the table
of ERCC ids and mix concentrations from ThermoFisher Scientific. As of this revision, this can be
found at [https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).

Make sure the relevant line in `resources.yml` points to this file:

```{text}
ercc_concentrations: path/to/file/from/TFS/filename.txt
```
