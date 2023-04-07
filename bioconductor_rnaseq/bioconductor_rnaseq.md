## RNA-seq workflow: gene-level exploratory analysis and differential expression

In this tutorial, we use a [bioconductor docker image](https://www.bioconductor.org/help/docker/) that includes [asciinema](https://asciinema.org/) functionality to run the [RNA-seq tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html). Here we walk through an end-to-end gene-level RNA-seq differential expression workflow using Bioconductor packages. We will start from the FASTQ files, show how these were quantified to the reference transcripts, and prepare gene-level count datasets for downstream analysis. We will perform exploratory data analysis (EDA) for quality assessment and to explore the relationship between samples, perform differential gene expression analysis, and visually explore the results.

In the examples below, `$` indicates the command line prompt within the container.

<!-- blank line -->
----
<!-- blank line -->

## Learning Objectives:
 - pull a [bioconductor docker image](https://hub.docker.com/r/bioconductor/bioconductor_docker) from DockerHub
 - run the bioconductor/bioconductor_docker container via docker
 - perform exploratory data analysis (EDA) on RNA-seq data

## Assignment 
1. Complete the assignment described below.
2. Upload a link to your dockerhub account.
3. Upload a html or word document with plots.
4. Upload a link with asciinema screen-cast.


### Prerequisites
* create an asciinema account using email at [asccinema.org](https://asciinema.org/login/new) 
* navigate to the your project directory: ~/path/to/project/
<!-- blank line -->
----
<!-- blank line -->

 ### Assignment Points
|  Rubric        | Points | 
|----------------|-------|
| Screencast     |  -/50  |
| Plots          |  -/50 |
| On Time        |  -/30  |
*Total Points: -/130*

## Getting Started

### 1. open docker teminal

![asciinema_auth](https://github.com/GMS6804-master/assignment/blob/main/images/terminal_start.png)
<!-- blank line -->
----
<!-- blank line -->

### 2. pull a docker image from DockerHub
```
docker pull [YOUR DOCKERHUB ID]/bioconductor_rnaseq:[month_year]
```

### 3. boot into container as bash while also mounting a "dropbox-style" directory that will link your docker container to your local machine
```
docker run -it -v [path-to-working-directory]:/projeect dominicklemas/bioconductor_rnaseq:04_2023 bash
```
As an example: 
```
docker run -it -v C:/Users/djlemas/OneDrive/Documents/rna-seq:/project dominicklemas/bioconductor_rnaseq:04_2023 bash
```
<!-- blank line -->
----
<!-- blank line -->

### 4. link your container to your asciinema.org account by opening the URL in a web browser 
```
asciinema auth
```
![asciinema_auth](https://github.com/GMS6804-master/assignment/blob/main/images/asciinema_auth.png)
<!-- blank line -->
----
<!-- blank line -->

### 5. add screen-cast headers 
```
asciinema rec
# Name: 
# Date: 
# bioconductor:: bioconductor_rnaseq
```
<!-- blank line -->
----
<!-- blank line -->

### 6. start R, 
```
R
```

### 7. Begin the tutorial at 2.3 Reading in data with tximeta

First 7 lines of code:
```
library("airway")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata
```
continue through the rest of the tutorial. 

### 8. Exporting Figures

Exporting figures via command-line requires a few more R functions. Note- the example below has a function called png(file="file_name.png") that comes before the plotting function and the dev.off() directly after the plotting function. 

```
# 4 Exploratory analysis and visualization
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)

# Figure_4.2.1
png(filename = "Figure_4.2.1.png");
meanSdPlot(cts, ranks = FALSE);
dev.off()
```
