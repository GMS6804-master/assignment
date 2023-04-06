# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

RUN apt-get update

### Bioconductor: Part 1
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("DESeq2", \
    "tximeta", "airway", ask = F))' \
	&& R -e 'tinytex::install_tinytex()' \
	&& R -e 'install.packages("dplyr")'

### Bioconductor: Part 2
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("vsn", \
    "Gviz", "pheatmap", "PoiClaClu", "genefilter", "IHW", "org.Hs.eg.db", \
	"ReportingTools", "sva", "RUVSeq", "fission", \
	"glmpca", "ggbeeswarm", "apeglm", ask = F))' 

### Install asciinema
RUN apt-get install -y  \
asciinema 

### Set Working Directory
WORKDIR /project