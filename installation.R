requirements <- c(
    "XML",
    "base64enc",
    "png",
    "plyr",
    "ggplot2",
    "dplyr",
    "jsonlite",
    "httr",
    "summarytools",
    "RColorBrewer",
    "ggsignif",
    "tidyverse",
    "MBESS",
    "readxl",
    "gridExtra",
    "ggprism",
    "compareGroups",
    "stringr",
    "ggpubr"
)

install.packages("https://github.com/jwfoley/bioanalyzeR/releases/download/v0.9.0/bioanalyzeR_0.9.0.tar.gz", repos = NULL)
install.packages(requirements)
# Tested with R 4.2.0
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
install.packages("dplyr")
install.packages("ggprism")
