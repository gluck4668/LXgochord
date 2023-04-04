\name{LXgochord}
\alias{LXgochord}
\title{LXgochord}
\usage{
LXgochord(gene_file,species,term_number)
}
\description{
LXgochord can be used to show the chord diagram of GO enrichment.
}
\examples{

if(!requireNamespace("devtools"))
  install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgochord")

library(LXgochord)

??LXgochord

#------------------------------------
data(gene_data_example)

#------------------------------------

gene_file <- "gene_data.xlsx" # The file should include three columns of gene_symbol,log2FC (or log10FC) and pvalue.

species <- "Mouse" # It should be "human","rat", or "mouse“.

term_number=c("T1157","T1272","t572") # "term_number" can be selected from "term number.xlsx", and it should be more than two.

LXgochord(gene_file,species,term_number)

}
