# R script with functions used in Define Housekeeping Genes script

# Libraries required for analysis
library(reldist)
library(tidyverse)
library(ComplexHeatmap)
library(UpSetR)

# Function to calculate gini for gene expression data

gene.Gini <- function(x){
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x!=0))
      {
        result <- gini(x)
      } else {
        result <- 0
      }
    } else {
      result <- NA
    }	 		
  } else {
    result <- NA
  }
  return(result)
} 

# Function to standardize TPM score decimal place between data sets

round_number <- function(x,decimal){
  numeric_values <- sapply(x, mode) == 'numeric'
  x[numeric_values] <-  round(x[numeric_values], decimal)
  x
}

# Function to calculate percent of genes with TPM >10 and gini score 

# Genes with >10 TPM in over 90% of tissues and gini score <= 0.4 are considered housekeeping

define.HK.genes <- function(expression.file,no.columns,no.tissues){
  Expression <- read.csv(expression.file)
  
  Expression[is.na(Expression)] <- 0
  
  Expression <- round_number(Expression,1)
  
  genes <- Expression$Name
  
  expression.score <- Expression[,3:no.columns]  
  
  inverted <- as.data.frame(t(Expression))
  
  gene.exp <- inverted[c(-1,-2),]
  
  names(gene.exp) <- genes
  
  percentage <- data.frame()
  
  for(i in genes) {
    intermediate <- gene.exp[,i]
    inter.sum <- sum(as.numeric(intermediate) > 10)/no.tissues*100 # Calculate percent of tissues with TPM >10 
    percentage <- rbind(percentage,inter.sum) }
  
  names(percentage) <- "percent"
  
  with.percent <- cbind(Expression,percentage)
  
  score <- as.data.frame(t(expression.score))
  
  names(score) <- genes
  
  gini.score <- data.frame()
  
  for(i in 1:ncol(score)) {       
    gene.x <- score[i] 
    lc.order <- gene.x %>% arrange(gene.x) 
    names(lc.order) <- "mid"
    lc.data <- lc.order$mid
    gin <- gene.Gini(lc.data) #calculate gini score
    gini.score <- rbind(gini.score,gin)
  }
  
  names(gini.score) <- "equality"
  
  with.equality <- cbind(with.percent,gini.score) # Table of gene expression data with the percent of tissues >10 TPM calculated and gini score values
  
  housekeeping <- with.equality[with.equality$percent >= 90 & with.equality$equality <= 0.4,] # Select genes with 90% of tissues with >10TPM and <=0.4 gini score
  
  Housekeeping.genes <- housekeeping$Name
  
  Housekeeping.genes # Return list of housekeeping genes
}
