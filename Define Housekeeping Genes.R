# R script to determine Housekeeping genes in 4 species and categorise ancestral, ruminant specific, and cattle specific housekeeping genes

source("Housekeeping gene functions.R")


# Filter gene expression tables so they contain only similar tissue types

# Create list of housekeeping genes for each species:

human <- define.HK.genes("path/to/filtered.human.expression.table.csv",26,24) # Human housekeeping genes list 

sheep <- define.HK.genes("path/to/filtered.sheep.expression.table.csv",15,13) # Sheep housekeeping genes list

pig <- define.HK.genes("path/to/filtered.pig.expression.table.csv",19,17) # Pig housekeeping genes list 

cow <- define.HK.genes("path/to/cow.expression.table.csv",11,9) # Cow housekeeping genes list 



# Human GTEx data uses gene ID version - to standardise remove the version and use only gene stable ID 

human.hk <- gsub("[[:punct:]]", " ", human)

human.hk <- as.data.frame(stringr::str_split_fixed(human.hk," ", 2))

human.hk <- human.hk$V1



# Orthologue tables from Ensembl Biomart - homologue type = ortholog_one2one

human.orthologues <- read.csv("path/to/human.orthologues.table.txt") # Ensembl orthologue table with human as reference

cow.orthologues <- read.csv("path/to/cow.orthologues.table.txt") # Ensembl orthologue table with cow as reference



# Expression atlas and Ensembl sheep gene IDs do not match, the following code fixes this: 

human.orthologues$Sheep.gene.stable.ID <- sub('020','000',human.orthologues$Sheep.gene.stable.ID)

cow.orthologues$Sheep.gene.stable.ID <- sub('020','000',cow.orthologues$Sheep.gene.stable.ID)



# Defining ancestral Housekeeping genes

shared.orthologues <- na.omit(human.orthologues) # Orthologues shared in all 4 species 

ancestral.HK.genes <- shared.orthologues[shared.orthologues$Gene.stable.ID %in% human.hk & 
                                           shared.orthologues$Cow.gene.stable.ID %in% cow & 
                                           shared.orthologues$Sheep.gene.stable.ID %in% sheep & 
                                           shared.orthologues$Pig.gene.stable.ID %in% pig,]  # Orthologues that are housekeeping in all species 



# Defining ruminant specific  Housekeeping genes 

ruminant.HK.genes <-  cow.orthologues[cow.orthologues$Gene.stable.ID %in% cow &
                                        cow.orthologues$Sheep.gene.stable.ID %in% sheep &
                                        !(cow.orthologues$Pig.gene.stable.ID %in% pig) &
                                        !(cow.orthologues$Human.gene.stable.ID %in% human.hk),] # Orthologues that are housekeeping in cow and sheep but not human or pig 



# Defining cattle specific: 

# Any cattle housekeeping genes not in the ancestral or ruminant specific housekeeping genes and not housekeeping in any other species is cattle specific
# Use similar code as ancestral/ruminant and go through each combination of shared housekeeping genes -> removing any housekeeping shared with other species 



# Upset R plot 

# Required a common identifier for shared orthologues - as species gene IDs don't overlap a new identifier needs creating:

# For each orthologue combination create a new common identifier

# Create a list of all species new common identifiers - list.common.identifiers

matrix.identifiers <- as.data.frame(list_to_matrix(list.common.identifiers))

comb.identifiers <- make_comb_mat(matrix.identifiers)

# Make upsetR plot

upset(comb.identifiers, keep.order = T, 
      sets = c("human set name in list","pig set name in list","sheep set name in list","cow set name in list"),
      order.by = "freq", query.legend = "top",queries = list(list(query = intersects, params = list("sheep set name in list","cow set name in list"), color = "darkgreen", active = T, query.name = "Ruminant specific"),
                                                             (list(query = intersects, params = list("cow set name in list","sheep set name in list","human set name in list","pig set name in list"), color = "maroon", active = T, query.name = "Ancestral"))),
      set_size.show = T)

