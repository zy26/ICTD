# Description

This method is described in the publication from Biorxiv, 2018 available at https://www.biorxiv.org/content/10.1101/426593v2

The web application demo is available at https://wnchang.shinyapps.io/ICTD_server/


# Installation


```
#install dependent pkg
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("impute", version = "3.8")
BiocManager::install("GO.db", version = "3.8")
BiocManager::install("sva", version = "3.8")
BiocManager::install("preprocessCore", version = "3.8")

rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

#install ICTD
install.packages("devtools")
devtools::install_github("zy26/ICTD")
```

# Example

```
library(ICTD)

data_bulk = GSE72056_diri_example[[1]]
ictd_result <- ICTD(data_bulk)

#Return value is a list, which the first element is the predicted proportion and 
#the second element is the predicted markers of ICTD
```


