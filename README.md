# Description

This method is described in the publication from Biorxiv, 2018 available at https://www.biorxiv.org/content/10.1101/426593v2

The web application demo is available at https://ictd.ccbb.iupui.edu

![image](https://github.com/zy26/ICTD/blob/master/img/web_app.png)

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



Note : For old R version which cannot install 'BiocManager', please use below command to install the dependency.
```
source("https://bioconductor.org/biocLite.R")
biocLite("impute")
biocLite("GO.db")
biocLite("sva")
biocLite("preprocessCore")

rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
```

# Example

```
library(ICTD)

data_bulk = GSE72056_diri_example[[1]]
ictd_result <- ICTD(data_bulk)

#Return value is a list, which the first element is the predicted proportion and 
#the second element is the predicted markers of ICTD
```
# Dependencies

We also provide a Docker image to recreate the compute environment. See the Dockerfile for more details.

https://cloud.docker.com/u/wnchang/repository/docker/wnchang/ictd

Using the Docker image could void the conflict issue that R version and several R packages version confict. 

For more details about the Docker, please see Docker documentation page https://docs.docker.com/.
