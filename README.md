# Description

This method is described in the publication from Biorxiv, 2018 available at https://www.biorxiv.org/content/10.1101/426593v1


# Installation


```
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


