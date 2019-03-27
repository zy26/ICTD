

#2018-01-23
#wnchang@iu.edu


coexpression_WGCNA_0.2 <- function(data_m, cancerType)
{

print( dim(data_m) )

cancerStr <- cancerType
print(cancerStr)

#==========================first, log data===========================================

datExpr <- t(data_m) # gene should be colnames, sample should be rownames
#datExpr <- log(datExpr + 1)
dim(datExpr)


#=====================================================================================
#
#  Code chunk 1 : basic setting of WGCNA
#
#=====================================================================================

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.1
workingDir = ".";
#setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()

#=====================================================================================
#
#  Code chunk 2 : choose the power
# 
#  Here, we choose a fixed power, thus omit this part
#=====================================================================================


## Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
## Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
## Plot the results:
sizeGrWindow(12, 9)
##pdf(file = "COAD.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
## Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
## this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
## Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3 : cal net
#                 in this part, we need to pay attention to parameters which are
#                 'power', 'minModuliSize', and 'mergeCutHeight'
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,     #30,
                       reassignThreshold = 0, mergeCutHeight =  0.25,    # 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 4 : plot clustring
#
#  Here, annotate plotDendro because we want to turn off visual function
#=====================================================================================

title_str <- paste(cancerStr, ", power=9,minModuliSize=30,mergeCutHeight=0.25", sep="")

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main =title_str)

#=====================================================================================
#
#  Code chunk 5 : output and check
#
#=====================================================================================

table(net$colors)
#table(mergedColors)

colorNumber <- length(table(mergedColors))

#'1' means turquoise, biggest cluster	
#'0' more like noise, we ignore it. The rest of group(different color) is 	

data_coexp <- list()
length(data_coexp)

for (i in 1: (colorNumber-1)  )
{		
	#print(i)
	data_class_trans <- datExpr[ ,which(net$colors == i)]
	data_class <- t(data_class_trans)
	####print( length(data_coexp) )
	data_coexp[[length(data_coexp)+1]] <- data_class	
}
for (i in 1: (colorNumber-1) )
{
	list_str <- paste('Top_', i, sep = "")
	names(data_coexp)[i] <- list_str
}

#data_class1_trans <- datExpr[,which(net$colors == 1)]
#data_class1 <- t(data_class1_trans)

#data_class2_trans <- datExpr[,which(net$colors == 2)]
#data_class2 <- t(data_class2_trans)

#data_class3_trans <- datExpr[,which(net$colors == 3)]
#data_class3 <- t(data_class3_trans)

#data_class4_trans <- datExpr[,which(net$colors == 4)]
#data_class4 <- t(data_class4_trans)

#data_class5_trans <- datExpr[,which(net$colors == 5)]
#data_class5 <- t(data_class5_trans)

##using rbind
#data_5class <- rbind(data_class1, data_class2, data_class3, data_class4, data_class5)
#####print( length(data_coexp) )
#data_coexp[[length(data_coexp)+1]] <- data_5class
#names(data_coexp)[6] <- "mix5"

#file_str <- paste(cancerStr, "_coexp_data.RData", sep="")
#save(data_coexp, file = file_str)

####print("saved top * cluster and mix5 coexpression data to local workspace")

return(data_coexp)


}	#parenthesis for function