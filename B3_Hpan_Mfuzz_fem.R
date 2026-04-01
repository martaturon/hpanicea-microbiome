
#title: "Mfuzz tutorial"
#subtitle: "Whole transcriptome data -- Hpan reproduction"
#First set working directory

#Create in your working directory a folder for images and another for stats where you can save all the files you generate
setwd(file.path(wd))
imgdir  <- "./images/"
statdir <- "./stats/"

## Intalling and loading libraries

#We use mfuzz through bioconductor. If you don't have bioconductor installed install it using the following code. Do not run the code if you have it installed since it can be little bit slow.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#If not, install mfuzz. If you have it already installed please skip this step because it can be time consuming.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Mfuzz")
#Load the necessary libraries. If they are not installed install them using `install.packages`
library("Mfuzz")
library("Biobase")
library("kableExtra")
library("dplyr")

#if needed:
#install.packages("svglite")
#install.packages("kableExtra") 

## MFuzz in RNAseq data

### Import data

#We need a table of normalized counts. We just load a table delimited with tabulator to R with the data. Mfuzz cannot deal with replicates so we need a single column per condition, and ordered however you want them to appear (chronologically if it is a time series, etc)
rnaseq_counts <-read.table("Hpan_final.isoform.counts.matrix_BLASTIDMetazoa_females_elements.txt", sep="\t", header = TRUE, row.names = 1)
head(rnaseq_counts, 5)
str(rnaseq_counts)

#We define some of the stages that we have in the analysis in a vector for easy reproducibility of the analysis.
#For females
stages <- c("Mar_OO","Apr_OO","May_lateOO", "May_EE", "May_ME", "Jun_Lar")


#For males
#stages <- c("Mar_NR","Mar_OO","Apr_OO","Apr_SP","May_OOE", "May_SP", "Jun_NR", "Jun_Lar", "Jun_FL")

### Prepare data

#We need transform the data to a matrix because mfuzz can not work with a data frame, and needs the `ExpressionSet` to work
matrix_rnaseq_counts <- as.matrix(rnaseq_counts)

#We save the data in a new object using `ExpressionSet`
normalizedSet <- ExpressionSet(assayData=matrix_rnaseq_counts)
normalizedSet
#When looking at the normalizedSet:
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 179328 features, 6 samples 
#element names: exprs 
#protocolData: none
#phenoData: none
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation:  

#Mfuzz don't like NA or missing values, so we have to filter any possible NA values. We also standardize the expression values so that the average expression value for each gene is zero and the standard deviation of its expression profile is one. This step ensures that vectors of genes with similar changes in expression are close in Euclidean space
data.r <- filter.NA(normalizedSet, thres=0.5) #based on NA missing values
data.k <- fill.NA(data.r,mode="knn") #replace remaining missing values by the average values expression value of the corresponding gene. k-nearest neighbour method can be used (mode=`knn'/`wknn').

#Note that you won't get NA values when using edgeR, those are DESeq2 outputs. However, you may want to remove genes which are expressed at low levels or show only small changes in expression. For that include the filtering step: filter.std --\> exclude genes with low standard deviation
data.f <- filter.std(data.k, min.std=0.001, visu=TRUE) #based on standard deviation. 
#If visu is set to TRUE, the ordered standard deviations of genes' expression values will be plotted
#Here the result is that 131511 genes are excluded

#Finally standarize to make transcripts (or genes or proteins) comparable to Mfuzz. Note that the data.s is an object of the class `ExpressionSet`
data.s <- standardise(data.f)
data.s

#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 74023 features, 6 samples  
#element names: exprs 
#protocolData: none
#phenoData: none
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation: 

#### Setting the parameters for FCM clustering

#Usually we will optimize the parameters for the MFuzz analysis for each dataset. This process is memory demanding and for the sake of time we will not be ruining it in this practical.
#The parameters used in Mfuzz are:
  
 # -   **fuzzifier value (m)**: We will determine the optimal setting of fuzzifie value using `mestimate`. This will be only an estimation since the final result tuned depending on the performance of the function.

#As this is a subclustering it can be soft if one gene can be in more than one cluster, or hard if you limit the gene to be only in one cluster. This is determined by the fuzzifier value (m). Under a biological perspective, this soft clustering is more appropiate as one gene can be involved in different functions.

#-   **number of clusters (cluster)**: final number of clusters. This will be determined using repeated soft clustering for a range of cluster numbers. This will report the minimum centroid distance. The minimum centroid distance can be used as cluster validity index. For an optimal cluster number, we may see a 'drop' of minimum centroid distance wh plotted versus a range of cluster number and a slower decrease of the minimum centroid distance for higher cluster number.

#For parameter selection use your `ExpressionSet` object, i.e., data.s:
m1 <- mestimate(data.s)
m1
#m1=1.700821 (for METAZOA_MTX_female)
#m1=1.723651 (for BCT)
#For election of c, either cselection or:
tmp <- Dmin(data.s,m=m1,crange=seq(4,40,4),repeats=3,visu=TRUE) # Note: This calculation might take some time

#The suitable number of clusters might be those that show the decrease for c

#The minimum centroid distance is defined as the minimum distance between two cluster centers produced by the c-means clusterings. The average minimum centroid distance for the given range of cluster number is returned.

#The minimum centroid distance can be used as cluster validity index. For an optimal cluster number, we may see a 'drop' of minimum centroid distance when plotted versus a range of cluster number and a slower decrease of the minimum centroid distance for higher cluster number. Thus, the optimal number of clusters is often chosen at the point where the minimum centroid distance starts to level off, and adding more clusters doesn't provide much improvement in terms of separation

m1=1.700821 #m1=1.694431 (from mestimate)
cluster=7


### Run MFuzz

#We perform mfuzz clustering and represent the plot.

#Note: every time you rerun the mfuzz the order of the clusters may change (but not the structure)
mf <- mfuzz(data.s,c=cluster,m=m1)
mfuzz.plot(data.s,cl=mf, time.labels=stages,new.window=FALSE)


#Colours in the plot represent the belonginess ("membership") to the pattern, e.e., how good or bad each gene (line) fits to the cluster pattern. Green or yellow means they belong less to the clustering. Blue and red colored genes better fit to the cluster dynamics.

#> As you modify the recommended parameters the red parts in the clustering is lost, so we loose reliability.
# Note: there is always important to have red genes that indicate the belonging to the cluster is strong, so good. Thus, it is important to find a balance between m and cluster:
# As it is a soft clustering you always expect some deviation so you need a balance between being so soft (same gene in many clusters) or hard (one gene only in one cluster).
# Also with the number of clusters you need to reach a compromise: less clusters will show a lot of genes assigned to one cluster but with less belonginess whereas many clusters will get genes with low belonginess to the cluster.

#You can try to modify the m and cluster parameters for FCM clustering and how mfuzz clustering change
#For instance:
#m2=1.6
#cluster2=30

#mf2 <- mfuzz(data.s,c=cluster2,m=m2)
#mfuzz.plot(data.s,cl=mf2, time.labels=stages,new.window=FALSE)
```

```{r}
m3=1.5
cluster3=30

mf3 <- mfuzz(data.s,c=cluster3,m=m3)
mfuzz.plot(data.s,cl=mf3, time.labels=stages,new.window=FALSE)
```

```{r}
m4=1.4
cluster4=30

mf4 <- mfuzz(data.s,c=cluster4,m=m4)
mfuzz.plot(data.s,cl=mf4, time.labels=stages,new.window=FALSE)
```

```{r}
m5=1.3
cluster5=30

mf5 <- mfuzz(data.s,c=cluster5,m=m5)
mfuzz.plot(data.s,cl=mf5, time.labels=stages,new.window=FALSE)
```
#to make mfuzz more restrictive you can use lower m values and then get only genes that belong to one cluster
cl <- mfuzz(data.s, c=6, m=1.1)
strong_genes <- acore(data.s, cl, min.acore=0.7)
#### **Save the generated clusters in a PDF**

#NB: I would add in the file name the parameters for FCM clustering

pdf(paste0("mfuzz_rnaseq_m1.7_c7_sponge.pdf"), height=10, width=15)
mfuzz.plot(data.s,cl=mf, time.labels=stages,new.window=FALSE)
dev.off()
```

#### **Extract the genes in each expression cluster**

#The results are stored in the `mf` object. Some interesting metrics are: - centers: the final cluster centers. - size: the number of genes in each cluster of the closest hard clustering. - cluster: indices of the clusters where the genes are assigned to for the closest hard clustering, as obtained by assigning points to the (first) class with maximal membership. - membership: a matrix with the membership values of the data points to the clusters.

#Note that genes are assigned to the cluster for which they had a maximal membership value.

#We can observe this data.

centers <- head(mf$centers)
kable(centers, caption = "Cluster centers", row.names = TRUE) %>%
  kable_styling(full_width = FALSE)

membership <- head(mf$membership)
kable(membership, caption = "Membership", row.names = TRUE) %>%
  kable_styling(full_width = FALSE)

sizes <- mf$size
names(sizes) <- 1:cluster

size <- head(sizes,10)
kable(t(size), caption = "Cluster size", row.names = F) %>%
  kable_styling(full_width = FALSE)

hcluster <- head(mf$cluster,10)
kable(t(hcluster), caption = "Cluster values", row.names = F) %>%
  kable_styling(full_width = FALSE)

#We save the data in a table for downstream analysis.
write.table( mf$centers, file=paste0("Center_mfuzz_rnaseq_c7_sponge.txt"), sep="\t")
write.table( mf$membership, file=paste0("Membership_mfuzz_rnaseq_c7_sponge.txt"), sep="\t")
write.table( mf$size, file=paste0("Size_mfuzz_rnaseq_c7_sponge.txt"), sep="\t")
write.table( mf$cluster, file=paste0("Cluster_mfuzz_rnaseq_c7_sponge.txt"), sep="\t")

### Profiling the cluster cores: acore function

#This function extracts genes forming the alpha cores of soft clusters. It produces an list of alpha cores including genes and their membership values for the corresponding cluster.

#Genes can be differentiated by examining whether they are included in a certain α- core. The increase of the the α-threshold led to a decreased number of genes included in the α-core. α-core of 0 is equivalent to hard clustering. increasing the α-core decrease the within-cluster variation.

#The use of the α-threshold can therefore act as a posteriori filtering of genes. 0.7 is the α-core suggested by the developers.

#-   cl = an object of class flcust as produced by mfuzz.

#-   min.acore = minimum membership values of gene belonging to the cluster core --\> α-core threshold


#mf <- mfuzz(data.s,c=cluster,m=m1)
mf
acore.list <- acore(data.s,cl=mf,min.acore=0.5)


#Save the data in a table (do the same for each cluster)

#The number in brackets is the number of clusters detected in the data
write.table(acore.list[[5]],"cluster5_BCT.txt",sep="\t")


## Session info

#We show the information of the session to know which software we have been using. This is a important steep to ensure reproducibility.

sessionInfo()
#save the info in a txt file
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#NOW FOR SUSY
rnaseq_counts <-read.table("Hpan_final.isoform.counts.matrix_SUSY_females_elements.txt", sep="\t", header = TRUE, row.names = 1)
head(rnaseq_counts, 5)
str(rnaseq_counts)

#We define some of the stages that we have in the analysis in a vector for easy reproducibility of the analysis.
#For females
stages <- c("Mar_OO","Apr_OO","May_lateOO", "May_EE", "May_ME", "Jun_Lar")


#For males
#stages <- c("Mar_NR","Mar_OO","Apr_OO","Apr_SP","May_OOE", "May_SP", "Jun_NR", "Jun_Lar", "Jun_FL")

### Prepare data

#We need transform the data to a matrix because mfuzz can not work with a data frame, and needs the `ExpressionSet` to work
matrix_rnaseq_counts <- as.matrix(rnaseq_counts)

#We save the data in a new object using `ExpressionSet`
normalizedSet <- ExpressionSet(assayData=matrix_rnaseq_counts)
normalizedSet
#When looking at the normalizedSet:
#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 179328 features, 6 samples 
#element names: exprs 
#protocolData: none
#phenoData: none
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation:  

#Mfuzz don't like NA or missing values, so we have to filter any possible NA values. We also standardize the expression values so that the average expression value for each gene is zero and the standard deviation of its expression profile is one. This step ensures that vectors of genes with similar changes in expression are close in Euclidean space
data.r <- filter.NA(normalizedSet, thres=0.5) #based on NA missing values
data.k <- fill.NA(data.r,mode="knn") #replace remaining missing values by the average values expression value of the corresponding gene. k-nearest neighbour method can be used (mode=`knn'/`wknn').

#Note that you won't get NA values when using edgeR, those are DESeq2 outputs. However, you may want to remove genes which are expressed at low levels or show only small changes in expression. For that include the filtering step: filter.std --\> exclude genes with low standard deviation
data.f <- filter.std(data.k, min.std=0.001, visu=TRUE) #based on standard deviation. 
#If visu is set to TRUE, the ordered standard deviations of genes' expression values will be plotted
#Here the result is that 131511 genes are excluded

#Finally standarize to make transcripts (or genes or proteins) comparable to Mfuzz. Note that the data.s is an object of the class `ExpressionSet`
data.s <- standardise(data.f)
data.s

#ExpressionSet (storageMode: lockedEnvironment)
#assayData: 196650 features, 6 samples 
#element names: exprs 
#protocolData: none
#phenoData: none
#featureData: none
#experimentData: use 'experimentData(object)'
#Annotation: 

#### Setting the parameters for FCM clustering

#Usually we will optimize the parameters for the MFuzz analysis for each dataset. This process is memory demanding and for the sake of time we will not be ruining it in this practical.
#The parameters used in Mfuzz are:

# -   **fuzzifier value (m)**: We will determine the optimal setting of fuzzifie value using `mestimate`. This will be only an estimation since the final result tuned depending on the performance of the function.

#As this is a subclustering it can be soft if one gene can be in more than one cluster, or hard if you limit the gene to be only in one cluster. This is determined by the fuzzifier value (m). Under a biological perspective, this soft clustering is more appropiate as one gene can be involved in different functions.

#-   **number of clusters (cluster)**: final number of clusters. This will be determined using repeated soft clustering for a range of cluster numbers. This will report the minimum centroid distance. The minimum centroid distance can be used as cluster validity index. For an optimal cluster number, we may see a 'drop' of minimum centroid distance wh plotted versus a range of cluster number and a slower decrease of the minimum centroid distance for higher cluster number.

#For parameter selection use your `ExpressionSet` object, i.e., data.s:
m1 <- mestimate(data.s)
m1
#m1=1.694431 (for ALL_MTX_female)
#m1=1.723651 (for BCT)
#For election of c, either cselection or:
tmp <- Dmin(data.s,m=m1,crange=seq(4,40,4),repeats=3,visu=TRUE) # Note: This calculation might take some time

#The suitable number of clusters might be those that show the decrease for c

#The minimum centroid distance is defined as the minimum distance between two cluster centers produced by the c-means clusterings. The average minimum centroid distance for the given range of cluster number is returned.

#The minimum centroid distance can be used as cluster validity index. For an optimal cluster number, we may see a 'drop' of minimum centroid distance wh plotted versus a range of cluster number and a slower decrease of the minimum centroid distance for higher cluster number. Thus, the optimal number of clusters is often chosen at the point where the minimum centroid distance starts to level off, and adding more clusters doesn't provide much improvement in terms of separation

m1=1.700821 #m1=1.694431 (from mestimate)
cluster=15


### Run MFuzz

#We perform mfuzz clustering and represent the plot.

#Note: every time you rerun the mfuzz the order of the clusters may change (but not the structure)
mf <- mfuzz(data.s,c=cluster,m=m1)
mfuzz.plot(data.s,cl=mf, time.labels=stages,new.window=FALSE)

#HE HECHO HASTA AQUÍ: MAÑANA PROBAR QUITANDO FL
#Colours in the plot represent the belonginess ("membership") to the pattern, e.e., how good or bad each gene (line) fits to the cluster pattern. Green or yellow means they belong less to the clustering. Blue and red colored genes better fit to the cluster dynamics.

#> As you modify the recommended parameters the red parts in the clustering is lost, so we loose reliability.
# Note: there is always important to have this red genes that indicate the belonginess to the cluster is strong, is good. Thus, it is important to find a balance between m and cluster:
# As it is a soft clustering you always expect some deviation so you need a balance between being so soft (same gene in many clusters) or hard (one gene only in one cluster).
# Also with the number of clusters you need to reach a compromise: less clusters will show a lot of genes assigned to one cluster but with less belonginess whereas many clusters will get genes with low belonginess to the cluster.

#You can try to modify the m and cluster parameters for FCM clustering and how mfuzz clustering change
#For instance:
#m2=1.6
#cluster2=30

#mf2 <- mfuzz(data.s,c=cluster2,m=m2)
#mfuzz.plot(data.s,cl=mf2, time.labels=stages,new.window=FALSE)
```

```{r}
m3=1.5
cluster3=30

mf3 <- mfuzz(data.s,c=cluster3,m=m3)
mfuzz.plot(data.s,cl=mf3, time.labels=stages,new.window=FALSE)
```

```{r}
m4=1.4
cluster4=30

mf4 <- mfuzz(data.s,c=cluster4,m=m4)
mfuzz.plot(data.s,cl=mf4, time.labels=stages,new.window=FALSE)
```

```{r}
m5=1.3
cluster5=30

mf5 <- mfuzz(data.s,c=cluster5,m=m5)
mfuzz.plot(data.s,cl=mf5, time.labels=stages,new.window=FALSE)
```

#### **Save the generated clusters in a PDF**

#NB: I would add in the file name the parameters for FCM clustering

pdf(paste0("mfuzz_rnaseq_m1.7_c15_susy.pdf"), height=10, width=15)
mfuzz.plot(data.s,cl=mf, time.labels=stages,new.window=FALSE)
dev.off()
```

#### **Extract the genes in each expression cluster**

#The results are stored in the `mf` object. Some interesting metrics are: - centers: the final cluster centers. - size: the number of genes in each cluster of the closest hard clustering. - cluster: indices of the clusters where the genes are assigned to for the closest hard clustering, as obtained by assigning points to the (first) class with maximal membership. - membership: a matrix with the membership values of the data points to the clusters.

#Note that genes are assigned to the cluster for which they had a maximal membership value.

#We can observe this data.

centers <- head(mf$centers)
kable(centers, caption = "Cluster centers", row.names = TRUE) %>%
  kable_styling(full_width = FALSE)

membership <- head(mf$membership)
kable(membership, caption = "Membership", row.names = TRUE) %>%
  kable_styling(full_width = FALSE)

sizes <- mf$size
names(sizes) <- 1:cluster

size <- head(sizes,10)
kable(t(size), caption = "Cluster size", row.names = F) %>%
  kable_styling(full_width = FALSE)

hcluster <- head(mf$cluster,10)
kable(t(hcluster), caption = "Cluster values", row.names = F) %>%
  kable_styling(full_width = FALSE)

#We save the data in a table for downstream analysis.
write.table( mf$centers, file=paste0("Center_mfuzz_rnaseq_c15_susy.txt"), sep="\t")
write.table( mf$membership, file=paste0("Membership_mfuzz_rnaseq_c15_susy.txt"), sep="\t")
write.table( mf$size, file=paste0("Size_mfuzz_rnaseq_c15_susy.txt"), sep="\t")
write.table( mf$cluster, file=paste0("Cluster_mfuzz_rnaseq_c15_susy.txt"), sep="\t")

### Profiling the cluster cores: acore function

#This function extracts genes forming the alpha cores of soft clusters. It produces an list of alpha cores including genes and their membership values for the corresponding cluster.

#Genes can be differentiated by examining whether they are included in a certain α- core. The increase of the the α-threshold led to a decreased number of genes included in the α-core. α-core of 0 is equivalent to hard clustering. increasing the α-core decrease the within-cluster variation.

#The use of the α-threshold can therefore act as a posteriori filtering of genes. 0.7 is the α-core suggested by the developers.

#-   cl = an object of class flcust as produced by mfuzz.

#-   min.acore = minimum membership values of gene belonging to the cluster core --\> α-core threshold


#mf <- mfuzz(data.s,c=cluster,m=m1)
mf
acore.list <- acore(data.s,cl=mf,min.acore=0.5)


#Save the data in a table (do the same for each cluster)

#The number in brackets is the number of clusters detected in the data
write.table(acore.list[[5]],"cluster5_BCT.txt",sep="\t")


## Session info

#We show the information of the session to know which software we have been using. This is a important steep to ensure reproducibility.

sessionInfo()
#save the info in a txt file
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
