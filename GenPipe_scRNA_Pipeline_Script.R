
## ..[] Write a single step GenPipes pipeline that launches this R script

##----install packages & load libraries
#install.packages('Seurat')

##load libraries 
library('Seurat')
library('dplyr')

##import data using Read10x function 
#Data: 1k Brain cells from E18 Mouse, obtained from: https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.0.1/neurons_900
neural_E18_data <- Read10X(data.dir = "../GenPipe_scRNA_Pipeline/data")
neural_E18_data

#create seurat object 
neural_E18 <- CreateSeuratObject(counts = neural_E18_data, project = "Neural_E18_1k")
neural_E18 

##filter cells with more than 20% mitochondrial genes 
#calculate mitochondrial QC using PercentageFeatureSet function & output to column named mQC_Percent
neural_E18[["mQC_Percent"]] <- PercentageFeatureSet(object = neural_E18, pattern = "^MT-")
#subset to cells with mQC_Percent > 20%
neural_E18 <- subset(x = neural_18, subset = mQC_Percent > 20)
neural_E18

##run a PCA on the most variable genes 
#identify top 10 most variable genes using head
neural_E18_variable <- head(x = VariableFeatures(object = neural_E18), 10)
neural_E18_variable

#prep for PCA using ScaleData function 
all_genes <- rownames(x = neural_E18)
neural_E18 <- ScaleData(object = neural_E18, features = all_genes)

#run PCA analysis
neural_E18 <- RunPCA(object = neural_E18, features = VariableFeatures(object = neural_E18))
#visualize 
pca.plot(neural_E18_variable, 1, 2, pt.size = 2)

##run a tSNE analysis 
neural_E18 <- RuntSNE(object = neural_E18, dims = 1:10)
#plot 
PlottSNE(neural_E18, pt.size = 1)

#save to pdf 
save(neural_E18, file = "~/GenPipe_scRNA_Pipeline/neural_E18_tSNE_plot.pdf")
