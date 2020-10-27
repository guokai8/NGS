
# RNA-centric Analysis of CITE-seq data with ADT, HTO and RNA libraries, with Seurat v4.0 (This code works for Seurat v3, too)


### DataSets
- Human PBMC 
- Hashed with 2 hashtags HTO1 and 2, 
- Stained with TotalSeq-A CD3, CD4 and CD8 antibodies 
- ADT, HTO and RNA libraries generated using 10x Genomics Chromium 3' gene expression v3 kits. 
- Sequenced on 2 lanes of a Illumina Sequencer 

### Softwares
You will need the following softwares to process the fastq files
- [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation) (v3.0 or above)
-- reference files should be downloaded prior to analysis 
- [Seurat v4](https://satijalab.org/seurat/install.html) (same code can work on v3.0

### Folder layout
The folder layout of this guide is
```
/path/to/data : the main project folder
/path/to/data/fastq : the fastq files 
/path/to/data/ref: the 10x Genomics reference 
/path/to/data/output     
``` 

### Fastq Files
A typical bcl2fastq run will produce separate fastq files for the ADT, HTO and RNA fractions. the fastq files will have names like "ADT_S1_L001_R1_001.fastq.gz" ( in a format of \<SampleID\>\_\<Name\>\_\<Lane\>\_\<Read\>_001.fastq.gz ) 
Please also take a look at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis for more information about naming of the fastq files

You should get a of fastqs this guide assume:
ADT:
> ADT_S1_L001_R1_001.fastq.gz, ADT_S1_L001_R2_001.fastq.gz, ADT_S1_L002_R1_001.fastq.gz, ADT_S1_L002_R2_001.fastq.gz

HTO:
> HTO_S2_L001_R1_001.fastq.gz, HTO_S2_L001_R2_001.fastq.gz, HTO_S2_L002_R1_001.fastq.gz, HTO_S2_L002_R2_001.fastq.gz

RNA:
>RNA_S3_L001_R1_001.fastq.gz, RNA_S3_L001_R2_001.fastq.gz, RNA_S2_L002_R1_001.fastq.gz, RNA_S3_L002_R2_001.fastq.gz

It's important to follow the naming convention for Cell Ranger to find the fastq files.


## Count UMI per cell with Cell Ranger

### Prepare CSV files for Cell Ranger
There are 2 required csv files for Cell Ranger, to describe fastq files (Libraries.csv) and antibodies (Features.csv) repectively. 
1. Libraries.csv 
```
fastqs,sample,library_type
/path/to/data/fastq,ADT,Antibody Capture
/path/to/data/fastq,HTO,Antibody Capture
/path/to/data/fastq,RNA,Gene Expression
```

Note that the "sample" column should match the SampleID of fastq file name( i.e. for a set of ADT_S1_L001_R1_001.fastq.gz, should specify "ADT" in the "sample" column) 

2. Features.csv

Here's a list that contain 3 TotalSeq-A antibodies and 2 Hashtags 
```
id,name,read,pattern,sequence,feature_type
ADT_A0034,A0034_hCD3,R2,5P(BC),CTCATTGTAACTCCT,Antibody Capture
ADT_A0072,A0072_hCD4,R2,5P(BC),TGTTCCCGCTCAACT,Antibody Capture
ADT_A0046,A0046_hCD8,R2,5P(BC),GCGCAACTTGATGAT,Antibody Capture
HTO_A0251,A0251_Hash1,R2,5P(BC),GTCAACTCTTTAGCG,Antibody Capture
HTO_A0252,A0252_Hash2,R2,5P(BC),TGATGGCCTATTGGG,Antibody Capture
```

If using TotalSeq-B or C antibodies, the "pattern" column should contains the leading Ns, for example
```
id,name,read,pattern,sequence,feature_type
ADT_B0034,B0034_hCD3,R2,5PNNNNNNNNNN(BC),CTCATTGTAACTCCT,Antibody Capture
```
For the ease of executing the Cell Ranger, please put these csv files at the main project folder `/path/to/data` 

### Run Cell Ranger
Go to the project folder
```
$ cd /path/to/data/output
```

Run the Cell Ranger with required files, given that we run a **chromium 3' v3 gene expression** lane with **10000** cells.
```
$ cellranger count --id=PBMC \
  --transcriptome=/path/to/data/ref/GRCh38 \
  --libraries=/path/to/data/Libraries.csv \
  --feature-ref=/path/to/data/Features.csv \
  --expect-cells=10000 --chemistry=SC3Pv3
```

This will create a output folder in `/path/to/data/output/PBMC`, you can find the mex format of UMI count in `/path/to/data/output/PBMC/outs/filtered_feature_bc_matrix`

## Demultiplexing using Seurat

### Load data into R
**The following commands are run in R**

Set the working directory into main project folder
```
> setwd('/path/to/data/')
```

Load the required library
```
> library(Seurat)
> library(Matrix)
```

Load the mex format output into R
```
> matrix_dir = "output/PBMC/outs/filtered_feature_bc_matrix/"

> barcode.path = paste0(matrix_dir, "barcodes.tsv.gz")
> features.path = paste0(matrix_dir, "features.tsv.gz")
> matrix.path = paste0(matrix_dir, "matrix.mtx.gz")
> mat = readMM(file = matrix.path)
> feature.names = read.delim(features.path, 
                            header = FALSE,
                           stringsAsFactors = FALSE)
> barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
> colnames(mat) = substr(barcode.names$V1, 1, 16)
> rownames(mat) = feature.names$V1
```
Find the ADT and HTO features
```
> ADT = grepl("ADT_", rownames(mat))
> HTO = grepl("HTO_", rownames(mat))
```
### Create Seurat object

Create a Seurat object, using RNA expression only
```
> pbmc = CreateSeuratObject(mat[!(ADT | HTO),])
```

Normalize and scale the data
```
> pbmc = NormalizeData(pbmc)
> pbmc = FindVariableFeatures(pbmc)
> pbmc = ScaleData(pbmc)
```

Load the ADT and HTO data 
```
> pbmc[["ADT"]] = CreateAssayObject(counts = mat[feat_ADT,])
> pbmc = NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
> pbmc = ScaleData(pbmc, assay = "ADT")

> pbmc[["HTO"]] = CreateAssayObject(counts = mat[feat_HTO,])
> pbmc = NormalizeData(pbmc, assay = "HTO", normalization.method = "CLR")
> pbmc = ScaleData(pbmc, assay = "HTO")

```

Demultiplexing with Seurat
```
> pbmc = HTODemux(pbmc, assay = "HTO", positive.quantile =  0.99)
```

Look at the heatmap of HTO
```
> HTOHeatmap(pbmc, assay = "HTO", ncells = 5000)
```

Please note that the quantile threshold should be modified according to your dataset. Values from 0.99 to 0.80 are common, while some even go lower than 0.6
We recommend to look at the number of result Negative vs. Doublets. As the threshold going down, you will have less Negatives but more Doublets. The best threshold should be where you first get more Doublets than Negatives.

For more information, please visit [HTODemux](https://satijalab.org/seurat/hashing_vignette.html) page of Seurat package  

## Analyzing data with Seurat
```
> table(pbmc$HTO_classification.global)
```
## Visualize enrichment for selected HTOs with ridge plots
```
Idents(pbmc) <- "HTO_maxID"
RidgePlot(pbmc, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:2], ncol = 2)
```

## Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
```
FeatureScatter(pbmc.hashtag, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
```
## Generate a two dimensional tSNE embedding for HTOs.Here we are grouping cells by singlets and doublets for simplicity.

```
pbmc.subset <- subset(pbmc, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO

hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
pbmc.subset <- RunTSNE(pbmc.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(pbmc.subset)
```
## Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.
```
pbmc.singlet <- subset(pbmc, idents = "Singlet")

# Select the top 1000 most variable features
pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

# Run PCA
pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

```
