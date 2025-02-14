# 2 sample analysis- all cells



# 1. Install package, load libraries and data
# 2.

## If used please cite Luskin and Li et al peri-LC paper

# 1. preprocess (including removing doublet cells)

library(Seurat)
library(dplyr)
library(magrittr)
library(grid)
IRdisplay::display_html("<style> .container { width:95% !important; } </style>")
library("ggplot2")
library("xlsx")
library(tidyr)
library(ggplot2)
library(cowplot)
library(svglite)
library(clustree)
library(scales)
library(reshape2)
library(metap)
#Set working directory
setwd('/FILL_IN_DIR/')
print("Library Loaded.")

SAMPLE_A.data <-Read10X(data.dir = "/FILL_IN_DIR_A/raw_feature_bc_matrix")
colnames(SAMPLE_A.data) = paste0(colnames(SAMPLE_A.data),"SAMPLE_A")
SAMPLE_A<- CreateSeuratObject(counts = SAMPLE_A.data, min.cells = 3, min.features = 200, project = "SAMPLE_A")
new_SAMPLE_A_barcodes <- gsub("-1","",colnames(SAMPLE_A))
SAMPLE_A <- RenameCells(SAMPLE_A, new.names = new_SAMPLE_A_barcodes)
SAMPLE_A@meta.data$stim <- "SAMPLE_A"

SAMPLE_A_CLEAN<-rownames(read.table(file="/FILL_IN_DIR_A/Final_nondoublets_groups_cntl.txt"))
SAMPLE_A_DOUBLET<-rownames(read.table(file="/FILL_IN_DIR_A/Final_doublets_groups_cntl.txt"))
Doubletrate_SAMPLE_A<-100*length(SAMPLE_A_DOUBLET)/(length(SAMPLE_A_CLEAN)+length(SAMPLE_A_DOUBLET)) #calculate doublet percentage
Doubletrate_SAMPLE_A
SAMPLE_A_CLEAN_PROJECTNAME <- gsub(".1Seurat","SAMPLE_A",(SAMPLE_A_CLEAN))
SAMPLE_A<-subset(x=SAMPLE_A,cells=SAMPLE_A_CLEAN_PROJECTNAME)
length(SAMPLE_A_CLEAN_PROJECTNAME)

SAMPLE_B.data <-Read10X(data.dir = "/FILL_IN_DIR_B/raw_feature_bc_matrix")
colnames(SAMPLE_B.data) = paste0(colnames(SAMPLE_B.data),"SAMPLE_B")
SAMPLE_B<- CreateSeuratObject(counts = SAMPLE_B.data, min.cells = 3, min.features = 200, project = "SAMPLE_B")
new_SAMPLE_B_barcodes <- gsub("-1","",colnames(SAMPLE_B))
SAMPLE_B <- RenameCells(SAMPLE_B, new.names = new_SAMPLE_B_barcodes)
SAMPLE_B@meta.data$stim <- "SAMPLE_B"

SAMPLE_B_CLEAN<-rownames(read.table(file="/FILL_IN_DIR_B/Final_nondoublets_groups_cntl.txt"))
SAMPLE_B_DOUBLET<-rownames(read.table(file="/FILL_IN_DIR_B/Final_doublets_groups_cntl.txt"))
Doubletrate_SAMPLE_B<-100*length(SAMPLE_B_DOUBLET)/(length(SAMPLE_B_CLEAN)+length(SAMPLE_B_DOUBLET)) #calculate doublet percentage
Doubletrate_SAMPLE_B
SAMPLE_B_CLEAN_PROJECTNAME <- gsub(".1Seurat","SAMPLE_B",(SAMPLE_B_CLEAN))
SAMPLE_B<-subset(x=SAMPLE_B,cells=SAMPLE_B_CLEAN_PROJECTNAME)
length(SAMPLE_B_CLEAN_PROJECTNAME)

SAMPLE_C.data<-Read10X(data.dir = "/FILL_IN_DIR_C/raw_feature_bc_matrix")
colnames(SAMPLE_C.data) = paste0(colnames(SAMPLE_C.data),"SAMPLE_C")
SAMPLE_C<-CreateSeuratObject(counts = SAMPLE_C.data, min.cells = 3, min.features = 200, project = "SAMPLE_C")
new_SAMPLE_C_barcodes <- gsub("-1","",colnames(SAMPLE_C))
SAMPLE_C <- RenameCells(SAMPLE_C, new.names = new_SAMPLE_C_barcodes)
SAMPLE_C@meta.data$stim <- "SAMPLE_C"

SAMPLE_C_CLEAN<-rownames(read.table(file="/FILL_IN_DIR_C/Final_nondoublets_groups_cntl.txt"))
SAMPLE_C_DOUBLET<-rownames(read.table(file="/FILL_IN_DIR_C/Final_doublets_groups_cntl.txt"))
Doubletrate_SAMPLE_C<-100*length(SAMPLE_C_DOUBLET)/(length(SAMPLE_C_CLEAN)+length(SAMPLE_C_DOUBLET)) #calculate doublet percentage
Doubletrate_SAMPLE_C
SAMPLE_C_CLEAN_PROJECTNAME <- gsub(".1Seurat","SAMPLE_C",(SAMPLE_C_CLEAN))
SAMPLE_C<-subset(x=SAMPLE_C,cells=SAMPLE_C_CLEAN_PROJECTNAME)
length(SAMPLE_C_CLEAN_PROJECTNAME)

mito.features <- grep(pattern = "^mt-", x = rownames(x =SAMPLE_A), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = SAMPLE_A, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = SAMPLE_A, slot = 'counts'))
SAMPLE_A[['percent.mito']] <- percent.mito

mito.features <- grep(pattern = "^mt-", x = rownames(x =SAMPLE_B), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = SAMPLE_B, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = SAMPLE_B, slot = 'counts'))
SAMPLE_B[['percent.mito']] <- percent.mito

mito.features <- grep(pattern = "^mt-", x = rownames(x =SAMPLE_C), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = SAMPLE_C, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = SAMPLE_C, slot = 'counts'))
SAMPLE_C[['percent.mito']] <- percent.mito

# you do not need to run this, just make sure to note what cut off you use
SAMPLE_A <- subset(x = SAMPLE_A, subset = nCount_RNA > 700 & nCount_RNA < 15000 & percent.mito < 0.1)
SAMPLE_B <- subset(x = SAMPLE_B, subset = nCount_RNA > 700 & nCount_RNA < 15000 & percent.mito < 0.1)
SAMPLE_C<- subset(x = SAMPLE_C, subset = nCount_RNA > 700 & nCount_RNA < 15000 & percent.mito < 0.1)



# Integration from here

SAMPLE_A<- NormalizeData(object = SAMPLE_A,verbose = FALSE) 
SAMPLE_B<- NormalizeData(object = SAMPLE_B,verbose = FALSE) 
SAMPLE_C<- NormalizeData(object = SAMPLE_C,verbose = FALSE) 

SAMPLE_A<- FindVariableFeatures(object =SAMPLE_A,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = SAMPLE_A))
SAMPLE_B<- FindVariableFeatures(object =SAMPLE_B,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object =SAMPLE_B))
SAMPLE_C<- FindVariableFeatures(object =SAMPLE_C,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
length(x = VariableFeatures(object = SAMPLE_C))

PROJECTNAME.list<-objects()
PROJECTNAME.list$SAMPLE_A<-SAMPLE_A
PROJECTNAME.list$SAMPLE_B<-SAMPLE_B
PROJECTNAME.list$SAMPLE_C<-SAMPLE_C

# 2. Integration of multiple dataset

# Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input. Here, we integrate three of the objects into a reference 
# 
# We use all default parameters here for identifying anchors, including the ‘dimensionality’ of the dataset (30; feel free to try varying this parameter over a broad range, for example between 10 and 50).

reference.list <- PROJECTNAME.list[c("SAMPLE_A","SAMPLE_B","SAMPLE_C")]
PROJECTNAME.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) # this dims is the which CCAs you use )

PROJECTNAME.integrated <- IntegrateData(anchorset = PROJECTNAME.anchors, dims = 1:30) # in this case the dims is the nmber of PCs to use for weight (correction vector)

# After running IntegrateData, the Seurat object will contain a new Assay with the integrated expression matrix. Note that the original (uncorrected values) are still stored in the object in the “RNA” assay, so you can switch back and forth.
# 
# We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP. The integrated datasets cluster by cell type, instead of by technology.

DefaultAssay(object = PROJECTNAME.integrated) <- "integrated" # you can switch between RNA and integrated any time

PROJECTNAME.integrated <- ScaleData(object = PROJECTNAME.integrated, , vars.to.regress = c("nCount_RNA", "percent.mito"),verbose = FALSE)

PROJECTNAME.integrated <- RunPCA(object = PROJECTNAME.integrated, npcs = 30, verbose = FALSE)# npcs 50 default

PROJECTNAME.integrated <- RunUMAP(object = PROJECTNAME.integrated, reduction = "pca", dims = 1:30)

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23", "24") #as many clusters as you have - change as necessary every time this comes up
names(new.cluster.ids) <- levels(PROJECTNAME.integrated)
PROJECTNAME.integrated <- RenameIdents(PROJECTNAME.integrated, new.cluster.ids)
#DimPlot(PROJECTNAME.integrated, reduction="umap", label=TRUE, pt.size=0.5)+NoLegend()
p1 <- DimPlot(object = PROJECTNAME.integrated, reduction = "umap", group.by = "stim")
plot_grid(p1)
#ggsave(file="/FILL_IN_DIR/UMAP_FILEPROJECTNAME.pdf",plot=image,width=10,height=10)
#ggsave(file="UMAP_FILEPROJECTNAME.svg",plot=image,width=6,height=6)

# integrated analysis

str(PROJECTNAME.integrated)

PROJECTNAME.integrated <- FindNeighbors(object = PROJECTNAME.integrated, dims=1:30)

PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.3, dims.use = 1:30, graph.name = "integrated_snn")

#help(FindClusters)

PROJECTNAME.integrated<-readRDS(file="DATE_PROJECTNAME_integrated.rds")
DefaultAssay(PROJECTNAME.integrated) <- "integrated"

PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.1, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.2, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.3, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.4, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.5, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.6, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.7, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.8, dims = 1:30)
PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.9, dims = 1:30)

clustree(PROJECTNAME.integrated@meta.data, prefix = "integrated_snn_res.", legend.position='bottom', edge_width = 0.5, node_size_range=c(3,10),node_text_size = 2, layout='tree') +
  guides(edge_colour = FALSE, edge_alpha = FALSE, fill=guide_legend(title="resolution")) +
  theme(plot.title=element_text(size=12))+
  ggtitle('Resolution analysis')
ggsave(file="All_cells_resolution_clustree.svg",width=4,height=5)
ggsave(file="All_cells_resolution_clustree.png",width=5,height=5)

## these plots will help you figure out what resolution to use; you should also consider using something like ChooseR to figure out the optimal resolution

DefaultAssay(PROJECTNAME.integrated) <- "integrated"

PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.3, dims = 1:30)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23", "24")
names(new.cluster.ids) <- levels(PROJECTNAME.integrated)
PROJECTNAME.integrated <- RenameIdents(PROJECTNAME.integrated, new.cluster.ids)
DimPlot(object = PROJECTNAME.integrated, reduction = "umap", group.by = "integrated_snn_res.0.3", label = TRUE, repel = TRUE) 
#ggsave(file="/FILL_IN_DIR/All_cells/allcell_umap_number_2000_all20_int.pdf",width=10,height=10)

PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.4, dims = 1:30)
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23", "24", "25","26","27","28")
names(new.cluster.ids) <- levels(PROJECTNAME.integrated)
PROJECTNAME.integrated <- RenameIdents(PROJECTNAME.integrated, new.cluster.ids)
DimPlot(object = PROJECTNAME.integrated, reduction = "umap", group.by = "integrated_snn_res.0.4", label = TRUE, repel = TRUE) 
#ggsave(file="/FILL_IN_DIR/All_cells/allcell_umap_number_2000_all20_int.pdf",width=10,height=10)

DefaultAssay(PROJECTNAME.integrated) <- "RNA"

PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.3, dims = 1:30)
DimPlot(object = PROJECTNAME.integrated, reduction = "umap", group.by = "integrated_snn_res.0.3", label = TRUE, repel = TRUE) 
#ggsave(file="All_umap_features_2000_res0_3.svg",width=10,height=10)

#?DimPlot

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23", "24")
names(new.cluster.ids) <- levels(PROJECTNAME.integrated)
PROJECTNAME.integrated <- RenameIdents(PROJECTNAME.integrated, new.cluster.ids)
DimPlot(object = PROJECTNAME.integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.65) +
  theme_classic() + 
  NoLegend() + 
  theme(panel.background = element_rect(colour = "black", size=1, fill=NA),
        axis.line=element_blank())
#ggsave(file="DATE_All_cells_names.svg",width=4,height=4)
#ggsave(file="/FILL_IN_DIR/DATE_All_cells_names.png",width=4,height=4)

PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.3, dims = 1:30)
DimPlot(object = PROJECTNAME.integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size=2.5) +
  theme_classic() + 
  NoLegend() + 
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank())
#ggsave(file="DATE_All_cells_labeled_res0_3.svg",width=3,height=3, units="in")

DimPlot(object = PROJECTNAME.integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.5, label.size=2.5) +
  theme_classic() + 
  NoLegend() + 
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank())
#ggsave(file="DATE_All_cells_unlabeled_res0_3.svg",width=3,height=3, units="in")


PROJECTNAME.integrated <- FindClusters(PROJECTNAME.integrated, resolution = 0.4, dims = 1:30)
DimPlot(object = PROJECTNAME.integrated, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5, label.size=2.5) +
  theme_classic() + 
  NoLegend() + 
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank())
#ggsave(file="DATE_All_cells_labeled_res0_4.svg",width=3,height=3, units="in")

DimPlot(object = PROJECTNAME.integrated, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.5, label.size=2.5) +
  theme_classic() + 
  NoLegend() + 
  theme(axis.line=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank())
#ggsave(file="DATE_All_cells_unlabeled_res0_4.svg",width=3,height=3, units="in")


saveRDS(PROJECTNAME.integrated, file = "DATE_PROJECTNAME_integrated.rds")

saveRDS(SAMPLE_A, file = "DATE_SAMPLE_A.rds")
saveRDS(SAMPLE_B, file = "DATE_SAMPLE_B.rds")
saveRDS(SAMPLE_C, file = "DATE_SAMPLE_C.rds")





PROJECTNAME.integrated<-readRDS(file="DATE_PROJECTNAME_integrated.rds")

DefaultAssay(PROJECTNAME.integrated) <- "RNA"

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23", "24")
names(new.cluster.ids) <- levels(PROJECTNAME.integrated)
PROJECTNAME.integrated <- RenameIdents(PROJECTNAME.integrated, new.cluster.ids)
table(Idents(PROJECTNAME.integrated))

F<-FeaturePlot(object =PROJECTNAME.integrated, features = c("Stmn2","Thy1","Aldoc","Aqp4","Tmem119","C1qc","Mog","Opalin","Pdgfra","Gpr17","Cldn5","Tagln","Vtn","Flt1","Spp1","Foxj1","Fam216b","Cspg4"))
F
#ggsave(file="DATE_PROJECTNAME_int_feature_all_genes.pdf",width=10,height=10)

F<-FeaturePlot(object =PROJECTNAME.integrated, features = c("Stmn2","Thy1","Dbh","Th","Chat","Slc6a4","Slc6a3","Slc17a7","Slc17a6","Slc17a8","Slc32a1","Gad1","Gad2"))
F
#ggsave(file="DATE_PROJECTNAME_int_feature_key_neuronal_genes.pdf",width=10,height=10)

F<-FeaturePlot(object =PROJECTNAME.integrated, features = c("Fos","Egr1","Arc","Jun","Fosb"))
#ggsave(file="DATE_PROJECTNAME_int_feature_IEG_genes.pdf",width=10,height=10)


image=VlnPlot(PROJECTNAME.integrated, features = c("Stmn2","Thy1","Aldoc","Aqp4","Tmem119","C1qc","Mog","Opalin","Pdgfra","Gpr17","Cldn5","Tagln","Vtn","Flt1","Spp1","Foxj1","Fam216b","Cspg4"), pt.size = 0)
ggsave(file="All_cells_violin_cell_subtypes.svg", plot=image, width=14, height=14, units="in")
ggsave(file="All_cells_violin_cell_subtypes.png", plot=image, width=14, height=14, units="in")

image=VlnPlot(PROJECTNAME.integrated, features = c("Dbh","Th", "Chat", "Slc6a4", "Slc6a3", "Slc17a7", "Slc17a6", "Slc17a8","Slc32a1","Gad2","Gad1"), pt.size = 0)
ggsave(file="All_cells_violin_neuron_subtypes.svg", plot=image, width=14, height=14, units="in")
ggsave(file="All_cells_violin_neuron_subtypes.png", plot=image, width=14, height=14, units="in")



new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20", "21", "22","23", "24")
names(new.cluster.ids) <- levels(PROJECTNAME.integrated)
PROJECTNAME.integrated <- RenameIdents(PROJECTNAME.integrated, new.cluster.ids)
new.ident=new.cluster.ids
for (i in 1:length(new.ident)){
  assign(paste("barcode_",new.ident[i],sep=""),colnames(PROJECTNAME.integrated@assays$RNA@data[,which(Idents(object=PROJECTNAME.integrated) %in% new.ident[i])]))# this gives all barcodes in cluster
}

#go through the clusters and based on the marker genes expressed, categorize the clusters into cell identities

Neuron_barcode<-c(barcode_2,barcode_3,barcode_6,barcode_10,barcode_11,barcode_12,barcode_13,barcode_14,barcode_15,barcode_16,barcode_17,barcode_18,barcode_19,barcode_21,barcode_23)
#NSC_barcode<-c()
Astrocyte_barcode<-c(barcode_5)
Oligo_barcode<-c(barcode_1,barcode_4)
Microglia_barcode<-c(barcode_7)
OPC_barcode<-c(barcode_8,barcode_24)
#Mural_barcode<-c()
Endothelial_barcode<-c(barcode_22)
Ependymal_barcode<-c()
Ambiguous_barcode<-c(barcode_9,barcode_20)

type<-numeric()
for (i in 1:dim(PROJECTNAME.integrated@meta.data)[1]){
  if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Neuron_barcode){type[i]<-"Neurons"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Astrocyte_barcode){type[i]<-"Astrocytes"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Oligo_barcode){type[i]<-"Oligodendrocytes"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Microglia_barcode){type[i]<-"Microglia"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% OPC_barcode){type[i]<-"OPCs"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Endothelial_barcode){type[i]<-"Endothelial cells"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Ependymal_barcode){type[i]<-"Ependymal cells"}
  else if(rownames(PROJECTNAME.integrated@meta.data)[i] %in% Ambiguous_barcode){type[i]<-"Ambiguous"}
  else{type[i]<-"NA"}
}
PROJECTNAME.integrated@meta.data$type<-type
PROJECTNAME.integrated@meta.data$type<-factor(PROJECTNAME.integrated@meta.data$type,levels=c("Neurons","Astrocytes","Oligodendrocytes","OPCs","Microglia", "Endothelial cells", "Ependymal cells", "Ambiguous"))
colors<-c("Astrocytes"="#8a8a8a","Oligodendrocytes"="#8a8a8a","Neurons"="#fc9403","OPCs"="#8a8a8a","Microglia"="#8a8a8a","Endothelial cells"="#8a8a8a", "Ependymal cells"="#8a8a8a", "Ambiguous"="#8a8a8a")
subset_meta<-subset(PROJECTNAME.integrated@meta.data,type=="Neurons"|type=="Astrocytes"|type=="Oligodendrocytes"|type=="OPCs"|type=="Microglia"|type=="Endothelial cells"|type=="Ependymal cells"|type=="Ambiguous")


ggplot(subset_meta,aes_string(x="type",y="nCount_RNA",fill="type"))+
  geom_violin(scale = "width", color="#8a8a8a")+scale_fill_manual(values=colors)+
  stat_summary(fun=median, geom="point", size=1, color="white")+ylab("UMIs")+ 
  theme_classic()+
  theme(axis.title.x=element_blank(),
        ,axis.text.y=element_text(size=15, angle=90),axis.title.y=element_text(size=20,angle=90,face="bold",margin = margin(t = 10, r = 10, b = 0, l = 10),vjust=0.5),axis.text.x=element_text(size=20, face="bold", angle=90, vjust=0.5)
        ,axis.title=element_text(size=20,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.5),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")) 
#ggsave(file="FILL_IN_DIR/DATE/All_cells/UMI_celltype_int.pdf",height=4, width=4 , paper = "letter")
#ggsave(file="DATE_PROJECTNAME_integrated_QI_UMI_cell_types.png", height=4, width=4)
ggsave(file="All_cells_QC_UMI.svg", width=5, height=6, units="in")
ggsave(file="All_cells_QC_UMI.png", width=5, height=6, units="in")

ggplot(subset_meta,aes_string(x="type",y="nFeature_RNA",fill="type"))+
  geom_violin(scale = "width", color="#8a8a8a")+scale_fill_manual(values=colors)+
  stat_summary(fun=median, geom="point", size=1, color="white")+ylab("Genes")+ 
  theme_classic()+
  theme(axis.title.x=element_blank(),
        ,axis.text.y=element_text(size=15, angle=90),axis.title.y=element_text(size=20,angle=90,face="bold",margin = margin(t = 10, r = 10, b = 0, l = 10),vjust=0.5),axis.text.x=element_text(size=20, face="bold", angle=90, vjust=0.5)
        ,axis.title=element_text(size=20,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 1,0, 0), "cm"))
#ggsave(file="DATE_PROJECTNAME_integrated_QI_gene_cell_types.png", height=4, width=4)
ggsave(file="All_cells_QC_gene.svg", width=5, height=6, units="in")
ggsave(file="All_cells_QC_gene.png", width=5, height=6, units="in")

Cell_type<-c("2","3","6","10","11","12","13","14","15","16","17","18","19","21","23", "5","20","1","4","8","24","7","22","9")
gene_list<-c("Stmn2" ,"Thy1", "Aldoc","Aqp4", "Mog", "Opalin", "Pdgfra","Gpr17","Tmem119","C1qc", "Cldn5","Flt1", "Foxj1","Fam216b")
colors<-c("#8a8a8a","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#fc9403","#8a8a8a","#fc9403","#8a8a8a","#fc9403","#8a8a8a","#fc9403","#8a8a8a","#8a8a8a","#fc9403","#8a8a8a","#8a8a8a","#8a8a8a")

#initialize empty data frame
Cell_number<-data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)
for (i in 1:length(Cell_type)){
  L<-length(eval(parse(text = paste("barcode_", Cell_type[i],sep=""))))
  #print(L)
  Cell_number_t<-data.frame("type"=c(rep(Cell_type[i],L)))
  #used normalized uncorrected data
  for (p in 1:length(gene_list)){
    Cell_number_t[gene_list[p]]<-as.vector(PROJECTNAME.integrated@assays$RNA@data[gene_list[p],eval(parse(text = paste("barcode_", Cell_type[i],sep="")))])
  } #end for loop
  
  #print(Cell_number_t)
  #    Cell_number["color"]<-as.character(colors[i])    
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for loop

Cell_types<-factor(Cell_number$type, levels=Cell_type)

#In dictated order
for (k in 1:length(gene_list)){
  if (k==length(gene_list)){
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x=Cell_types,y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      ,axis.text.y=element_text(size=10),axis.title.y=element_text(size=13,angle=0,face="bold",margin = margin(t = 10, r = 14, b = 0, l = 10),vjust=0.5),axis.text.x=element_text(size=13,face="bold",angle=50,vjust=0.5)
                                                                                                      ,axis.title=element_text(size=10,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")))}
  else{
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x=Cell_types,y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      axis.text.x=element_blank(),axis.text.y=element_text(size=10),axis.title.y=element_text(size=13,angle=0,face="bold",margin = margin(t = 10, r = 14, b = 0, l = 10),vjust=0.5),axis.ticks.x=element_blank(),axis.text=element_text(size=13)
                                                                                                      ,axis.title=element_text(size=10,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")))}
}#end for

#In natural order
for (k in 1:length(gene_list)){
  if (k==length(gene_list)){
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x="type",y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      ,axis.text.y=element_text(size=8),axis.title.y=element_text(size=10,angle=0,face="bold",margin = margin(t = 10, r = 14, b = 0, l = 10),vjust=0.5),axis.text.x=element_text(size=10,face="bold",angle=50,vjust=0.5)
                                                                                                      ,axis.title=element_text(size=8,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")))}
  else{
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x="type",y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      axis.text.x=element_blank(),axis.text.y=element_text(size=8),axis.title.y=element_text(size=10,angle=0,face="bold",margin = margin(t = 10, r = 14, b = 0, l = 10),vjust=0.5),axis.ticks.x=element_blank(),axis.text=element_text(size=10)
                                                                                                      ,axis.title=element_text(size=8,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")))}
}#end for

library(grid)
merge<-list()
for (i in length(gene_list):1){   
  if (length(merge)==0){
    merge<-ggplotGrob(eval(parse(text=paste("P",i,sep = ""))))
  }
  else{
    merge<-rbind(ggplotGrob(eval(parse(text=paste("P",i,sep = "")))),merge,size="last")
  }
}#end for

#all genes
#pdf(paste("Markers_stacked_violinplots_all_clusters",".pdf",sep=""),height=6, width=10)
grid.newpage()
grid.draw(merge)
ggsave(file="All_cells_violin_all_clusters.svg", plot=merge, width=8, height=8, units="in")
ggsave(file="All_cells_violin_all_clusters.png", plot=merge, width=8, height=8, units="in")
#dev.off()

Cell_type<-c("13","23","3","15","16","17","21","6","10","12","14","18","19","2","11")
gene_list<-c("Dbh","Th", "Chat", "Slc6a4", "Slc6a3", "Slc17a7", "Slc17a6", "Slc17a8","Slc32a1","Gad1","Gad2")
gray = "#8a8a8a"
orange = "#ff5500"
green = "#00aa33"
red =  "#aa0000"
blue = "deepskyblue4"
#blue = "#0000ee"
colors<-c(blue,gray,blue,orange,blue,red,red,red,blue,blue,gray,red,green,red,blue)

#initialize empty data frame
Cell_number<-data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)
for (i in 1:length(Cell_type)){
  L<-length(eval(parse(text = paste("barcode_", Cell_type[i],sep=""))))
  #print(L)
  Cell_number_t<-data.frame("type"=c(rep(Cell_type[i],L)))
  #used normalized uncorrected data
  for (p in 1:length(gene_list)){
    Cell_number_t[gene_list[p]]<-as.vector(PROJECTNAME.integrated@assays$RNA@data[gene_list[p],eval(parse(text = paste("barcode_", Cell_type[i],sep="")))])
  } #end for loop
  
  #print(Cell_number_t)
  #    Cell_number["color"]<-as.character(colors[i])    
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for loop

Cell_types<-factor(Cell_number$type, levels=Cell_type)
for (k in 1:length(gene_list)){
  if (k==length(gene_list)){
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x=Cell_types,y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      ,axis.text.y=element_text(size=10),axis.title.y=element_text(size=13,angle=0,face="bold",margin = margin(t = 10, r = 14, b = 0, l = 10),vjust=0.5),axis.text.x=element_text(size=13,face="bold",angle=50,vjust=0.5)
                                                                                                      ,axis.title=element_text(size=10,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")))}
  else{
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x=Cell_types,y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      axis.text.x=element_blank(),axis.text.y=element_text(size=10),axis.title.y=element_text(size=13,angle=0,face="bold",margin = margin(t = 10, r = 14, b = 0, l = 10),vjust=0.5),axis.ticks.x=element_blank(),axis.text=element_text(size=13)
                                                                                                      ,axis.title=element_text(size=10,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 0,0, 0), "cm")))}
}#end for

merge<-list()
for (i in length(gene_list):1){   
  if (length(merge)==0){
    merge<-ggplotGrob(eval(parse(text=paste("P",i,sep = ""))))
  }
  else{
    merge<-rbind(ggplotGrob(eval(parse(text=paste("P",i,sep = "")))),merge,size="last")
  }
}#end for

grid.newpage()
grid.draw(merge)
ggsave(file="Neurons_only_violin_all_clusters.svg", plot=merge, width=8, height=8, units="in")
ggsave(file="Neurons_only_violin_all_clusters.png", plot=merge, width=8, height=8, units="in")
#dev.off()

Cell_number<- data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)

for (i in 1:length(gene_list)){
  L<-length(Cell_type)
  Cell_number_t<- data.frame("cluster" =Cell_type, "gene"=(rep(gene_list[i],L))) # do not use c if the gene is factorized
  #used normalized uncorrected data
  for (p in 1:length(Cell_type)){
    Cell_number_t$pct[p]<-100*sum(PROJECTNAME.integrated@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_",Cell_type[p],sep="")))]>0)/length(eval(parse("barcode_",text=paste(Cell_type[p],sep=""))))
    Cell_number_t$avg[p]<-(mean(PROJECTNAME.integrated@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_",Cell_type[p],sep="")))])-mean(PROJECTNAME.integrated@assays$RNA@data[gene_list[i],]))/sd(PROJECTNAME.integrated@assays$RNA@data[gene_list[i],])
  }#end for p
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for i
#as.factor(Cell_number$cluster)

Cell_number$cluster<-factor(Cell_number$cluster,levels=c("13","23","3","15","16","17","21","6","10","12","14","18","19","2","11")) #neural clusters in the order that you want
Cell_number$gene<-factor(Cell_number$gene,levels=rev(c("Dbh","Th", "Chat", "Slc6a4", "Slc6a3", "Slc17a7", "Slc17a6", "Slc17a8","Slc32a1","Gad1","Gad2"))) #put whatever genes separate your neurons of interest here

#dev.off()
p<-ggplot(Cell_number, aes(cluster, gene)) + geom_point(aes(size = pct, colour=avg)) +  scale_x_discrete(limits = (levels(Cell_number$cluster)))+scale_y_discrete(limits =(levels(Cell_number$gene)))+
  scale_color_gradient(low = "white", high = "darkblue",limits = c(-1.5,1.5),oob=squish) + 
  geom_point(aes(size = pct), pch=21,, lwd=0,stroke=0)+ scale_size_continuous(range = c(0,6))+
  theme(axis.title.y=element_text(size=15),axis.text.y=element_text(size=15,colour = "black"),axis.text.x=element_text(size=15,angle = 50, hjust =1,colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("% cells")
png("Neurons_only_dotplot_clusters.png")
print(p)
ggsave(file="Neurons_only_dotplot_clusters.svg", plot=p, width=6, height=6, units="in")
#ggsave(file="Neurons_only_dotplot_clusters.png", plot=merge, width=6, height=6, units="in")

dev.off()

saveRDS(PROJECTNAME.integrated, file = "DATE_PROJECTNAME_integrated.rds")

F<-FeaturePlot(object =PROJECTNAME.integrated, features = c("Cnp","Foxo1","Foxo4"))
ggsave(file="/FILL_IN_DIR/DATE_PROJECTNAME_int_feature_OPC_diff.pdf",width=10,height=10)

sum(GetAssayData(object = PROJECTNAME.integrated, slot = "data")["Uty",]>0)

sum(GetAssayData(object = PROJECTNAME.integrated, slot = "data")["Xist",]>0)

sum(GetAssayData(object = PROJECTNAME.integrated, slot = "data")["Dbh",]>0)

#GetAssayData(object = PROJECTNAME.integrated, slot = "data")["Dbh",]>0 & ["Xist",]>0)

#help(FeaturePlot)

F<-FeaturePlot(object =PROJECTNAME.integrated, features = c("Stmn2","Thy1"), ncol=2,cols=c('gray','red'))
ggsave(file="/FILL_IN_DIR/DATE_PROJECTNAME_int_feature_neurons.png",width=6,height=3)

F<-FeaturePlot(object =PROJECTNAME.integrated, features = c("Aldoc","Aqp4","Tmem119","C1qc","Mog","Opalin","Pdgfra","Gpr17"), ncol=4)
ggsave(file="/FILL_IN_DIR/DATE_PROJECTNAME_int_feature_remaining_types.png",width=12,height=6)



# downstream analysis: you can start from here

PROJECTNAME.integrated<-readRDS(file = "DATE_PROJECTNAME_integrated.rds")
SAMPLE_A<-readRDS(file = "DATE_SAMPLE_A.rds")
SAMPLE_B<-readRDS(file = "DATE_SAMPLE_B.rds")
SAMPLE_C<-readRDS(file = "DATE_SAMPLE_C.rds")
DefaultAssay(PROJECTNAME.integrated) <- "RNA"

png("All_cells_QC_percent_mito.png")
FeatureScatter(object = PROJECTNAME.integrated, feature1 = "nCount_RNA", feature2 = "percent.mito",pt.size=0.2)
dev.off()
svg("All_cells_QC_percent_mito.svg")
FeatureScatter(object = PROJECTNAME.integrated, feature1 = "nCount_RNA", feature2 = "percent.mito",pt.size=0.2)
dev.off()
#ggsave(file="All_cells_QC_percent_mito.svg", width=5,height=5,units="in")
png("All_cells_QC_percent_mito_hist.png")
hist(PROJECTNAME.integrated@meta.data$percent.mito,breaks=seq(0,0.7,0.01))
dev.off()
svg("All_cells_QC_percent_mito_hist.svg")
hist(PROJECTNAME.integrated@meta.data$percent.mito,breaks=seq(0,0.7,0.01))
dev.off()
#ggsave(file="All_cells_QC_percent_mito_hist.svg", width=5,height=5,units="in")
png("All_cells_QC_RNA_count.png")
hist(PROJECTNAME.integrated@meta.data$nCount_RNA,breaks=seq(0,1000000,300),xlim = c(0,10000))
dev.off()
svg("All_cells_QC_RNA_count.svg")
hist(PROJECTNAME.integrated@meta.data$nCount_RNA,breaks=seq(0,1000000,300),xlim = c(0,10000))
dev.off()
#ggsave(file="All_cells_QC_RNA_count.svg", width=5,height=5,units="in")

# assgin clusters to cell types

head(PROJECTNAME.integrated@meta.data)

new.ident <- c("Oligo_1","Neuron_1","Neuron_2","Oligo_2","Astrocyte","Neuron_3","Microglia","OPC_1","Ambiguous_1","Neuron_4","Neuron_5","Neuron_6","Neuron_7_Dbh","Neuron_8","Neuron_9","Neuron_10","Neuron_11","Neuron_12","Neuron_13","Ambiguous_2","Neuron_14","Endothelium","Neuron_15_Chat","OPC_2")
names(x = new.ident) <- levels(x =PROJECTNAME.integrated)
PROJECTNAME.integrated<- RenameIdents(object =PROJECTNAME.integrated, new.ident)

color<-c("#19647e","#9B0C1E","#9B0C1E","#19647e","#ffc857","#9B0C1E","#9B0C1E","#9B0C1E","#4b3f72","#4b3f72","#ffc857","#9B0C1E","#676833","#9B0C1E","#b7b7b7","#E8C6C7","#4b3f72","#ffc857","#19647e","#B99A69", "#B99A69","#9B0C1E","#b7b7b7","#E8C6C7","#4b3f72","#ffc857","#19647e","#B99A69", "#B99A69")

#help(DimPlot)

#DimPlot(object = PROJECTNAME.integrated, reduction = 'umap')
DimPlot(PROJECTNAME.integrated, reduction="umap", label=TRUE, pt.size=0.5)+NoLegend()
image=DimPlot(PROJECTNAME.integrated, reduction="umap", label=TRUE, pt.size=0.5)+NoLegend()
ggsave(file="All_cells_umap_clusters_labeled.svg", plot=image, width=7, height=7, units="in")
ggsave(file="All_cells_umap_clusters_labeled.png", plot=image, width=7, height=7, units="in")



# Cell type analysis in each condition

#identifying cluster barcodes in each condition
# exp. Oligo_1_barcode_ctrl gives all barcodes in it
for (i in 1:length(new.ident)){
  assign(paste(new.ident[i],"_barcode",sep=""),colnames(PROJECTNAME.integrated@assays$RNA@data[,which(Idents(object=PROJECTNAME.integrated) %in% new.ident[i])]))# this gives all barcodes in cluster
  assign(paste(new.ident[i],"_barcode_SAMPLE_A",sep=""),intersect(colnames(SAMPLE_A@assays$RNA@data),eval(parse(text = paste(new.ident[i],"_barcode",sep="")))))
  assign(paste(new.ident[i],"_barcode_SAMPLE_B",sep=""),intersect(colnames(SAMPLE_B@assays$RNA@data),eval(parse(text = paste(new.ident[i],"_barcode",sep="")))))
  assign(paste(new.ident[i],"_barcode_SAMPLE_C",sep=""),intersect(colnames(SAMPLE_C@assays$RNA@data),eval(parse(text = paste(new.ident[i],"_barcode",sep="")))))}


DefaultAssay(PROJECTNAME.integrated) <- "integrated"

pdf(file="//FILL_IN_DIR/DATE_cluster_tree_int.pdf",width=20,height=8,paper='special') 
#https://www.rdocumentation.org/packages/ape/versions/5.2/topics/plot.phylo
PROJECTNAME.integrated <- BuildClusterTree(PROJECTNAME.integrated, verbose = FALSE, reorder = FALSE)
plot(PROJECTNAME.integrated@tools$BuildClusterTree, type = "phylogram", edge.color = "black", edge.width = 4,edge.lty = 1, srt = 0, label.offset = 10, direction = "downwards", tip.color = "black")
dev.off()

PROJECTNAME.integrated@tools

# barcode for simpler cell type
# re-create barcode for each cell type
Neuron_barcode<-c(Neuron_1_barcode,Neuron_2_barcode,Neuron_3_barcode, Neuron_4_barcode,Neuron_5_barcode,Neuron_6_barcode,Neuron_7_Dbh_barcode,Neuron_8_barcode,Neuron_9_barcode,Neuron_10_barcode,Neuron_11_barcode,Neuron_12_barcode,Neuron_13_barcode,Neuron_14_barcode,Neuron_15_Chat_barcode)
Neuron_barcode_SAMPLE_A<-intersect(Neuron_barcode,rownames(SAMPLE_A@meta.data))
Neuron_barcode_SAMPLE_B<-intersect(Neuron_barcode,rownames(SAMPLE_B@meta.data))
Neuron_barcode_SAMPLE_C<-intersect(Neuron_barcode,rownames(SAMPLE_C@meta.data))

# barcode for simpler cell type
# re-create barcode for each cell type
Astrocyte_barcode<-c(Astrocyte_barcode)
Astrocyte_barcode_SAMPLE_A<-intersect(Astrocyte_barcode,rownames(SAMPLE_A@meta.data))
Astrocyte_barcode_SAMPLE_B<-intersect(Astrocyte_barcode,rownames(SAMPLE_B@meta.data))
Astrocyte_barcode_SAMPLE_C<-intersect(Astrocyte_barcode,rownames(SAMPLE_C@meta.data))

# barcode for simpler cell type
# re-create barcode for each cell type
Microglia_barcode<-c(Microglia_barcode)
Microglia_barcode_SAMPLE_A<-intersect(Microglia_barcode,rownames(SAMPLE_A@meta.data))
Microglia_barcode_SAMPLE_B<-intersect(Microglia_barcode,rownames(SAMPLE_B@meta.data))
Microglia_barcode_SAMPLE_C<-intersect(Microglia_barcode,rownames(SAMPLE_C@meta.data))

length(Neuron_barcode)

saveRDS(Neuron_barcode_SAMPLE_A, file="/FILL_IN_DIR/DATE_Neuron_SAMPLE_A_id.rds")
saveRDS(Neuron_barcode_SAMPLE_B, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Neuron_SAMPLE_B_id.rds")
saveRDS(Neuron_barcode_SAMPLE_C, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Neuron_SAMPLE_C_id.rds")


saveRDS(Astrocyte_barcode_SAMPLE_A, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Astrocyte_SAMPLE_A_id.rds")
saveRDS(Astrocyte_barcode_SAMPLE_B, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Astrocyte_SAMPLE_B_id.rds")
saveRDS(Astrocyte_barcode_SAMPLE_C, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Astrocyte_SAMPLE_C_id.rds")

saveRDS(Microglia_barcode_SAMPLE_A, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Microglia_SAMPLE_A_id.rds")
saveRDS(Microglia_barcode_SAMPLE_B, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Microglia_SAMPLE_B_id.rds")
saveRDS(Microglia_barcode_SAMPLE_C, file="/FILL_IN_DIR/Lab_scRNAseq/DATE_Microglia_SAMPLE_C_id.rds")



## now you can take subsets of the data - here, we chose LC neurons based on expression of Dbh/Th and GABA neurons based on expression of GABA markers

print(PROJECTNAME.integrated)

DefaultAssay(object = PROJECTNAME.integrated) <- "RNA"

DBH_or_TH<- subset(PROJECTNAME.integrated, (Dbh>0 | Th>0))
DBH_or_TH<- NormalizeData(object=DBH_or_TH ,normalization.method="LogNormalize", scale.factor=10000, verbose = FALSE) 
DBH_or_TH<- FindVariableFeatures(object=DBH_or_TH,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(DBH_or_TH), 10)
plot1<-VariableFeaturePlot(DBH_or_TH)
plot1
plot2<-LabelPoints(plot=plot1, points=top10, repel=TRUE)
plot2

print(DBH_or_TH)

DBH_or_TH<- ScaleData(object=DBH_or_TH, features = rownames(x=DBH_or_TH), vars.to.regress = c("nCount_RNA", "percent.mito"), verbose=TRUE)
DBH_or_TH<- RunPCA(object=DBH_or_TH, features = VariableFeatures(object=DBH_or_TH), verbose = TRUE)
DBH_or_TH <- JackStraw(DBH_or_TH,num.replicate=100)
DBH_or_TH <-ScoreJackStraw(DBH_or_TH, dims=1:20)
JackStrawPlot(DBH_or_TH, dims=1:20)
image=JackStrawPlot(DBH_or_TH, dims=1:20)
ggsave(file="JackStrawPlot_All_DBH_or_TH_cells.svg", plot=image, width=6, height=6, units="in")
ElbowPlot(DBH)
image=ElbowPlot(DBH)
ggsave(file="Elbowplot_All_DBH_or_TH_cells.svg", plot=image, width=6, height=6, units="in")

DBH_or_TH<- FindNeighbors(object=DBH_or_TH, dims = 1:30)
DBH_or_TH<- FindClusters(object=DBH_or_TH, resolution = 0.8)
DBH_or_TH<- RunUMAP(DBH_or_TH, reduction='pca', dims=1:30, verbose=TRUE)

saveRDS(DBH_or_TH,file = "DATE_All_DBH_or_TH_dataset.rds")





DBH_or_TH<-readRDS(file = "DATE_All_DBH_or_TH_dataset.rds")

print(DBH_or_TH)

DefaultAssay(object = DBH_or_TH) <- "integrated" 

#green = "#227700"
#red =  "#aa0000"
#blue = "#0011ff"
#orange = "#ff5500"
new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11")
names(new.cluster.ids) <- levels(DBH_or_TH)
DBH_or_TH <- RenameIdents(DBH_or_TH, new.cluster.ids)
p<-DimPlot(DBH_or_TH, reduction="umap", label=TRUE, pt.size=2)+NoLegend()
ggsave(file="All_DBH_or_TH_umap_labeled.svg", plot=p, width=5, height=5, units="in")
ggsave(file="All_DBH_or_TH_umap_labeled.png", plot=p, width=5, height=5, units="in")
p<-DimPlot(DBH_or_TH, reduction="umap", label=FALSE, pt.size=2)+NoLegend()
ggsave(file="All_DBH_or_TH_umap_unlabeled.svg", plot=p, width=4, height=4, units="in")
ggsave(file="All_DBH_or_TH_umap_unlabeled.png", plot=p, width=4, height=4, units="in")
DimPlot(DBH_or_TH, reduction="umap",pt.size=2)

DefaultAssay(object = DBH_or_TH) <- "RNA" 

table(Idents(DBH_or_TH))

counts<-as.matrix(DBH_or_TH@assays$RNA@data)
write.table(data.frame("GENE"=rownames(counts),counts),file="counts_nondoublets_DBH_or_TH.txt",row.names=FALSE,sep="\t")
markers <- FindAllMarkers(object = DBH_or_TH, only.pos = TRUE, min.pct = 0.25)
top_50<-markers %>% group_by(cluster) %>% top_n(50)
write.table(data.frame("test"=as.character(rownames(top_50)),top_50),file="Top50Genes_nondoublets_DBH_or_TH.txt",row.names=FALSE,col.names=c("",colnames(top_50)),sep="\t",eol="\n")
cluster<-Idents(object=DBH_or_TH)
cluster<-as.matrix(cluster)
cluster[,1]<-as.character(cluster[,1])
cluster[,0]<-as.character(cluster[,0])
cluster<-data.frame("x"=rownames(cluster),cluster)
write.table(cluster,file="Cluster_nodoublets_DBH_or_TH.txt",row.names=FALSE,col.names=c("","x"),sep="\t",eol="\n")

#HVFInfo(DBH_or_TH)[VariableFeatures(DBH_or_TH),]

## in these next few plots, the genes were chosen based on differentially expressed genes between clusters

image=VlnPlot(DBH_or_TH, features = c("Dbh","Th","Zmat4","Ttc3","Trpc3","Trf","Ubxn4","Tshz2","Tril","Zic1","Vcan","Vmn1r209","Trem2"))
image
ggsave(file="All_int_DBH_or_TH_Markers_neurons.svg", plot=image, width=14, height=6, units="in")

new.ident <- c("1","2","3","4","5","6","7","8","9","10","11")
#new.ident <- c("DBH_TH_1","DBH_TH_2","DBH_TH_3","DBH_TH_4","DBH_TH_5","DBH_TH_6","DBH_TH_7","DBH_TH_8","DBH_TH_9","DBH_TH_10","DBH_TH_11")
names(new.cluster.ids) <- levels(DBH_or_TH)
DBH_or_TH <- RenameIdents(DBH_or_TH, new.cluster.ids)
new.ident=new.cluster.ids
for (i in 1:length(new.ident)){
  assign(paste("barcode_DBH_TH_",new.ident[i],sep=""),colnames(DBH_or_TH@assays$RNA@data[,which(Idents(object=DBH_or_TH) %in% new.ident[i])]))# this gives all barcodes in cluster
}#endfor

#Cell_type <- c("1","2","3","4","5","6","7","8","9","10","11")
Cell_type<-c("DBH_TH_1","DBH_TH_2","DBH_TH_3","DBH_TH_4","DBH_TH_5","DBH_TH_6","DBH_TH_7","DBH_TH_8","DBH_TH_9","DBH_TH_10","DBH_TH_11")
gene_list<-c("Dbh","Th","Zmat4","Ttc3","Trpc3","Trf","Ubxn4","Tshz2","Tril","Zic1","Vcan","Vmn1r209","Trem2")
colors<-c(blue,red,blue,red,blue,red,blue,red,blue,red,blue,red,blue,red,blue,red)

#initialize empty data frame
Cell_number<-data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)
for (i in 1:length(Cell_type)){
  L<-length(eval(parse(text = paste("barcode_", Cell_type[i],sep=""))))
  #print(L)
  Cell_number_t<-data.frame("type"=c(rep(Cell_type[i],L)))
  #used normalized uncorrected data
  for (p in 1:length(gene_list)){
    #print(p)
    Cell_number_t[gene_list[p]]<-as.vector(DBH_or_TH@assays$RNA@data[gene_list[p],eval(parse(text = paste("barcode_", Cell_type[i],sep="")))])
  } #end for loop
  
  #print(Cell_number_t)
  #    Cell_number["color"]<-as.character(colors[i])    
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for loop

Cell_types<-factor(Cell_number$type, levels=Cell_type)
for (k in 1:length(gene_list)){
  if (k==length(gene_list)){
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x=Cell_types,y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      ,axis.text.y=element_text(size=10),axis.title.y=element_text(size=15,angle=0,face="bold",margin = margin(t = 10, r = 20, b = 0, l = 2),vjust=0.5),axis.text.x=element_text(size=13,face="bold",angle=50,vjust=0.5)
                                                                                                      ,axis.title=element_text(size=10,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 2,0, 1), "cm")))}
  else{
    assign(paste("P",k,sep=""),ggplot(Cell_number,aes_string(x=Cell_types,y=gene_list[k],fill="type"))+geom_violin(scale = "width")+scale_fill_manual(values=colors)
           +stat_summary(fun=median, geom="point", size=0.6, color="white")+ylab(gene_list[k])+ theme(axis.title.x=element_blank(),
                                                                                                      axis.text.x=element_blank(),axis.text.y=element_text(size=10),axis.title.y=element_text(size=15,angle=0,face="bold",margin = margin(t = 10, r = 20, b = 0, l = 2),vjust=0.5),axis.ticks.x=element_blank(),axis.text=element_text(size=10)
                                                                                                      ,axis.title=element_text(size=10,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="none",plot.margin = unit(c(0, 2, 0, 1), "cm")))}
}#end for

merge<-list()
for (i in length(gene_list):1){   
  if (length(merge)==0){
    merge<-ggplotGrob(eval(parse(text=paste("P",i,sep = ""))))
  }
  else{
    merge<-rbind(ggplotGrob(eval(parse(text=paste("P",i,sep = "")))),merge,size="last")
  }
}#end for

#all genes
#pdf(paste("Markers_stacked_violinplots_all_clusters",".pdf",sep=""),height=6, width=10)
grid.newpage()
grid.draw(merge)
ggsave(file="All_DBH_or_TH_subtype_clusters.svg", plot=merge, width=5, height=6, units="in")
#dev.off()

Cell_number<- data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)

for (i in 1:length(gene_list)){
  L<-length(Cell_type)
  Cell_number_t<- data.frame("cluster" =Cell_type, "gene"=(rep(gene_list[i],L))) # do not use c if the gene is factorizsed
  #used normalized uncorrected data
  for (p in 1:length(Cell_type)){
    Cell_number_t$pct[p]<-100*sum(DBH_or_TH@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_",Cell_type[p],sep="")))]>0)/length(eval(parse(text=paste("barcode_",Cell_type[p],sep=""))))
    Cell_number_t$avg[p]<-(mean(DBH_or_TH@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_", Cell_type[p],sep="")))])-mean(DBH_or_TH@assays$RNA@data[gene_list[i],]))/sd(DBH_or_TH@assays$RNA@data[gene_list[i],])
  }#end for p
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for i
#as.factor(Cell_number$cluster)

Cell_number$cluster<-factor(Cell_number$cluster,levels=c("DBH_TH_1","DBH_TH_2","DBH_TH_3","DBH_TH_4","DBH_TH_5","DBH_TH_6","DBH_TH_7","DBH_TH_8","DBH_TH_9","DBH_TH_10","DBH_TH_11"))
Cell_number$gene<-factor(Cell_number$gene, levels=rev(c("Dbh","Th","Zmat4","Ttc3","Trpc3","Trf","Ubxn4","Tshz2","Tril","Zic1","Vcan","Vmn1r209","Trem2")))

ggplot(Cell_number, aes(cluster, gene)) + geom_point(aes(size = pct, colour=avg)) +  scale_x_discrete(limits = (levels(Cell_number$cluster)))+scale_y_discrete(limits =(levels(Cell_number$gene)))+
  scale_color_gradient(low = "white", high = "darkorange1",limits = c(-1.5,1.5),oob=squish) + 
  geom_point(aes(size = pct), pch=21,, lwd=0,stroke=0)+ scale_size_continuous(range = c(0,6))+
  theme(axis.title.y=element_text(size=20, face='bold'),axis.text.y=element_text(size=15,colour = "black"),axis.title.x=element_text(size=20, face='bold'),axis.text.x=element_text(size=15,angle = 50, hjust =1,colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("% cells")
ggsave(file="All_DBH_or_TH_top_genes_dotplot.svg",height=7, width=7)
ggsave(file="All_DBH_or_TH_top_genes_dotplot.png",height=7, width=7)



#example dotplot - neuropeptides
Cell_type<-c("DBH_TH_1","DBH_TH_2","DBH_TH_3","DBH_TH_4","DBH_TH_5","DBH_TH_6","DBH_TH_7","DBH_TH_8","DBH_TH_9","DBH_TH_10","DBH_TH_11")
gene_list<-c("Calca", "Cck", "Gal", "Hcrt", "Npy", "Nts", "Pdyn", "Penk", "Pnoc", "Pomc", "Sst", "Tac1","Tac2", "Vip", "Trh")

Cell_number<- data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)

for (i in 1:length(gene_list)){
  L<-length(Cell_type)
  Cell_number_t<- data.frame("cluster" =Cell_type, "gene"=(rep(gene_list[i],L))) # do not use c if the gene is factorizsed
  #used normalized uncorrected data
  for (p in 1:length(Cell_type)){
    Cell_number_t$pct[p]<-100*sum(DBH_or_TH@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_",Cell_type[p],sep="")))]>0)/length(eval(parse(text=paste("barcode_",Cell_type[p],sep=""))))
    Cell_number_t$avg[p]<-(mean(DBH_or_TH@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_", Cell_type[p],sep="")))])-mean(DBH_or_TH@assays$RNA@data[gene_list[i],]))/sd(DBH_or_TH@assays$RNA@data[gene_list[i],])
  }#end for p
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for i
#as.factor(Cell_number$cluster)
Cell_number$cluster<-factor(Cell_number$cluster,levels=c("DBH_TH_1","DBH_TH_2","DBH_TH_3","DBH_TH_4","DBH_TH_5","DBH_TH_6","DBH_TH_7","DBH_TH_8","DBH_TH_9","DBH_TH_10","DBH_TH_11"))
Cell_number$gene<-factor(Cell_number$gene, levels=rev(gene_list))
ggplot(Cell_number, aes(cluster, gene)) + geom_point(aes(size = pct, colour=avg)) +  scale_x_discrete(limits = (levels(Cell_number$cluster)))+scale_y_discrete(limits =(levels(Cell_number$gene)))+
  scale_color_gradient(low = "white", high = "darkorange1",limits = c(-1,1),oob=squish) + 
  geom_point(aes(size = pct), pch=21,, lwd=0,stroke=0)+ scale_size_continuous(range = c(0,6))+
  theme(axis.title.y=element_text(size=20, face='bold'),axis.text.y=element_text(size=15,colour = "black"),axis.title.x=element_text(size=20, face='bold'),axis.text.x=element_text(size=15,angle = 50, hjust =1,colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle("% cells")
#ggsave(file="All_int_DBH_or_TH_top_genes_dotplot.pdf",height=6, width=16 , paper = "letter")
ggsave(file="All_DBH_or_TH_neuropeptides_dotplot.svg",height=7, width=7)
ggsave(file="All_DBH_or_TH_neuropeptides_dotplot.png",height=7, width=7)



#example dotplot - load from .xlsx file...these are convenient to download from gene ontology database and edit, but long lists can take a long time to run
#load gene names
all_genes_of_interest<-read.xlsx2(file="/FILL_IN_DIR/genes_of_interest.xlsx", 1,stringsAsFactors = FALSE)
gene_list<-intersect(unique(all_genes_of_interest$Symbol),rownames(DBH_or_TH@assays$RNA@data)) #takes out duplicates and genes that aren't in dataset


fig_header<-"genes_of_interest"

#confirm gene list
fig_header
gene_list
length(gene_list)
sizefactor<-(length(gene_list)*0.25)
sizefactor

#
Cell_number<- data.frame(Date=as.Date(character()),File=character(),User=character(),stringsAsFactors=FALSE)
list_size<-length(gene_list)

L<-length(Cell_type)

for (i in 1:length(gene_list)){
  Cell_number_t<- data.frame("cluster"=Cell_type, "gene"=(rep(gene_list[i],L))) # do not use c if the gene is factorized
  #used normalized uncorrected data
  for (p in 1:length(Cell_type)){
    Cell_number_t$pct[p]<-100*sum(DBH_or_TH@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_",Cell_type[p],sep="")))]>0)/length(eval(parse(text=paste("barcode_",Cell_type[p],sep=""))))
    Cell_number_t$avg[p]<-(mean(DBH_or_TH@assays$RNA@data[gene_list[i],eval(parse(text=paste("barcode_", Cell_type[p],sep="")))])-mean(DBH_or_TH@assays$RNA@data[gene_list[i],]))/sd(DBH_or_TH@assays$RNA@data[gene_list[i],])
    #print(Cell_number_t$pct[p])
    #print(Cell_number_t$avg[p])
  }#end for p
  Cell_number<-rbind(Cell_number_t,Cell_number)
}#end for i
Cell_number$cluster<-factor(Cell_number$cluster,levels=rev(Cell_type))
Cell_number$gene<-factor(Cell_number$gene, levels=rev(gene_list))



ggplot(Cell_number, aes(gene, cluster)) + geom_point(aes(size = pct, colour=avg)) +  scale_x_discrete(limits = rev(levels(Cell_number$gene)))+scale_y_discrete(limits =rev(levels(Cell_number$cluster)))+
  scale_color_gradient(low = "white", high = "springgreen4",limits = c(-1,1),oob=squish) + 
  geom_point(aes(size = pct), pch=21,, lwd=0,stroke=0)+ scale_size_continuous(range = c(0,7),limits=c(0,100),breaks=seq(0,100,25))+
  theme(axis.title.y=element_text(size=15),axis.text.y=element_text(size=10,colour = "black"),axis.text.x=element_text(size=8,angle = 90, hjust =1,colour = "black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ggtitle(fig_header)
ggsave((filename = paste0(fig_header,".pdf")),height=7, width=sizefactor+2, limitsize = FALSE)