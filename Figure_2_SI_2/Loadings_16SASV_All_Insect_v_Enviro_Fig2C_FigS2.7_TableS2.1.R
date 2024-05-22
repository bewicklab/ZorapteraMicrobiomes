library('vegan')
library('phyloseq')
library('stringr')
library('microshades')
library('speedyseq')
library('forcats')
library('cowplot')
library('knitr')
library('ggplot2')
library('PERMANOVA')
library('pairwiseAdonis')
library('ggsignif')
library('tidyverse')
library('ggpubr')
library('betapart')
library('stats')
library('ape')
library('HybridMicrobiomes')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make OTU Table

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Read in CSV file for OTU table (table.csv for ASVs, tablelevel6.csv for genera)
data <- import_biom('16S/ASV_Table.biom')

#Find ASV table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('16S/Zoraptera_metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('16S/rooted_tree_out/tree.nwk'))

#Rename the tree tips to match the ASV table names
cut_names<-c()
for (k in 1:length(taxa_names(tree_file))){
  #Depending on the Qiime2 output, you may need to split the names with a space or with an underscore
  #cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],' ')[[1]][1])
  cut_names<-c(cut_names,strsplit(taxa_names(tree_file)[k],'_')[[1]][1])
}
taxa_names(tree_file)<-cut_names

#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP,tree_file)


#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('16S/ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)
chloroplasts<-rownames(tax_table(ZNC_physeq)[,1][which(tax_table(ZNC_physeq)[,3]!='c__Chloroplast')])
ZNC_physeq = prune_taxa(chloroplasts, ZNC_physeq)

#Just look at samples with more than 10000 reads
ZNC_physeq<-prune_samples(colSums(otu_table(ZNC_physeq))>7000,ZNC_physeq)

#Rarefy the OTU table specifically for this bootstrap sample
rtempZOTU<-rarefy_even_depth(ZNC_physeq,verbose=FALSE,rngseed = 1)

#Remove taxa that drop out of the rarefied OTU table (this makes UniFrac run more smoothly... and you don't need them anyhow)
rZOTU<-prune_taxa(row.names(tax_table(rtempZOTU)),rtempZOTU)

#Change the sample names from the Plate IDs to the Colony IDs (unique to Zoraptera... because of our naming system... you may need to do something else)
sample_names(rZOTU) <- sample_data(rZOTU)$Colony

#Find the treatment group for each of your samples
type<-sample_data(rZOTU)$Type

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('Adult','Environment')

#Convert names of groups to numbers
no_type<-rep(1,length(type))
for (k in 1:length(types)){
  no_type[which(type==types[k])]<-k
}

colvec<-c("wheat2", "brown")

set.seed(1)

#Find the euclidean, jaccard, bray, UniFrac and weighted-UniFrac indices
euclidean<-vegdist(t(otu_table(rZOTU)),method='euclidean',upper=TRUE,diag=TRUE,binary=FALSE)
jaccard<-vegdist(t(otu_table(rZOTU)),method='jaccard',upper=TRUE,diag=TRUE,binary=TRUE)
bray<-vegdist(t(otu_table(rZOTU)),method='bray',upper=TRUE,diag=TRUE, binary=FALSE)
unifrac<-UniFrac(rZOTU,weighted=FALSE)
wunifrac<-UniFrac(rZOTU,weighted=TRUE)

#Convert the 'dist' type of object that you get from vegdist to a matrix with labelled rows and columns
euclidean_matrix<-as.matrix(euclidean,labels=TRUE)
bray_matrix<-as.matrix(bray,labels=TRUE)
jaccard_matrix<-as.matrix(jaccard,labels=TRUE)
unifrac_matrix<-as.matrix(unifrac,labels=TRUE)
wunifrac_matrix<-as.matrix(wunifrac,labels=TRUE)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform 2D PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Run PCA analysis
pca <- prcomp(t(otu_table(rZOTU)))

#Find variation explained by axes
variance <- (pca$sdev)^2
varPercent <- variance/sum(variance) * 100

#Find loadings on axes
loadings <- pca$rotation

#Find five largest loadings for PCA1
loads1<-loadings[,1]*loadings[,1]*100
toploadings1<-sort(loadings[,1]*loadings[,1],decreasing=TRUE)[1:5]*100
topnames1<-names(toploadings1)
toptaxa1<-c()
topdirection1<-c()
for (k in 1:5){
  temp<-which(assigned_taxa[,1]==topnames1[k])
  genus<-strsplit(strsplit(assigned_taxa[temp,2],';')[[1]][6],'g__')[[1]][2]
  species<-strsplit(strsplit(assigned_taxa[temp,2],';')[[1]][7],'s__')[[1]][2]
  if (genus=="NA" || genus == 'unidentified'){
    genus<-assigned_taxa[temp,1]
    species<-assigned_taxa[temp,2]
  }
  else{
    print(genus)
  }
  toptaxa1<-c(toptaxa1,paste(genus,species))
  temp<-which(names(loads1)==topnames1[k])
  topdirection1<-c(topdirection1,(loadings[temp,1]))
}

#Find five largest loadings for PCA2
loads2<-loadings[,2]*loadings[,2]*100
toploadings2<-sort(loadings[,2]*loadings[,2],decreasing=TRUE)[1:5]*100
topnames2<-names(toploadings2)
toptaxa2<-c()
topdirection2<-c()
for (k in 1:5){
  temp<-which(assigned_taxa[,1]==topnames2[k])
  genus<-strsplit(strsplit(assigned_taxa[temp,2],';')[[1]][6],'g__')[[1]][2]
  species<-strsplit(strsplit(assigned_taxa[temp,2],';')[[1]][7],'s__')[[1]][2]
  if (genus=="NA" || genus == 'unidentified'){
    genus<-assigned_taxa[temp,1]
    species<-assigned_taxa[temp,2]
  }
  toptaxa2<-c(toptaxa2,paste(genus,species))
  temp<-which(names(loads2)==topnames2[k])
  topdirection2<-c(topdirection2,(loadings[temp,2]))
  
}

tops<-data.frame(toptaxa1,topdirection1,toploadings1,toptaxa2,topdirection2,toploadings2)
rownames(tops)<-c('1','2','3','4','5')

write.csv(x=tops,file='DominantLoadingsAll16S.csv')
