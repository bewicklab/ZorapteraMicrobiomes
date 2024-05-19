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
data <- import_biom('ITS/ASV_Table.biom')

#Find ASV table in biom file
OTU_biom<-otu_table(data)

#Find the taxonomy table in the biom file (this exists for the ASV biom file!)
TAX<-tax_table(data)

#Read in metadata
meta_csv<-read.csv('ITS/Zoraptera_metadata.csv', row.names=1,header= TRUE)

#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Read in the phylogenetic tree that you made using Qiime2 (from the sv.seqs.fna file given to you by Zymo)
tree_file<-multi2di(read.tree('ITS/rooted_tree_out/tree.nwk'))

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
assigned_taxa<-read.csv('ITS/ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
ZNC_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)

#Just look at samples with more than 10000 reads
ZNC_physeq<-prune_samples(colSums(otu_table(ZNC_physeq))>7000,ZNC_physeq)

#Rarefy the OTU table specifically for this bootstrap sample
rtempZOTU<-rarefy_even_depth(ZNC_physeq,verbose=FALSE,rngseed = 1)

#Remove taxa that drop out of the rarefied OTU table (this makes UniFrac run more smoothly... and you don't need them anyhow)
rZOTU<-prune_taxa(row.names(tax_table(rtempZOTU)),rtempZOTU)

#Change the sample names from the Plate IDs to the Colony IDs (unique to Zoraptera... because of our naming system... you may need to do something else)
sample_names(rZOTU) <- sample_data(rZOTU)$Colony

#Find the treatment group for each of your samples
typeA<-sample_data(rZOTU)$Type
typeA[typeA=="Adult"]<-"Insect"
typeA[typeA=="Juvenile"]<-"Insect"
sample_data(rZOTU)$TypeA<-typeA
type<-sample_data(rZOTU)$TypeA

#A list of your treatment groups (this should be the order you intend to use for plotting and the color vector as well)
types<-c('Insect','Environment')

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

#Find coordinates for plotting (identical to PCoA on Euclidean distances)
scores <- pca$x
scores_2<-scores[,1:2]

pcoa_euclidean<-pcoa(euclidean)$vectors[,1:2]
pcoa_jaccard<-pcoa(jaccard)$vectors[,1:2]
pcoa_bray<-pcoa(bray)$vectors[,1:2]
pcoa_unifrac<-pcoa(unifrac)$vectors[,1:2]
pcoa_wunifrac<-pcoa(wunifrac)$vectors[,1:2]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

par(mar = c(5, 5,5, 5))
plot(pcoa_euclidean,type='n',cex.lab=1.75,xlab='PCA1',ylab='PCA2')
ordihull(pcoa_euclidean,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_euclidean, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_euclidean[which(!(rownames(pcoa_euclidean) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_euclidean<-c()    #Value of the centroid along the first NMDS axis
ymeans_euclidean<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_euclidean<-c(xmeans_euclidean,mean(pcoa_euclidean[pts,1]))
  ymeans_euclidean<-c(ymeans_euclidean,mean(pcoa_euclidean[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_euclidean[k],ymeans_euclidean[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_euclidean[k],ymeans_euclidean[k],col='black',pch=1,cex=2)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA Jaccard

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Jaccard PCoA (equivalent to PCA)
plot(pcoa_jaccard,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_jaccard,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_jaccard, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_jaccard[which(!(rownames(pcoa_jaccard) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_jaccard<-c()    #Value of the centroid along the first NMDS axis
ymeans_jaccard<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_jaccard<-c(xmeans_jaccard,mean(pcoa_jaccard[pts,1]))
  ymeans_jaccard<-c(ymeans_jaccard,mean(pcoa_jaccard[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_jaccard[k],ymeans_jaccard[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_jaccard[k],ymeans_jaccard[k],col='black',pch=1,cex=2)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA Bray-Curtis

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Bray PCoA (equivalent to PCA)
plot(pcoa_bray,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_bray,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_bray, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_bray[which(!(rownames(pcoa_bray) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")


#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_bray<-c()    #Value of the centroid along the first NMDS axis
ymeans_bray<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_bray<-c(xmeans_bray,mean(pcoa_bray[pts,1]))
  ymeans_bray<-c(ymeans_bray,mean(pcoa_bray[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_bray[k],ymeans_bray[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_bray[k],ymeans_bray[k],col='black',pch=1,cex=2)
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the Unifrac PCoA (equivalent to PCA)
plot(pcoa_unifrac,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_unifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_unifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_unifrac[which(!(rownames(pcoa_unifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")

#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_unifrac<-c()    #Value of the centroid along the first NMDS axis
ymeans_unifrac<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_unifrac<-c(xmeans_unifrac,mean(pcoa_unifrac[pts,1]))
  ymeans_unifrac<-c(ymeans_unifrac,mean(pcoa_unifrac[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_unifrac[k],ymeans_unifrac[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_unifrac[k],ymeans_unifrac[k],col='black',pch=1,cex=2)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Plot the PCoA weighted-UniFrac

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Plot the weighted Unifrac PCoA (equivalent to PCA)
plot(pcoa_wunifrac,type='n',cex.lab=1.75,xlab='PCoA1',ylab='PCoA2')
ordihull(pcoa_wunifrac,groups=no_type,draw="polygon",col= colvec[1:length(types)], label=F,alpha=0.25)
points(pcoa_wunifrac, display = "sites", pch=16, col = colvec[no_type])
points(pcoa_wunifrac[which(!(rownames(pcoa_wunifrac) %in% c("G.Env","P.Env","L.Env","Q.Env","I.Env","J.Env"))),],  pch=1, col = "black")


#Find the centroid (based on mean) of each treatment group in the NMDS plot
xmeans_wunifrac<-c()    #Value of the centroid along the first NMDS axis
ymeans_wunifrac<-c()    #Value of the centroid along the second NMDS axis

#For each treatment group...
for (k in 1:length(types)){
  #Find the points corresponding to that group
  pts<-which(type==types[k])
  #Find the averages values of the points in that treatment group along each axis
  xmeans_wunifrac<-c(xmeans_wunifrac,mean(pcoa_wunifrac[pts,1]))
  ymeans_wunifrac<-c(ymeans_wunifrac,mean(pcoa_wunifrac[pts,2]))
}

#Plot the centroid for each group overtop of your NMDS scatterplot
for (k in 1:length(types)){
  points(xmeans_wunifrac[k],ymeans_wunifrac[k],col=colvec[k],pch=16,cex=2)
  points(xmeans_wunifrac[k],ymeans_wunifrac[k],col='black',pch=1,cex=2)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Euclidean PERMANOVA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_euclidean<-list(Data=t(otu_table(rZOTU)),D=euclidean_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_euclidean=PERMANOVA(inputpermanova_euclidean, factor(type),nperm=1000)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#           Perform Jaccard PERMANOVA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_jaccard<-list(Data=t(otu_table(rZOTU)),D=jaccard_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_jaccard=PERMANOVA(inputpermanova_jaccard, factor(type),nperm=1000)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Bray Curtis Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_bray<-list(Data=t(otu_table(rZOTU)),D=bray_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_bray=PERMANOVA(inputpermanova_bray, factor(type),nperm=1000)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on UniFrac Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_unifrac<-list(Data=t(otu_table(rZOTU)),D=unifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_unifrac=PERMANOVA(inputpermanova_unifrac, factor(type),nperm=1000)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Perform PERMANOVA Tests on Bray Curtis Distances

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Make a list of the OTU table and the distance matrix (input into PERMANOVA package)
inputpermanova_wunifrac<-list(Data=t(otu_table(rZOTU)),D=wunifrac_matrix,Coefficient="Other")

#Perform PERMANOVA comparing treatment groups
permanova_by_group_wunifrac=PERMANOVA(inputpermanova_wunifrac, factor(type),nperm=1000)

Fstats<-c(permanova_by_group_euclidean$Initial$Global[5],permanova_by_group_jaccard$Initial$Global[5],permanova_by_group_bray$Initial$Global[5],permanova_by_group_unifrac$Initial$Global[5],permanova_by_group_wunifrac$Initial$Global[5])
pvalues<-c(permanova_by_group_euclidean$pvalue,permanova_by_group_jaccard$pvalue,permanova_by_group_bray$pvalue,permanova_by_group_unifrac$pvalue,permanova_by_group_wunifrac$pvalue)

permanova_stats<-data.frame(Fstats,pvalues)
rownames(permanova_stats)<-c('Euclidean','Jaccard','Bray-Curtis','Unifrac','weighted Unifrac')

write.csv(x=permanova_stats,file='AllITSPermanova.csv')

#Test for homogeneity of variances 
bd<-betadisper(euclidean,as.factor(type))
dispersion_anova_euclidean<-anova(bd)

bd<-betadisper(jaccard,as.factor(type))
dispersion_anova_jaccard<-anova(bd)

bd<-betadisper(bray,as.factor(type))
dispersion_anova_bray<-anova(bd)

bd<-betadisper(unifrac,as.factor(type))
dispersion_anova_unifrac<-anova(bd)

bd<-betadisper(wunifrac,as.factor(type))
dispersion_anova_wunifrac<-anova(bd)

Fs<-c(dispersion_anova_euclidean$`F value`[1],dispersion_anova_jaccard$`F value`[1],dispersion_anova_bray$`F value`[1],dispersion_anova_unifrac$`F value`[1],dispersion_anova_wunifrac$`F value`[1])
ps<-c(dispersion_anova_euclidean$`Pr(>F)`[1],dispersion_anova_jaccard$`Pr(>F)`[1],dispersion_anova_bray$`Pr(>F)`[1],dispersion_anova_unifrac$`Pr(>F)`[1],dispersion_anova_wunifrac$`Pr(>F)`[1])

dispersion_stats<-data.frame(Fs,ps)
rownames(dispersion_stats)<-c('Euclidean','Jaccard','Bray-Curtis','Unifrac','weighted Unifrac')

write.csv(x=dispersion_stats,file='AllITSDispersion.csv')





