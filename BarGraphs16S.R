library(vegan)
library(phyloseq)
library(ggplot2)
library(forcats)
library(RColorBrewer)
library(microshades)
library(speedyseq)
library(cowplot)
library(stringr)
library(knitr)

#Installation of microshades and speedyseq is a bit different; all the other packages can be installed as usual
#remotes::install_github("KarstensLab/microshades")
#remotes::install_github("mikemc/speedyseq")


#%%%%%%%%%%%%%%%%%%%%%%   GET THE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Read in biom file
zoraptera_biom<-import_biom('16S/Composition_Summary_L6.biom')
#Find OTU table in biom file
OTU_biom<-otu_table(zoraptera_biom)

#Read in the readcounts in each sample (To get this, open the txt file and then save the file as a csv)
zoraptera_read_processing<-read.csv('16S/Read_Processing_Summary.csv', row.names=1,header= TRUE)

#Find the sample names in the biom file
biom_label<-colnames(OTU_biom)

#Find the sample names in the readcount file (in this case they were numbers so I converted them to strings to match the biom and csv sample labels)
read_summary_label<-as.character(zoraptera_read_processing$customer_label)

#Find the reads per sample in the readcount file
read_summary_readcount<-zoraptera_read_processing$seqs.after_size_filtration

#For every column in the biom OTU table...
for (k in 1:length(biom_label)){
  #...find the corresponding label in the readcount file...
  find_read_summary_row<-which(read_summary_label==biom_label[k])
  #...find the corresponding number of reads from that sample, multiply the column by the total readcount and round to the nearest integer
  OTU_biom[,k]<-round(OTU_biom[,k]*read_summary_readcount[find_read_summary_row])
}

#Remove the non-bacteria/archaea (make sure the None;Other;etc. row is the FIRST row in the table... it usually is)
badTaxa = "None;Other;Other;Other;Other;Other"
allTaxa = taxa_names(OTU_biom)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
OTU_biom = prune_taxa(allTaxa, OTU_biom)

#Remove samples with reads below a threshold (here it is 10000)
OTU_biom<-OTU_biom[,-which(colSums(OTU_biom)<7000)]

#Here I'm showing rarification and that it's the same... we won't make the phyloseq object from the rarified tables... but we will rarify later on...
#This is more just for illustration so you see that the biom and the csv files give the same thing
set.seed(1)


#Create a taxonomy table (the zymo biom file doesn't seem to include this...)
domain<-c()
phylum<-c()
class<-c()
order<-c()
family<-c()
genus<-c()
species<-c()
#Extract the taxonomic name at each taxonomic level
for (k in 1:length(OTU_biom[,1])){
  temp<-strsplit(rownames(OTU_biom)[k],split=';')
  domain<-c(domain,substr(temp[[1]][1],4,50))
  if (temp[[1]][2]=="p__NA"){temp[[1]][2]<-paste0('unclassified_',substr(temp[[1]][1],4,50))}
  else {temp[[1]][2]<-substr(temp[[1]][2],4,50)}
  phylum<-c(phylum,temp[[1]][2])
  if (temp[[1]][3]=="c__NA"){temp[[1]][3]<-paste0('unclassified_',str_replace(temp[[1]][2],'unclassified_',''))}
  else {temp[[1]][3]<-substr(temp[[1]][3],4,50)}
  class<-c(class,temp[[1]][3])
  if (temp[[1]][4]=="o__NA"){temp[[1]][4]<-paste0('unclassified_',str_replace(temp[[1]][3],'unclassified_',''))}
  else {temp[[1]][4]<-substr(temp[[1]][4],4,50)}
  order<-c(order,temp[[1]][4])
  if (temp[[1]][5]=="f__NA"){temp[[1]][5]<-paste0('unclassified_',str_replace(temp[[1]][4],'unclassified_',''))}
  else {temp[[1]][5]<-substr(temp[[1]][5],4,50)}
  family<-c(family,temp[[1]][5])
  if (temp[[1]][6]=="g__NA"){temp[[1]][6]<-paste0('unclassified_',str_replace(temp[[1]][5],'unclassified_',''))}
  else {temp[[1]][6]<-substr(temp[[1]][6],4,50)}
  genus<-c(genus,temp[[1]][6])
  species<-c(species,temp[[1]][7])
}
#bind the taxonomic level names together into a table/matrix
taxmat = cbind(domain,phylum,class,order,family,genus,species)
#The rownames are the rownames from the OTU_biom table
rownames(taxmat) <- rownames(OTU_biom)
#The column names are the taxonomic levels
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus","Species")
#Convert the taxonomy table to the format required for a phyloseq object
TAX<-tax_table(taxmat)

#Read in metadata
meta_csv<-read.csv('16S/Zoraptera_metadata.csv', row.names=1,header= TRUE)
#Convert the metadata table to the format required for a phyloseq object
SAMP<-sample_data(meta_csv)

#Create a phyloseq object by combining the OTU table, taxonomy table and sample metadata (could include a tree if we had one)
Z_physeq<-phyloseq(OTU_biom,TAX,SAMP)

#Rarefy the OTU table
rZOTU<-rarefy_even_depth(Z_physeq,rngseed=sample(1000,1))

#Change the sample names from the Plate IDs to the Colony IDs (unique to Zoraptera... because of our naming system... you may need to do something else)
sample_names(rZOTU) <- sample_data(rZOTU)$Colony


#%%%%%%%%%%%%%%%%%%%%%%   MAKE FANCY BAR GRAPHS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Find the total abundance of each phylum
summed_phyla<-rowMeans(data.frame(otu_table(tax_glom(rZOTU,taxrank="Phylum"))))
#Find the names of the different phylum
phyla_list<-data.frame(tax_table(tax_glom(rZOTU,taxrank="Phylum")))[,2]
#Sort the phylum names so that the first is the most abundant, the second the second most abundat... we will explicitly show the four most abundant
sorted_phyla_list<-phyla_list[order(-summed_phyla)]

#Prepare the OTU table in the correct format for the plotting program (this must be done with level 6 data including genera!)
mdf_prep <- prep_mdf(rZOTU)
#Tell the function to plot the five most abundant phyla (you could pick something different... )
color_objs_GP <- create_color_dfs(mdf_prep,selected_groups = c("Actinobacteria", "Proteobacteria", "Firmicutes", "Acidobacteria"  ),cvd=TRUE)
#Extract the OTU table and color choices (again, putting things in the right format for the function)
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf
cdf_GP <- color_reassign(cdf_GP,group_assignment = c("Actinobacteria", "Proteobacteria", "Firmicutes", "Acidobacteria"),color_assignment = c("micro_cvd_green", "micro_cvd_orange","micro_cvd_blue","micro_cvd_purple"))
#Define the legend
GP_legend <-custom_legend(mdf_GP, cdf_GP)

#Define the plot that you will be making
plot <- plot_microshades(mdf_GP, cdf_GP)
#Define all of the formatting aspects of the plot
plot_diff <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_grid(~fct_relevel(Type,'Adult','Juvenile','Environment'), scale="free_x", space = "free_x") +
  theme(axis.text.x = element_text(size= 6)) +
  theme(plot.margin = margin(6,20,6,6))
#Make the plot
plot_grid(plot_diff, GP_legend,  rel_widths = c(1, .25))


#Expand the number of Proteobacteria shown
new_groups <- extend_group(mdf_GP, cdf_GP, "Phylum", "Genus", "Proteobacteria", existing_palette = "micro_cvd_orange", new_palette = "micro_orange", n_add = 5)
GP_legend_new <-custom_legend(new_groups$mdf, new_groups$cdf)
#Expand the number of Actinobacteria shown
new_groups2 <- extend_group(new_groups$mdf, new_groups$cdf, "Phylum", "Genus", "Actinobacteria", existing_palette = "micro_cvd_green", new_palette = "micro_green", n_add = 5)
#Expand the number of Firmicutes shown
new_groups3 <- extend_group(new_groups2$mdf, new_groups2$cdf, "Phylum", "Genus", "Firmicutes", existing_palette = "micro_cvd_blue", new_palette = "micro_blue", n_add = 5)
#Expand the number of Acidobacteria shown
new_groups4 <- extend_group(new_groups3$mdf, new_groups3$cdf, "Phylum", "Genus", "Acidobacteria", existing_palette = "micro_cvd_purple", new_palette = "micro_purple", n_add = 5)

#Define your new legend with the expanded groups
GP_legend_new <-custom_legend(new_groups4$mdf, new_groups4$cdf,legend_text_size = 13)

#Define the plot that you will be making
plot2 <- plot_microshades(new_groups4$mdf, new_groups4$cdf)
#Define all the formatting aspects of the plot
plot_diff2 <- plot2 + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none",text=element_text(size=28))  +
  theme(axis.text.x = element_text(size= 10)) +
  facet_grid(~fct_relevel(Type,'Adult','Juvenile','Environment'), scale="free_x", space = "free_x") +
  theme(axis.text.x = element_text(size= 10)) +
  theme(plot.margin = margin(6,20,6,6))
#Make the plot
plot_grid(plot_diff2, GP_legend_new,  rel_widths = c(1, .25))


#Make a phyloseq object of just the Zoraptera samples
zZOTU<-prune_samples(sample_data(rZOTU)$Type %in% c("Adult","Juvenile"),rZOTU)

#Dominant Phyla by Mean Abundance
z_summed_phyla<-rowMeans(data.frame(otu_table(tax_glom(zZOTU,taxrank="Phylum")))/colSums(data.frame(otu_table(tax_glom(zZOTU,taxrank="Phylum")))))
zphyla_list<-data.frame(tax_table(tax_glom(zZOTU,taxrank="Phylum")))[,2]
names(z_summed_phyla)<-zphyla_list
z_summed_phyla<-sort(z_summed_phyla,decreasing=TRUE)
z_summed_phyla[1:6]*100

#Dominant Phyla by Prevalence
z_summed_phyla<-rowSums(sign(data.frame(otu_table(tax_glom(zZOTU,taxrank="Phylum")))))/35
zphyla_list<-data.frame(tax_table(tax_glom(zZOTU,taxrank="Phylum")))[,2]
names(z_summed_phyla)<-zphyla_list
z_summed_phyla<-sort(z_summed_phyla,decreasing=TRUE)
z_summed_phyla[1:6]*100
z_summed_phyla[which(names(z_summed_phyla)=="Cyanobacteria")]*100

#Dominant Genera by Mean Abundance
z_summed_genus<-rowMeans(data.frame(otu_table(tax_glom(zZOTU,taxrank="Genus")))/colSums(data.frame(otu_table(tax_glom(zZOTU,taxrank="Phylum")))))
zgenus_list<-data.frame(tax_table(tax_glom(zZOTU,taxrank="Genus")))[,6]
names(z_summed_genus)<-zgenus_list
z_summed_genus<-sort(z_summed_genus,decreasing=TRUE)
z_summed_genus[1:11]*100
#Abundance of unclassified_Rhizobiales 
z_summed_genus[which(names(z_summed_genus)=="unclassified_Rhizobiales")]*100
#Abundance of Aquicella 
z_summed_genus[which(names(z_summed_genus)=="Aquicella")]*100
#Abundance of Legionella
z_summed_genus[which(names(z_summed_genus)=="Legionella")]*100
#Abundance of unclassified_Acidimicrobiales 
z_summed_genus[which(names(z_summed_genus)=="unclassified_Acidimicrobiales")]*100
#Abundance of unclassified_Solirubrobacterales 
z_summed_genus[which(names(z_summed_genus)=="unclassified_Solirubrobacterales")]*100
#Abundance of unclassified_Acetobacteraceae 
z_summed_genus[which(names(z_summed_genus)=="unclassified_Acetobacteraceae")]*100
#Abundance of Conexibacter 
z_summed_genus[which(names(z_summed_genus)=="Conexibacter")]*100



#Dominant Genera by Prevalence
z_summed_genus<-rowSums(sign(data.frame(otu_table(tax_glom(zZOTU,taxrank="Genus")))))/35
zgenus_list<-data.frame(tax_table(tax_glom(zZOTU,taxrank="Genus")))[,6]
names(z_summed_genus)<-zgenus_list
z_summed_genus<-sort(z_summed_genus,decreasing=TRUE)
z_summed_genus[1:10]*100
#Prevalence of Wolbachia
z_summed_genus[which(names(z_summed_genus)=="Wolbachia")]*100
#Prevalence of Holosporaceae
z_summed_genus[which(names(z_summed_genus)=="unclassified_Holosporaceae")]*100
#Prevalence of Rickettsiella
z_summed_genus[which(names(z_summed_genus)=="Rickettsiella")]*100
#Prevalence of Streptacidiphilus-Streptomyces 
z_summed_genus[which(names(z_summed_genus)=="Streptacidiphilus-Streptomyces")]*100
#Prevalence of Aeromonas 
z_summed_genus[which(names(z_summed_genus)=="Aeromonas")]*100
#Prevalence of Acetobacter 
z_summed_genus[which(names(z_summed_genus)=="Acetobacter")]*100
#Prevalence of Actinocorallia 
z_summed_genus[which(names(z_summed_genus)=="Actinocorallia")]*100



