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
library('picante')
library('SYNCSA')
library('pracma')
library('ggpubr')
library('tidyverse')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Read in Data and Make ASV Table

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
tree_file<-read.tree('ITS/rooted_tree_out/tree.nwk')

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

#Plot the phylogenetic tree - does it look crazy? Some 'bad stuff' may be in there
#plot.phylo(phy_tree(Z_physeq))

#Find the sequences that Zymo was actually able to find in their database... other stuff is probably weird and may explain your funky tree
assigned_taxa<-read.csv('ITS/ASV_tax_assignments.csv')
assigned_taxa_seqs<-assigned_taxa[,1]

#Cut out the crazy (all the stuff Zymo couldn't find in their database)
Z_physeq = prune_taxa(assigned_taxa_seqs, Z_physeq)

#Plot the phylogenetic tree - does it look better?
#plot.phylo(phy_tree(Z_physeq))

#Remove samples below a certain number of reads (I found that the prune_taxa command was somehow adding back in the low read samples... so I'm removing them here...)
highreads<-which(colSums(otu_table(Z_physeq))>7000)
ZNC_physeq<-prune_samples(sample_names(Z_physeq) %in% c(sample_names(Z_physeq)[highreads]),Z_physeq)

#Rarefy the (now pruned) ASV table
rZOTU<-rarefy_even_depth(ZNC_physeq,rngseed=1)

#Find the type of sample for each sample
typeA<-sample_data(rZOTU)$Type
typeA[typeA=='Adult']<-'Insect'
typeA[typeA=='Juvenile']<-'Insect'
sample_data(rZOTU)$TypeA<-typeA
type<-sample_data(rZOTU)$TypeA

#Do you want to display significant differences as stars or p-values? ('star' = stars, 'pvalue' = pvalues)
pshow<-'star'


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Calculate a variety of different ALPHA diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Find the richness for each sample
richness<-colSums(sign(otu_table(rZOTU)))

#Find the shannon diversity for each sample
shannon<-diversity(otu_table(t(rZOTU)),index='shannon')

#Find the simpson's index for each sample
simpson<-diversity(otu_table(t(rZOTU)),index='simpson')

#Find Faith's pd for each sample
faithpd<-pd(t(otu_table(rZOTU)), tree_file, include.root=TRUE)
faiths<-faithpd[,1]

#The Rao's entropy package requires an ultrametric tree.... make the ultrametric tree and plot it to see how it's different
ultrametric_tree<-chronoMPL(multi2di(phy_tree(Z_physeq)))
#plot.phylo(ultrametric_tree)

#Calculate Faith's pd on the ultrametric tree
#ufaithpd<-pd(t(otu_table(rZOTU)), ultrametric_tree, include.root=TRUE)
#ufaiths<-ufaithpd[,1]

#Calculate Rao's entropy on the ultrametric tree
raopd<-raoD(comm=t(otu_table(rZOTU)),phy=ultrametric_tree)
raos<-raopd$Dkk

#Put all of your different diversity metrics into a dataframe
diversity_df<-data.frame(type,richness,shannon,simpson,faiths, raos)
#Switch the ordering of the types (to the order you want them presented in your plots)
diversity_df$type<-factor(diversity_df$type,levels=c('Insect','Environment'))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Test whether there are differences in ALPHA diversity between groups

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Richness
kw_richness<-kruskal.test(richness ~ type, data = diversity_df)
if (kw_richness$p.value<0.05){
  print('Alpha Richness')
  pw_richness<-pairwise.wilcox.test(diversity_df$richness, diversity_df$type,p.adjust.method = "BH")
  print(pw_richness)
}

#Shannon Entropy
kw_shannon<-kruskal.test(shannon ~ type, data = diversity_df)
if (kw_shannon$p.value<0.05){
  print('Alpha Shannon')
  pw_shannon<-pairwise.wilcox.test(diversity_df$shannon, diversity_df$type,p.adjust.method = "BH")
  print(pw_shannon)
}

#Simpson's Index
kw_simpson<-kruskal.test(simpson ~ type, data = diversity_df)
if (kw_simpson$p.value<0.05){
  print('Alpha Simpson')
  pw_simpson<-pairwise.wilcox.test(diversity_df$simpson, diversity_df$type,p.adjust.method = "BH")
  print(pw_simpson)
}

#Faith's PD
kw_faiths<-kruskal.test(faiths ~ type, data = diversity_df)
if (kw_faiths$p.value<0.05){
  print('Alpha Faiths')
  pw_faiths<-pairwise.wilcox.test(diversity_df$faiths, diversity_df$type,p.adjust.method = "BH")
  print(pw_faiths)
}

#Ultrametric Faith's PD
#kw_ufaiths<-kruskal.test(ufaiths ~ type, data = diversity_df)
#if (kw_ufaiths$p.value<0.05){
#  print('Alpha Ultrametric Faiths')
#  pw_ufaiths<-pairwise.wilcox.test(diversity_df$ufaiths, diversity_df$type,p.adjust.method = "BH")
#  print(pw_ufaiths)
#}

#Rao's Entropy
kw_raos<-kruskal.test(raos ~ type, data = diversity_df)
if (kw_raos$p.value<0.05){
  print('Alpha Raos')
  pw_raos<-pairwise.wilcox.test(diversity_df$raos, diversity_df$type,p.adjust.method = "BH")
  print(pw_raos)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#                 Make Violin plots of the various diversity metrics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Do you want to show p-values for non-significant pairwise comparisons? (yes or no)
p_value_ns<-'no'

##########################Richness Violin Plot#########################################################################################################################################

#Define violin plot
p_richness <- ggplot(diversity_df, aes(x=type, y=richness)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_richness<-p_richness+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_richness<-p_richness+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_richness<-p_richness+scale_fill_manual(values=c("wheat2", "brown"))
#Add boxplots inside the violins
p_richness<-p_richness+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_richness<-p_richness+annotate("text", size=6,x=2, y=80, label= "Kruskal-Wallis")
p_richness<-p_richness+annotate("text", size=6,x=2, y=65, label= paste("p = ",round(kw_richness$p.value,4)))

#If there are significant differences in richness between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_richness$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-190
  y_step<-10
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_richness$p.value))){
    for (j in 1:k){
      if (rownames(pw_richness$p.value)[k]!=colnames(pw_richness$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_richness$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_richness$p.value)[k])
          group2<-c(group2,colnames(pw_richness$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_richness$p.value[k,j])),5)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  #Change p-values to stars for plotting
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_richness<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_richness<-p_richness+stat_pvalue_manual(stat.test_richness,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_richness

##########################Shannon Diversity Violin Plot#########################################################################################################################################

#Define violin plot
p_shannon <- ggplot(diversity_df, aes(x=type, y=shannon)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_shannon<-p_shannon+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_shannon<-p_shannon+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_shannon<-p_shannon+scale_fill_manual(values=c("wheat2", "brown"))
#Add boxplots inside the violins
p_shannon<-p_shannon+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_shannon<-p_shannon+annotate("text", size=6,x=2, y=1.8, label= "Kruskal-Wallis")
p_shannon<-p_shannon+annotate("text", size=6,x=2, y=1.5, label= paste("p = ",round(kw_shannon$p.value,6)))

#If there are significant differences in shannon between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_shannon$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-4
  y_step<-0.5
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_shannon$p.value))){
    for (j in 1:k){
      if (rownames(pw_shannon$p.value)[k]!=colnames(pw_shannon$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_shannon$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_shannon$p.value)[k])
          group2<-c(group2,colnames(pw_shannon$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_shannon$p.value[k,j])),6)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  #Change p-values to stars for plotting
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_shannon<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_shannon<-p_shannon+stat_pvalue_manual(stat.test_shannon,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_shannon

##########################Simpson's Diversity Violin Plot#########################################################################################################################################

#Define violin plot
p_simpson <- ggplot(diversity_df, aes(x=type, y=simpson)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_simpson<-p_simpson+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_simpson<-p_simpson+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_simpson<-p_simpson+scale_fill_manual(values=c("wheat2", "brown"))
#Add boxplots inside the violins
p_simpson<-p_simpson+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_simpson<-p_simpson+annotate("text", size=6,x=2, y=0.6, label= "Kruskal-Wallis")
p_simpson<-p_simpson+annotate("text", size=6,x=2, y=0.55, label= paste("p = ",round(kw_simpson$p.value,6)))

#If there are significant differences in simpson between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_simpson$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-1.1
  y_step<-0.1
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_simpson$p.value))){
    for (j in 1:k){
      if (rownames(pw_simpson$p.value)[k]!=colnames(pw_simpson$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_simpson$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_simpson$p.value)[k])
          group2<-c(group2,colnames(pw_simpson$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_simpson$p.value[k,j])),6)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  #Change p-values to stars for plotting
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_simpson<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_simpson<-p_simpson+stat_pvalue_manual(stat.test_simpson,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_simpson

##########################Faith's PD Violin Plot#########################################################################################################################################

#Define violin plot
p_faiths <- ggplot(diversity_df, aes(x=type, y=faiths)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_faiths<-p_faiths+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_faiths<-p_faiths+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_faiths<-p_faiths+scale_fill_manual(values=c("wheat2", "brown"))
#Add boxplots inside the violins
p_faiths<-p_faiths+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_faiths<-p_faiths+annotate("text", size=6,x=2, y=17, label= "Kruskal-Wallis")
p_faiths<-p_faiths+annotate("text", size=6,x=2, y=14, label= paste("p = ",round(kw_faiths$p.value,6)))

#If there are significant differences in faiths between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_faiths$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-35
  y_step<-3
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_faiths$p.value))){
    for (j in 1:k){
      if (rownames(pw_faiths$p.value)[k]!=colnames(pw_faiths$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_faiths$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_faiths$p.value)[k])
          group2<-c(group2,colnames(pw_faiths$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_faiths$p.value[k,j])),6)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  #Change p-values to stars for plotting
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_faiths<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_faiths<-p_faiths+stat_pvalue_manual(stat.test_faiths,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_faiths



##########################Rao's Entropy Violin Plot#########################################################################################################################################

#Define violin plot
p_raos <- ggplot(diversity_df, aes(x=type, y=raos)) + geom_violin(aes(fill=type))
#Choose the size of font for the axes titles and labels
p_raos<-p_raos+theme(axis.title = element_text(size = 25))+theme(axis.text = element_text(size = 20))
#Choose the size of font for the legend title and lables
p_raos<-p_raos+theme(legend.title = element_text(size = 25))+theme(legend.text = element_text(size = 20))
#Choose the violin colors for each group
p_raos<-p_raos+scale_fill_manual(values=c("wheat2", "brown"))
#Add boxplots inside the violins
p_raos<-p_raos+geom_boxplot(aes(fill=type),width=0.1)
#Add the p-value for the Kruskal-Wallis test somewhere on your figure (you may have to change the x and y positions of this label)
p_raos<-p_raos+annotate("text", size=6,x=2, y=0.25, label= "Kruskal-Wallis")
p_raos<-p_raos+annotate("text", size=6,x=2, y=0.2, label= paste("p = ",round(kw_raos$p.value,6)))

#If there are significant differences in raos between groups, make a tibble of the pairwise p-values for plotting and add the brackets to the plot
if (kw_raos$p.value<0.05){
  #Names of the groups you're comparing
  group1<-c()
  group2<-c()
  #p-value for the pairwise comparisons
  p.adj<-c()
  #locations of the p-value brackets
  ypos<-c()
  new_y<-0.73
  y_step<-0.05
  #For each pairwise comparison...
  for (k in 1:length(rownames(pw_raos$p.value))){
    for (j in 1:k){
      if (rownames(pw_raos$p.value)[k]!=colnames(pw_raos$p.value)[j]){
        #If there is a significant difference or you want to also show non-significant differences...
        if (pw_raos$p.value[k,j]<0.05 || p_value_ns=='yes'){
          #Add an entry to your tibble include the names of the two groups being compared and the p-value for the comparison
          group1<-c(group1,rownames(pw_raos$p.value)[k])
          group2<-c(group2,colnames(pw_raos$p.value)[j])
          p.adj<-round(c(p.adj,as.numeric(pw_raos$p.value[k,j])),6)
          #Add the y-position of the bracket and bump the y-location of the bracket up one for next time
          ypos<-c(ypos,new_y)
          new_y<-new_y+y_step
        }
      }
    }
  }
  #Change p-values to stars for plotting
  pstar<-c()
  for (k in 1:length(p.adj)){
    if (p.adj[k]<=0.001){
      pstar<-c(pstar,'***')
    }
    else if (p.adj[k]<=0.01){
      pstar<-c(pstar,'**')
    }
    else if (p.adj[k]<=0.05){
      pstar<-c(pstar,'*')
    }
    else if (p.adj[k]<=0.1){
      pstar<-c(pstar,'.')
    }
    else{
      pstar<-c(pstar,'ns')
    }
  }
  if (pshow=='star'){
    pdisplay<-"{pstar}"
  }
  else{
    pdisplay<-"p = {p.adj}"
  }
  #Create your tibble (what is needed for the stat_pvalue_manual function)
  stat.test_raos<-as_tibble(data.frame(group1,group2,p.adj))
  #Add the pairwise comparisons to your plot
  p_raos<-p_raos+stat_pvalue_manual(stat.test_raos,label=pdisplay,y.position=ypos,size=5)
}
#Make your plot
p_raos

