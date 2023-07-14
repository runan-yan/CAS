#Load packages
library(ape)
library(vegan)
library(ggplot2)
library(phyloseq)
library(cowplot)
library(tidyr)
library(dplyr)
library(ggpubr)
library(compositions)
library(zCompositions)
library(ALDEx2) #not installed
library(viridis)
library(dendextend)
library(readxl)
library(gplots)
library(BiocParallel)
library(pairwiseAdonis)
library(SpadeR)
library(psych)
library(plotly)
library(htmlwidgets)
library(glmnet)
library(Matrix)
set.seed(336)

#Import OTU table - 16S
otus_16s<-as.data.frame(import_mothur(mothur_shared_file = '16s.shared'))

#Import taxonomy table -16s
taxon_16s <- as.data.frame(import_mothur(mothur_constaxonomy_file = '16s.taxonomy'))
colnames(taxon_16s) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_16s =tax_table(as.matrix(taxon_16s))

#Import metadata -16s
metadata_16s <-read.csv("metadata.csv", header=TRUE, row.names=1)
META_16s = sample_data(metadata_16s)

#### OBTAIN READS INFORMATION AFTER MOTHUR ####
#These data are used for the reads analysis after OTU picking in the file Reads.R
reads_16s<-as.data.frame(colSums(otus_16s))
reads_16s$Sample<-row.names(reads_16s)
reads_16s$Origin<-rep("16s",40)  ## 
sum(reads)

#### RAREFACTION CURVES ####
phyloseq_rare_16s = phyloseq(otu_table(otus_16s, taxa_are_rows = TRUE), TAX_16s, META_16s)
## add ggrare function
rare_16s <- ggrare(phyloseq_rare_16s, step = 1000, se=TRUE, color="combine")


##ALPHA DIVERSITY
#For 16s rRNA data
phyloseq_16s = phyloseq(otu_table(otus_16s, taxa_are_rows = TRUE), TAX_16s, META_16s)
alpha16s_r <-estimate_richness(phyloseq_16s, measures=c("Shannon", "InvSimpson", "Chao1"))
estimate_richness(phyloseq_16s, split= TRUE, measures=c("Chao1", "Shannon", "InvSimpson"))

#write.csv(alpha16s_r, "alpha_16s_r.csv")

#Open the file created in Excel and add metadata information 

#Save and open in R
alpha_16s_r <- read.csv("alpha_16s_r.csv", sep = ",", header = T, row.names = 1)

#Pairwise.t.test for alpha diversity using Shannon and Inverse Simpson indices
pairwise.t.test(alpha_16s_r$Shannon, alpha_16s_r$combine, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_16s_r$InvSimpson, alpha_16s_r$combine, p.adjust.method = "bonferroni")

#### ESTIMATE % DISCOVERED DIVERSITY ####

##Calculate estimated richness for each sample.

## Original script adopted from "https://cran.r-project.org/web/packages/SpadeR/SpadeR.pdf"

## SpadeR::ChaoSpecies function is used to estimate richness.


richness_estimate = function(otu,...) {
  
  options(warn = -1)
  
  b = data.frame(matrix(nrow=as.matrix(dim(otu))[2,], ncol=3))
  
  colnames(b) <- c("Chao1 Estimates","Observed OTUs", "%Covered Species")
  
  for (i in 1:as.matrix(dim(otu))[2,]) {
    
    a =SpadeR::ChaoSpecies(otu[,i], datatype="abundance", k=10, conf=0.95)
    
    b[i,1]= as.numeric(a$Species_table[3,1])
    
    b[i,2]= apply(as.data.frame(a$Basic_data_information),2,as.numeric)[2,2]
    
    b[i,3]= (b[i,2]/b[i,1])*100
    
    rownames(b) <- colnames(otu) }
  
  print(b)
  
}

spadeR_16s_estimate <- richness_estimate(otus_16s)
spadeR_16s_estimate$Year<-rep("16s", 40)

#Merge with metadata
spadeR_16s<-bind_cols(spadeR_16s_estimate, metadata_16s)

#Statistical analysis
#Summary statistics
describeBy(spadeR_16s$`%Covered Species`, group=spadeR_16s$combine, mat = TRUE) #By meat category

#t-test

t.test_16s<-t.test(`%Covered Species`~ combine, data=spadeR_16s) # the comparison has to be made between two levels only


#### COMPOSITIONAL ANALYSIS OF MICROBIOME ####
#Based on Microbiome Analysis in R. Chap 10.
#At the family level

### Collapse OTU table to genus level
# Make phyloseq object
OTU_16s <- otu_table(otus_16s, taxa_are_rows = TRUE)
phyloseq_16s = phyloseq(OTU_16s, TAX_16s, META_16s)
TREE_16s = rtree(ntaxa(phyloseq_16s), rooted=TRUE, tip.label = taxa_names(phyloseq_16s))  
phyloseq16s <- phyloseq(OTU_16s, TAX_16s, TREE_16s, META_16s)

phyloseq16s_genus<-tax_glom(phyloseq16s,taxrank = "Genus")

#Extract OTU level at the genus level from Phyloseq object
otus_16s_genus<-as.data.frame(otu_table(phyloseq16s_genus))

#Step 1: Convert OTU table to appropriate format
#Following step requires samples on rows and OTUs in columns

head(t(otus_16s_genus)) #check that data in the correct format: samples on rows and OTUs in columns

#Step 2: Replace zero values before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
# Funciton: https://rdrr.io/cran/zCompositions/src/R/cmultRepl.R

otu.n0_16s<-t(cmultRepl(t(otus_16s_genus), label=0, method="CZM", output="p-counts")) #141 corrected values

head(otu.n0_16s)

#Step 3: Convert data to proportions
otu.n0_16s_prop<-apply(otu.n0_16s, 2, function(x) {x/sum(x)})

head(otu.n0_16s_prop)

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
otu.n0_16s_prop_f<-otu.n0_16s[apply(otu.n0_16s_prop, 1, min) > 0.0000001, ]

head(otu.n0_16s_prop_f)

#Step 5: perform CLR transformation
otu.n0.clr_16s<-t(apply(otu.n0_16s_prop_f, 2, function(x){log(x)-mean(log(x))}))
head(otu.n0.clr_16s)


#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16s<-prcomp(otu.n0.clr_16s)

png("Screeplot - PCA - by genus.png", width = 600, units = "px", height = 600)
par(mar=c(1,1,1,1))
par(mfrow=c(1,1))
screeplot(pc.clr_16s, type='lines', main= "16s - category")

dev.off()

#Calculate total variance of the data
# Function: https://rdrr.io/cran/mgm/src/R/mvar.R

# may need to re-install package "glmnet"
install.packages("glmnet")
install.packages('glmnet', dependencies=TRUE)

library("glmnet")

ncol(otu.n0.clr_16s) #obtain number of columns
mvar.clr_16s<-mvar(otu.n0.clr_16s, lags=1, type = rep(c("g"),1066))



#Display results - 16s
row_16s <-rownames(otu.n0.clr_16s)
pc_out_16s<-as.data.frame(pc.clr_16s$x[,1:3]) #Get PC1 and PC2 
pc_out_meta_16s<-as.data.frame(bind_cols(pc_out_16s,metadata_16s)) 
row.names(pc_out_meta_16s)<-row_16s 
pc_out_meta_16s$combine<-as.factor(pc_out_meta_16s$combine)

# Make PCA plot - First 2 axis 

PCA_16s <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=combine, shape=Packaging))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits =1), "%", sep="")) +
  ggtitle("16s")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s

#### PERMANOVA ####
#Using Aitchinson distances calculated previously for the dendrogram

#16s
dist_16s<-dist(otu.n0.clr_16s, method='euclidean') #Distance matrix based on Aitchinson simplex


permanova_16s_combine<-pairwise.adonis(dist_16s, factors=metadata_16s$combine, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_16s_combine$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')

#### MICROBIOTA AND MYCOBIOTA COMPOSITION PLOTS ####
#Note: used compositional approach to transform the sample counts to compositions. 

#Transform sample counts into compositions

otu.n0.acomp_16s<-as.data.frame(acomp(t(otu.n0_16s)), total=1)
rowSums(otu.n0.acomp_16s)

#Make Phyloseq object at genus/family level
OTU_16s_genus_prop <- otu_table(otu.n0.acomp_16s, taxa_are_rows = FALSE)
phyloseq_16s_genus_prop = phyloseq(OTU_16s_genus_prop, TAX_16s, META_16s)
TREE_16s_genus_prop = rtree(ntaxa(phyloseq_16s_genus_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_16s_genus_prop))  
phyloseq16s_genus_prop <- phyloseq(OTU_16s_genus_prop, TAX_16s, TREE_16s, META_16s)

OTU_16s_family_prop <- otu_table(otu.n0.acomp_16s, taxa_are_rows = FALSE)
phyloseq_16s_family_prop = phyloseq(OTU_16s_family_prop, TAX_16s, META_16s)
TREE_16s_family_prop = rtree(ntaxa(phyloseq_16s_family_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_16s_family_prop))  
phyloseq16s_family_prop <- phyloseq(OTU_16s_family_prop, TAX_16s, TREE_16s, META_16s)


otu_16s_otu_prop <- otu_table(otu.n0.acomp_16s, taxa_are_rows = FALSE)
phyloseq_16s_otu_prop = phyloseq(otu_16s_otu_prop, TAX_16s, META_16s)
TREE_16s_otu_prop = rtree(ntaxa(phyloseq_16s_otu_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_16s_otu_prop))  
phyloseq16s_otu_prop <- phyloseq(otu_16s_otu_prop, TAX_16s, TREE_16s, META_16s)



#Make long format table from Phyloseq object


otu_16s_obj1 <- phyloseq16s_otu_prop %>%
  transform_sample_counts(function(x){x*100}) %>%
  psmelt() %>%
  arrange(desc(Abundance))

genus_16s<- phyloseq16s_genus_prop %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

family_16s <- phyloseq16s_family_prop %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order


write.csv(otu_16s_obj1, "otu_16s_obj1.csv")
write.csv(family_16s, "family_16s.csv")
write.csv(genus_16s, "genus_16s.csv")

#Filter table to obtain only OTUs with over 0.5% in at least one sample
family_over05abund<- filter(family_16s, Abundance >0.5)
genus_over05abund <- filter(genus_16s, Abundance >0.5)
otu_over05abund <-filter(otu_16s_obj1,Abundance > 0.5)

#Calculate total abundance of each ASV by meat category
family_total<-as.data.frame(xtabs(Abundance ~ Sample, family_over05abund))
genus_total<-as.data.frame(xtabs(Abundance ~ Sample, genus_over05abund))
otu_total<-as.data.frame(xtabs(Abundance ~ Sample, otu_over05abund))

#Calculate "Others" total
family_total$Diff<-100-family_total$Freq
genus_total$Diff<-100-genus_total$Freq
otu_total$Diff<-100-otu_total$Freq

write.csv(otu_total, "otu_diff.csv")

##Replace abundance <0.5% with others in the local computer
otu_over05abund_final_obj1<-read.csv('otu_16s_obj1_modified.csv', header = TRUE, row.names = 1)

#Calculate % RA by genera per sample
genera_perc<-otu_over05abund_final_obj1%>%
  group_by(Sample, Genus)%>%
  tally(Abundance)

#Use ggplot2 to make the stacked barplot 

casPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#C99E8C","#8273E6","#59DBAF","#9DB366",
                "#E6CE9A","#66A7B3","#C71A6B","#34B5F7","#ED8ED4","#FF0000","#C5E6ED",
                "#C8C7C5","#98A594","#DBCCB1","#8A77A5","#DBE196","#467897","#E7CD79","#800020",
                "#002FA7","#BE98AA","#3E3F4C","#01847F","#F9D2E4","#FF770F","#80D1C8","#F8F5D6",
                "#FAEAD3","#492D22","#D8C7B5","#003153","#E5DDD7","#FFE76F","#95C194","#B1D1CF","#F57C69",
                "#665657","#A2E2C7","#FCB1AB","#779DE9","#F8C1B7","#FC9DA9","#BFDB36","#DB5A6C","#1A507E",
                "#805475","#BACAC7","#8AADA6","#BDD1B6","#45C4BE","#F4E6DD","#C2AFA1","#E3E1D2","#FED9B3","#DC6B88",
                "#0C03FF","#FFDDAB","#958CDD","#5A4590","#ECB42B","#E2D3B6","#6CA08B","#C98286","#E2DDB7")


#Genus  fill color
barplot_obj1.1<- ggplot(otu_over05abund_final_obj1, aes(x=reorder(Sample,SampleOrder), y = Abundance , fill = Family))  + 
  geom_bar(stat = "identity",color="black") + 
  #geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=3)+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=2)) +
  xlab("Chicken meat category") +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="right")+
  ggtitle("Microbiota composition of different chicken meat category", subtitle = "Family level")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  #Scale_fill_discrete(name="Family")
  scale_fill_manual(values=casPalette)

barplot_obj1.2<- ggplot(otu_over05abund_final_obj1, aes(x=reorder(Sample,SampleOrder), y = Abundance , fill = Genus))  + 
  geom_bar(stat = "identity",color="black") + 
  #geom_text(aes(label=Genus), position = position_stack(vjust=0.5), size=3)+
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=2)) +
  xlab("Chicken meat category") +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="right")+
  ggtitle("Microbiota composition of different chicken meat category", subtitle = "Genus level")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  #Scale_fill_discrete(name="Genus")
  scale_fill_manual(values=casPalette) +
  facet_wrap(~combine, strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")


### Core microbiota analysis ####
library(microbiome)

## subset data collection for BRRG, BRVA and WHRG
brrg_obj2 <- subset(otu_16s_obj1,combine=="BRRG")
brva_obj2 <-subset(otu_16s_obj1,combine=="BRVA")
whrg_obj2 <-subset(otu_16s_obj1,combine=="WHRG")
brma_obj2 <-subset(otu_16s_obj1,combine=="BRMA")

# make phyloseq project for each category
##subset phyloseq
phyloseq_brrg <-subset_samples(phyloseq16s_otu_prop,combine=="BRRG")
phyloseq_brva <-subset_samples(phyloseq16s_otu_prop,combine=="BRVA")
phyloseq_whrg <-subset_samples(phyloseq16s_otu_prop,combine=="WHRG")
phyloseq_brma <-subset_samples(phyloseq16s_otu_prop,combine=="BRMA")

## CORE TAXA definition: present in all 10 samples, relative abundance >0.5%
brrg.core<- core_members(phyloseq_brrg, detection = 0.0005, prevalence = 0.99)
brva.core<- core_members(phyloseq_brva, detection = 0.0005, prevalence = 0.99)
whrg.core<- core_members(phyloseq_whrg, detection = 0.0005, prevalence = 0.99)
brma.core<- core_members(phyloseq_brma, detection = 0.0005, prevalence = 0.99)

## Sum abundance of core microbiota
brrg.core.abundance <- sample_sums(core(phyloseq_brrg, detection = .0005, prevalence = .99))
brva.core.abundance <- sample_sums(core(phyloseq_brva, detection = .0005, prevalence = .99))
whrg.core.abundance <- sample_sums(core(phyloseq_whrg, detection = .0005, prevalence = .99))
brma.core.abundance <- sample_sums(core(phyloseq_brma, detection = .0005, prevalence = .99))

# Use the microbiome function add_besthit to get taxonomic identities of ASVs.
phyloseq_brrg2<- add_besthit(phyloseq_brrg)
phyloseq_brva2<- add_besthit(phyloseq_brva)
phyloseq_whrg2<- add_besthit(phyloseq_whrg)
phyloseq_brma2<- add_besthit(phyloseq_brma)

# Check 
taxa_names(phyloseq_brrg2)[1:10]
taxa_names(phyloseq_brva2)[1:10]
taxa_names(phyloseq_whrg2)[1:10]

##No add the best taxonomic classification to the core microbiota identification

brrg.core.taxa<- core_members(phyloseq_brrg2, detection = 0.0005, prevalence = 0.99)
brva.core.taxa<- core_members(phyloseq_brva2, detection = 0.0005, prevalence = 0.99)
whrg.core.taxa<- core_members(phyloseq_whrg2, detection = 0.0005, prevalence = 0.99)
brma.core.taxa<- core_members(phyloseq_brma2, detection = 0.0005, prevalence = 0.99)


# core visualization
# genus level
brrg.core.taxa.gen <- aggregate_taxa(phyloseq_brrg, "Genus")
brva.core.taxa.gen <- aggregate_taxa(phyloseq_brva, "Genus")
whrg.core.taxa.gen <- aggregate_taxa(phyloseq_whrg, "Genus")
brma.core.taxa.gen <- aggregate_taxa(phyloseq_whrg, "Genus")


## Manually prepare core taxa table 
core_taxa <-read.csv("core.csv",header = TRUE, row.names = 1)
core_meta <-read.csv("metadat_core.csv",header = TRUE, row.names = 1)

## Prepare phyloseq project for core taxa
taxon2 <-as.matrix(taxon_16s)
phyloseq_core<-phyloseq(otu_table(core_taxa, taxa_are_rows = FALSE), tax_table(taxon2), sample_data(core_meta))


core_CoDa <- phyloseq_core %>%  
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

write.csv(core_CoDa,"core_coda.csv")

#Calculate total abundance of each ASV by meat category
core_total<-as.data.frame(xtabs(Abundance ~ Sample, core_CoDa))

#Calculate "Others" total
core_total$Diff<-100-core_total$Freq

write.csv(core_total, "core_diff.csv")

##manually add non core abundance to the file
core_CoDa_modified <-read.csv("core_coda_modified.csv")

barplot_core1<- ggplot(core_CoDa_modified, aes(x=reorder(Sample,SampleOrder), y = Abundance , fill = Genus))  + 
  geom_bar(stat = "identity",color="black") + 
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(legend.text = element_text(face = "italic")) +
  theme(axis.title=element_text(size=15),axis.text = element_text(color='black', size=10)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Chicken meat category") +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position="right")+
  ggtitle("Core microbiota", subtitle = "Genus level")+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=2)) +
  #scale_fill_discrete(name="Genus")
  scale_fill_manual(values=casPalette) +
  facet_wrap(~combine, strip.position = "bottom", scales = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")
