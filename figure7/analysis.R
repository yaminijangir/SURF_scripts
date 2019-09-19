library("microbiome")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("plyr"); packageVersion("plyr")
library(RColorBrewer)
library(tidyverse)
library(magrittr)
library(dplyr)
library(stringr)
library(readr)
library(readxl)
library(vegan)

theme_set(theme_bw(base_size=18))

pal = "Set1"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

# remove the (*number*) in the final.an.0.02.cons.taxonomy file by replacing all (*) with blank.

rm(list=ls())
pseq <- read_phyloseq("final.an.0.02.filter.shared", "final.an.0.02.cons.taxonomy", "final.metadata.csv", type = "mothur")
meta <- meta(pseq)
taxonomy <- tax_table(pseq)
otu.absolute <- abundances(pseq)
# creating a new pseq after removing rows/columns with OTUs contaning no abundances across the samples
PS = prune_taxa(taxa_sums(pseq) > 0, pseq)

# #merging PS based on SampleType
# PS.merged = merge_samples(PS, "SampleType")
# 
# #splitting the PS into two based on metadata "Run"
# PS.merged.subset.run1 <- subset_samples(PS.merged, Run == "1")
# PS.merged.subset.run2 <- subset_samples(PS.merged, Run == "2")

#create OTU table withjust the blank control sample
PS.blank1 <- subset_samples(PS, SampleType =='Blank1')
PS.blank2 <- subset_samples(PS, SampleType =='Blank2')


#find OTUs which are non zero in blank sample and writing them to a file
PS.blank1.filter = filter_taxa(PS.blank1, function(x) sum(x) > 0, TRUE)
PS.blank2.filter = filter_taxa(PS.blank2, function(x) sum(x) > 0, TRUE)
blank1.contaminants = tax_table(PS.blank1.filter)
blank2.contaminants = tax_table(PS.blank2.filter)
blank2.contaminants.df = as.data.frame(blank2.contaminants)
write.table(blank2.contaminants.df, file = "blank2_contaminants_df.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
blank1.contaminants.df = as.data.frame(blank1.contaminants)
write.table(blank1.contaminants.df, file = "blank1_contaminants_df.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")

#remove bad otus (blank contaminants) from run1 and run2
bad_OTUs_run1 = rownames(tax_table(PS.blank1.filter))
bad_OTUs_run2 = c("Otu00002",  "Otu00019",  "Otu00088",  "Otu00183",  "Otu0075",  "Otu00961", "Otu01182","Otu01312", "Otu02096", "Otu02122")
bad_OTUs_contaminant = union(bad_OTUs_run1,bad_OTUs_run2)
# good_Taxa <- setdiff(taxa_names(PS.merged), bad_OTUs_contaminant)
# PS.merged.removed.contaminants.r1.r2 <- prune_taxa(good_Taxa, PS.merged)
# PS.merged.removed.contaminants.r1.r2.blank.removed <- prune_samples(sample_names(PS.merged.removed.contaminants.r1.r2) != "Blank1" & sample_names(PS.merged.removed.contaminants.r1.r2) != "Blank2", PS.merged.removed.contaminants.r1.r2)
# 

# process the duhumidifier contaminants in run1 
# algorithm:
# if OTU is present only in dehumidifier and not in DUSEL 3A - remove (most possibly origin is different) --> remove OTU's not present in DUSEL 3A
# if OTU is present in dehumidifier and DUSEL 3A - keep (origin could be dehumidifier or DUSEL 3A)
# if OTU is only present in DUSEL 3A - keep (orgin DUSEL 3A)

#create OTU table with just the dehumidifier and DUSEL3A sample
PS.dehumidifier <- subset_samples(PS, SampleType =='Dehumidifier')
PS.dusel3A <- subset_samples(PS, SampleType =='DUSEL 3A')

#find OTUs which are non zero in blank sample and writing them to a file
PS.dehumidifier.filter = filter_taxa(PS.dehumidifier, function(x) sum(x) > 0, TRUE)
PS.dusel3A.filter = filter_taxa(PS.dusel3A, function(x) sum(x) == 0, TRUE)
PS.dehumidifier.filter.tax = tax_table(PS.dehumidifier.filter)
PS.dusel3A.filter.tax = tax_table(PS.dusel3A.filter)
# PS.dehumidifier.filter.tax.df = as.data.frame(PS.dehumidifier.filter.tax)
# write.table(blank2.contaminants.df, file = "PS.dehumidifier.filter.tax.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")
# PS.dusel3A.filter.tax.df = as.data.frame(PS.dusel3A.filter.tax)
# write.table(blank1.contaminants.df, file = "PS.dusel3A.filter.tax.csv",row.names=FALSE, na="",col.names=TRUE, sep=",")

nonzero.OTUs.dehumidifier = rownames(tax_table(PS.dehumidifier.filter))
zero.OTUs.dusel3A = rownames(tax_table(PS.dusel3A.filter))

remove_OTUs = union(bad_OTUs_contaminant,zero.OTUs.dusel3A)
good_Taxa <- setdiff(taxa_names(PS), remove_OTUs)

PS.removed.OTUs.contaminants.r1.r2.nonzero.dusel3A <- prune_taxa(good_Taxa, PS)

PS.final = subset_samples(PS.removed.OTUs.contaminants.r1.r2.nonzero.dusel3A, SampleType != "Blank1" & SampleType !="Blank2")

PS.final.run1 <- subset_samples(PS.final, SampleType != "Dehumidifier" & SampleType !="ISEC Plank" & SampleType != "ISEC OC" & Location != "Lab Enrichment" & Location != "Lab Enrichment Control")
PS.final.run2 <- subset_samples(PS.final, SampleType != "Dehumidifier" & Potential !="none" & Potential !="OC" & Location != "ISEC")
PS.final.run1.run2 <- subset_samples(PS.final, Potential != "OC" & Potential != "none" & SampleType != "Lab Enrichment WE3")



z = plot_richness(PS.final.run1.run2, x="Location", color= "Potential",  measures=c("Observed","Chao1", "Shannon","Simpson"))
z <- z + geom_point(size=4, alpha=0.7) 
z<- z + scale_color_manual(values=c("red", "darkgreen","blue","purple","saddlebrown"))

# OTU relative abundances
# xt <- transform(x, 'compositional')

PS.final.run1.merged = merge_samples(PS.final.run1, "SampleType")
PS.final.run1.merged.relative <- transform(PS.final.run1.merged, "compositional", scale = 1)      
# PS.final.run1.relative.filter <- filter_taxa(PS.final.run1.merged.relative, function(x) var(x) > 1e-06, TRUE)

PS.final.run1.merged.relative.phylum <- aggregate_taxa(PS.final.run1.merged.relative , "Phylum") 
q = plot_bar(PS.final.run1.merged.relative.phylum, fill="Phylum") + facet_wrap(~Phylum, scales = "free", nrow =1)
PS.final.run1.merged.relative.proteobacteria = subset_taxa(PS.final.run1.merged.relative, Phylum == "Proteobacteria")
r = plot_bar(PS.final.run1.merged.relative.proteobacteria,fill="Order") + facet_wrap(~Order, scales = "free", nrow=1)




PS.final.run1.run2.merged = merge_samples(PS.final.run1.run2, "SampleType")

PS.final.run1.run2.merged.relative <- transform(PS.final.run1.run2.merged, "compositional", scale = 1)      
PS.final.run1.run2.merged.relative.filter <- filter_taxa(PS.final.run1.run2.merged.relative, function(x) var(x) > 1e-06, TRUE)

PS.final.run1.run2.merged.relative.phylum.filter <- aggregate_taxa(PS.final.run1.run2.merged.relative.filter , "Phylum") 
s = plot_bar(PS.final.run1.run2.merged.relative.phylum.filter, fill="Phylum") + facet_wrap(~Phylum, scales = "free", nrow=1)
PS.final.run1.run2.merged.relative.proteobacteria.filter = subset_taxa(PS.final.run1.run2.merged.relative.filter, Phylum == "Proteobacteria")
t = plot_bar(PS.final.run1.run2.merged.relative.proteobacteria.filter,fill="Order") + facet_wrap(~Order, scales = "free", nrow=1)


# PS.final.run1.run2.merged.relative.family.glom <- tax_glom(PS.final.run1.run2.merged.relative.filter , "Family") 
# s = plot_bar(PS.final.run1.run2.merged.relative.family.glom , fill="P") + facet_wrap(~Phylum, scales = "free")

PS.final.run1.run2.merged.relative.family.aggregate <- aggregate_taxa(PS.final.run1.run2.merged.relative.filter , "Family") 

# PS.final.run1.merged.relative.Bacilli = subset_taxa(PS.final.run1.merged.relative, Family == "Bacilli")
# PS.final.run1.merged.relative.comamonadaceae = subset_taxa(PS.final.run1.merged.relative, Family == "Comamonadaceae")
# 
# 
# PS.final.run1.merged.relative.family <- aggregate_taxa(PS.final.run1.merged.relative , "Family")



# # 
# # #### Glom plot box plot of phylum in run1
# # 
# # PS.final.run1.merged.relative
# # 
# # Extract abundance matrix from the phyloseq object
# otu_data_final = as(otu_table(PS.final.run1.run2.merged), "matrix")
# # transpose if necessary
# if(taxa_are_rows(PS.final.run1.run2.merged)){otu_data_final  <- t(otu_data_final)}
# # Coerce to data.frame
# otu_data_final.df = as.data.frame(otu_data_final)
# write.table(otu_data_final.df, file = "otu_data_final.shared", na="",col.names=TRUE, sep="\t")
# 
# tax_data_final = tax_table(PS.final.run1.run2.merged)
# tax_data_final.df = as.data.frame(tax_data_final)
# write.table(tax_data_final.df, file = "tax_data_final.taxonomy", na="",col.names=TRUE, sep="\t")



dist_methods = c("jsd","jaccard","chao","bray","w","gower")
plist <- vector("list", length(dist_methods))
names(plist) = dist_methods
for( i in dist_methods ){
  # Calculate distance matrix
  iDist <- distance(PS.final.run1.run2, method=i)
  # Calculate ordination
  iMDS  <- ordinate(PS.final.run1.run2, "MDS", distance=iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(PS.final.run1.run2, iMDS, color="Location", shape ="Potential")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  # Save the graphic to file.
  plist[[i]] = p
}


df = ldply(plist, function(x) x$data)
names(df)[1] <- "distance"
p = ggplot(df, aes(Axis.1, Axis.2, color=Location, shape=Potential))
p = p + geom_point(size=4, alpha=0.7) 
p = p + facet_wrap(~distance, scales="free")
p <- p + scale_color_manual(values=c("red", "darkgreen","blue"))
p <- p + scale_shape_manual(values=c(15,16,17,18,9))
p = p + ggtitle("Beta diversity measure via PCoA Ordination")
p <- p +  theme(plot.title = element_text(hjust = 0.5), base_size = 18) 


png(filename="betadiversity.png", 
    type="cairo",
    units="in", 
    width=15, 
    height=12, 
    pointsize=12, 
    res=200)
p
dev.off()

png(filename="alphadiversity.png", 
    type="cairo",
    units="in", 
    width=15, 
    height=12, 
    pointsize=12, 
    res=200)
z

dev.off()


alpha_diversity <- estimate_richness(PS.final.run1.run2)
H <- alpha_diversity$Shannon
S1 <- alpha_diversity$Observed
S <- log(S1)
evenness <- H/S
alpha_diversity$Evenness = evenness
alpha_diversity
df <- data.frame(alpha_diversity, sample_data(PS.final.run1.run2))
df
df2 <- tidyr::gather(df, key = "Measure", value = "Value", Observed, Shannon, Evenness)
df2
e <- ggplot(data = df2, aes(x = Location , y = Value, colour=Potential, shape=DateSampled)) + facet_wrap(~Measure, scale = "free") +
geom_point(size=4, alpha=0.7)  + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Alpha Diversity Measure")
e <- e + scale_color_manual(values=c("red", "darkgreen","blue","purple","saddlebrown"))

#----t -test for alphadiversity---
d = sample_data(PS.final.run1.run2)
result<- estimate_richness(PS.final.run1.run2, measures = 'Shannon')
DUSEL3A = result[d[,'Location'] == 'DUSEL 3A',]
ISEC = result[d[,'Location'] == 'ISEC',]
Lab = result[d[,'Location'] == 'Lab Enrichment',]
t.test(DUSEL3A, ISEC)
t.test(DUSEL3A, Lab)



png(filename="alphadiversityeithevenness.png", 
    type="cairo",
    units="in", 
    width=15, 
    height=12, 
    pointsize=12, 
    res=200)
e

dev.off()


##-------statistical tests ----

PS.final.run1.run2_bray <- phyloseq::distance(PS.final.run1.run2, method = "bray")
sampledf <- data.frame(sample_data(PS.final.run1.run2))
adonis(PS.final.run1.run2_bray ~ Location, data = sampledf)
beta <- betadisper(PS.final.run1.run2_bray, sampledf$Location)
permutest(beta)



PS.final.run1.run2_jaccard <- phyloseq::distance(PS.final.run1.run2, method = "jaccard")
sampledf <- data.frame(sample_data(PS.final.run1.run2))
adonis(PS.final.run1.run2_jaccard ~ Location, data = sampledf)
beta <- betadisper(PS.final.run1.run2_jaccard, sampledf$Location)
permutest(beta)
