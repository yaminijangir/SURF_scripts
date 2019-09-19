otus <- read.csv('SURF_ISECinLab_with_TSS_WE1_WE2_WE4.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm <- melt(cbind(otus, ind = rownames(otus)), id.vars = c('ind'))
write.table(datm, "mydata_WE1_WE2_WE4.csv", sep=",")

#Appended location criteria to the datm file. and stored the file as datm_SURF_ISECinLAb_no_TSS.csv

datm_new <- read.table('datm_SURF_ISECinLab_with_TSS_WE1_WE2_WE4.csv',sep=',', header=T,check=F,comment='')

datm_new$variable <- as.character(datm_new$variable)
datm_new$variable <- factor(datm_new$variable, levels=unique(datm_new$variable))
datm_new$location <- as.character(datm_new$location)
datm_new$location <- factor(datm_new$location, levels=unique(datm_new$location))


plog <- ggplot(datm_new, aes(variable, ind, color = factor(location)))+labs(size="Normalised Abundances", colour="Location")
plog <- plog + geom_point(aes(size = value))+ scale_color_manual(values=c("red", "darkgreen"))
plog <- plog + scale_size(trans="log10", range = c(-2, 7))
plog <- plog + theme_bw() + xlab('') + ylab('')
flevels <- levels (datm_new$ind)
flevels <- rev(flevels)
plog <- plog + scale_y_discrete(limits=flevels)
plog <- plog + theme(axis.text.x=element_text(size=14, angle = 45, hjust=1, colour="black"))
plog <- plog + theme(axis.text.y=element_text(angle=0, size=14, vjust=0.5, colour="black"))
plog <- plog + theme(legend.text=element_text(size=14, colour="black"))
plog <- plog + theme(legend.title = element_text(size=14, colour="black"))
plog <- plog + theme(legend.key.size = unit(1, "cm"), legend.key = element_rect(colour = "black"))
tiff(file="ISECinLAB_TSS_WE1_WE2_WE4_3.tiff",width=1000,height=800)
plog
dev.off()
