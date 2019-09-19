library(ggplot2)
library(reshape2)
library(scales)
library(RColorBrewer)
require(grid)
library(gridExtra)
library(lattice)

otus_a <- read.csv('1_WE1_max_TSS.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm_a <- melt(cbind(otus_a, ind = rownames(otus_a)), id.vars = c('ind'))
p_a <- ggplot(datm_a, aes(variable, value, group = ind, color = ind)) + geom_line(size = 3)
p_a <- p_a + geom_point(size=6)+ scale_y_log10(breaks=c(0.01,0.1,1)) + scale_color_manual(values=c("#80cdc1","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#4d4d4d"))
p_a <- p_a + theme_bw() + xlab('') + ylab('')
p_a <- p_a + theme(panel.grid.major = element_line(colour = "grey40"))
p_a <- p_a + theme(panel.grid.minor = element_line(colour = "grey40"))
p_a <- p_a + theme(axis.title.y=element_text(angle=90, size=22, vjust=1.0))
p_a <- p_a + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p_a <- p_a + theme(axis.text.y=element_text(angle=0, size=22, vjust=0.5))
p_a <- p_a + theme(legend.text=element_text(size=22))
p_a <- p_a + theme(legend.title = element_text(size=22))
p_a <- p_a + theme(legend.key.size = unit(1.5, "cm"))
p_a <- p_a + theme(legend.title=element_blank())
png(file="log10_a_ISEC_WE1_TSS.png",width=1000,height=800,res=100)
p_a
dev.off()


otus_b <- read.csv('2_WE2_max_TSS.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm_b <- melt(cbind(otus_b, ind = rownames(otus_b)), id.vars = c('ind'))
p_b <- ggplot(datm_b, aes(variable, value, group = ind, color = ind)) + geom_line(size = 3)
p_b <- p_b + geom_point(size=6)+ scale_y_log10(breaks=c(0.01,0.1,1)) + scale_color_manual(values=c("#1f78b4","#e31a1c","#33a02c","#6a3d9a","#b15928"))
p_b <- p_b + theme_bw() + xlab('') + ylab('')
p_b <- p_b + theme(panel.grid.major = element_line(colour = "grey40"))
p_b <- p_b + theme(panel.grid.minor = element_line(colour = "grey40"))
p_b <- p_b + theme(axis.title.y=element_text(angle=90, size=22, vjust=1.0))
p_b <- p_b + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p_b <- p_b + theme(axis.text.y=element_text(angle=0, size=22, vjust=0.5))
p_b <- p_b + theme(legend.text=element_text(size=22))
p_b <- p_b + theme(legend.title = element_text(size=22))
p_b <- p_b + theme(legend.key.size = unit(1.5, "cm"))
p_b <- p_b + theme(legend.title=element_blank())
png(file="log10_b_ISEC_WE2_TSS.png",width=1000,height=800,res=100)
p_b
dev.off()


otus_c <- read.csv('3_WE3_max_TSS.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm_c <- melt(cbind(otus_c, ind = rownames(otus_c)), id.vars = c('ind'))
p_c <- ggplot(datm_c, aes(variable, value, group = ind, color = ind)) + geom_line(size = 3)
p_c <- p_c + geom_point(size=6)+ scale_y_log10() + scale_color_manual(values=c("#1f78b4","#e31a1c","#33a02c","#6a3d9a"))
p_c <- p_c + theme_bw() + xlab('') + ylab('')
p_c <- p_c + theme(panel.grid.major = element_line(colour = "grey40"))
p_c <- p_c + theme(panel.grid.minor = element_line(colour = "grey40"))
p_c <- p_c + theme(axis.title.y=element_text(angle=90, size=22, vjust=1.0))
p_c <- p_c + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p_c <- p_c + theme(axis.text.y=element_text(angle=0, size=22, vjust=0.5))
p_c <- p_c + theme(legend.text=element_text(size=22))
p_c <- p_c + theme(legend.title = element_text(size=22))
p_c <- p_c + theme(legend.key.size = unit(1.5, "cm"))
p_c <- p_c + theme(legend.title=element_blank())
png(file="log10_c_ISEC_WE3_TSS.png",width=1000,height=800,res=100)
p_c
dev.off()


otus_d <- read.csv('4_WE4_max_TSS.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm_d <- melt(cbind(otus_d, ind = rownames(otus_d)), id.vars = c('ind'))
p_d <- ggplot(datm_d, aes(variable, value, group = ind, color = ind)) + geom_line(size = 3)
p_d <- p_d + geom_point(size=6)+ scale_y_log10() + scale_color_manual(values=c("#1f78b4","#e31a1c","#33a02c","#6a3d9a"))
p_d <- p_d + theme_bw() + xlab('') + ylab('')
p_d <- p_d + theme(panel.grid.major = element_line(colour = "grey40"))
p_d <- p_d + theme(panel.grid.minor = element_line(colour = "grey40"))
p_d <- p_d + theme(axis.title.y=element_text(angle=90, size=22, vjust=1.0))
p_d <- p_d + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
p_d <- p_d + theme(axis.text.y=element_text(angle=0, size=22, vjust=0.5))
p_d <- p_d + theme(legend.text=element_text(size=22))
p_d <- p_d + theme(legend.title = element_text(size=22))
p_d <- p_d + theme(legend.key.size = unit(1.5, "cm"))
p_d <- p_d + theme(legend.title=element_blank())
png(file="log10_d_ISEC_WE4_TSS.png",width=1000,height=800,res=100)
p_d
dev.off()


otus_e <- read.csv('5_WE1_WE3_max_TSS.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm_e <- melt(cbind(otus_e, ind = rownames(otus_e)), id.vars = c('ind'))
p_e <- ggplot(datm_e, aes(variable, value, group = ind, color = ind)) + geom_line(size = 3)
p_e <- p_e + geom_point(size=6)+ scale_y_log10() + scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))
p_e <- p_e + theme_bw() + xlab('') + ylab('')
p_e <- p_e + theme(panel.grid.major = element_line(colour = "grey40"))
p_e <- p_e + theme(panel.grid.minor = element_line(colour = "grey40"))
p_e <- p_e + theme(axis.title.y=element_text(angle=90, size=22, vjust=1.0))
p_e <- p_e + theme(axis.text.x=element_text(size=22, angle = 45, hjust=1))
p_e <- p_e + theme(axis.text.y=element_text(angle=0, size=22, vjust=0.5))
p_e <- p_e + theme(legend.text=element_text(size=22))
p_e <- p_e + theme(legend.title = element_text(size=22))
p_e <- p_e + theme(legend.key.size = unit(1.5, "cm"))
p_e <- p_e + theme(legend.title=element_blank())
png(file="log10_e_ISEC_WE1_WE3_TSS.png",width=1000,height=800,res=100)
p_e
dev.off()


otus_f <- read.csv('6_WE2_WE4_max_TSS.csv',sep=',', header=T,check=F,row.names=1,comment='')
datm_f <- melt(cbind(otus_f, ind = rownames(otus_f)), id.vars = c('ind'))
p_f <- ggplot(datm_f, aes(variable, value, group = ind, color = ind)) + geom_line(size = 3)
p_f <- p_f + geom_point(size=6)+  scale_y_log10() + scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"))
p_f <- p_f + theme_bw() + xlab('') + ylab('')
p_f <- p_f + theme(panel.grid.major = element_line(colour = "grey40"))
p_f <- p_f + theme(panel.grid.minor = element_line(colour = "grey40"))
p_f <- p_f + theme(axis.title.y=element_text(angle=90, size=22, vjust=1.0))
p_f <- p_f + theme(axis.text.x=element_text(size=22, angle = 45, hjust=1))
p_f <- p_f + theme(axis.text.y=element_text(angle=0, size=22, vjust=0.5))
p_f <- p_f + theme(legend.text=element_text(size=22))
p_f <- p_f + theme(legend.title = element_text(size=22))
p_f <- p_f + theme(legend.key.size = unit(1.5, "cm"))
p_f <- p_f + theme(legend.title=element_blank())
png(file="log10_f_ISEC_WE2_WE4_TSS.png",width=1000,height=800,res=100)
p_f
dev.off()


png(file="log10_all_TSS.png",width=2500,height=2500,res=100)

# Get the gtables
gA <- ggplotGrob(p_a)
gB <- ggplotGrob(p_b)
gC <- ggplotGrob(p_c)
gD <- ggplotGrob(p_d)
gE <- ggplotGrob(p_e)
gF <- ggplotGrob(p_f) 

# Set the widths
gB$widths <- gA$widths
gC$widths <- gA$widths
gD$widths <- gA$widths
gE$widths <- gA$widths
gF$widths <- gA$widths

# Arrange the two charts.
# The legend boxes are centered
grid.newpage()
grid.arrange(gA, gB, gC, gD, gE, gF, ncol=2, nrow = 3)
dev.off()
