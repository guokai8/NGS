######density plot for methylation 
library(circlize)
###prepare the data for hyper and hypo
###both these two data set should have three columns:chr,start,end
###chr   start     end
## chr1  919910  919910
## chr1  984384  984384
###
bfile<-list(hyper.hypo)
##now make the circle plot
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
###if you want include all chromosomes use the code below
#circos.initializeWithIdeogram()
#plot the methylation CpG site
circos.genomicRainfall(bfile, pch = 16, cex = 0.4, col = c("red","blue"))
##plot the density for hyper
circos.genomicDensity(hyper, col = c("red"), track.height = 0.1)
##plot the density for hypo
circos.genomicDensity(hypo, col = c("blue"), track.height = 0.1)
dev.print(pdf,file="cir_density.pdf")

#######Difference methylation change plot
##first got all chromosome information
chrom<-read.delim("hg19.chrom.sizes.txt",sep="\t",header=F)
chrom.length<-chrom$V2
names(chrom.length)<-chrom$V1
###build GRanges
require(GenomicRanges) 
require(ggbio)
mychr<-GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, width = chrom.length))
seqlevels(mychr)<-names(chrom.length)
seqlengths(mychr) = chrom.length
###prepare the hyper and hypo data 
###the data should at least have five columns
###the strand information could use star(*) to replace if there no strand information
# chr  start    end strand  diff
#   1 904137 904137      -  -20.16218
#   1 904198 904198      -  -22.44681
#   1 914743 914743      +  -22.71028
dhype<-makeGRangesFromDataFrame(hyper,keep.extra.columns=T)
dhypo<-makeGRangesFromDataFrame(hypo,keep.extra.columns=T)
value(dhype)$id<-"hyper"
value(dhypo)$id<-"hypo"
###make the figure
p <- ggplot() + layout_circle(mychr, geom = "ideo", fill = "grey60",radius = 39, trackWidth = 2)
p <- p + layout_circle(c(dhype,dhypo), geom = "point",size = 1, aes(x = midpoint,y = diff, color = id), radius = 25, trackWidth = 30) +
      scale_colour_manual(values = c("magenta","green")) 
        p + layout_circle(mychr, geom = "text", aes(label = seqnames),vjust = 0, radius = 55, trackWidth = 7)
print(p)
ggsave(file="cir_difference.pdf")


