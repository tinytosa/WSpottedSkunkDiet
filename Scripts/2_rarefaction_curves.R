#
# produce rarefraction curves for carnivore scat species diversity

# help with package from:
# https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html

#####
#load packages
require(iNEXT)
require(vegan)
require(ggplot2)
require(plyr)
require(tidyr) #for spread and gather

#####
#load data
data.all <- read.table("Data_Raw/dataall_meta_man_2022-04-18.txt", sep=",") #data with metabarcoding and manual sorting

sk.final <- read.table(file="Data_Raw/spilogale_gracilis_metabarcoding_manual_byphylum_2022-01-25.txt", sep=",", header=T) #NA values in rows for MT904 and F37-193
sk.final <- sk.final[!sk.final$ScatID %in% c("MT904","F37-193"),] #remove MT904 and F37-193 from sk.final, only clupea in samples

#separate data by taxon
#vertebrate data
vertdata <- data.all[data.all$carnivore == "Spilogale_gracilis" & data.all$p == "Chordata",]
vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),] #chordata here is Clupea
vertdata <- vertdata[!duplicated(vertdata[,c("lineage","sample")]),]
vert.wide <- spread(vertdata[,c("season","sample","lineage","numreads")], lineage, numreads, fill=0) #dims: 75 rows, 32 species, should be 128 rows
vert.wide <- merge(vert.wide, sk.final[,c("ScatID","scatseason")], by.x=c("sample"), by.y=c("ScatID"), all.y=T) #add in samples that don't have any verts
vert.wide[is.na(vert.wide)] <- 0 #replace NAs with 0
vert.wide[vert.wide > 0] <- 1

#plant data
plantdata <- data.all[data.all$carnivore == "Spilogale_gracilis" & data.all$p == "Streptophyta",]
plantdata$lineage <- paste(plantdata$p, plantdata$c, plantdata$o, plantdata$f, plantdata$g, sep=",") #only go up to genus, since so many congeners
plantdata <- plantdata[!duplicated(plantdata[,c("lineage","sample")]),]
plant.wide <- spread(plantdata[,c("season","sample","lineage","numreads")], lineage, numreads, fill=0) #dims: 37 rows, 44 species, should be 130 rows
plant.wide <- merge(plant.wide, sk.final[,c("ScatID","scatseason")], by.x=c("sample"), by.y=c("ScatID"), all.y=T) #add in samples that don't have any verts
plant.wide[is.na(plant.wide)] <- 0 #replace NAs with 0
plant.wide[plant.wide > 0] <- 1

#arthropod data
bugdata <- data.all[data.all$carnivore == "Spilogale_gracilis" & (data.all$p == "Arthropoda" | data.all$p == "Mollusca"),]
bugdata$lineage <- paste(bugdata$superk, bugdata$p, bugdata$c, bugdata$o, bugdata$f, bugdata$g, sep=",")
bugdata <- bugdata[!duplicated(bugdata[,c("lineage","sample")]),]
bug.wide <- spread(bugdata[,c("season","sample","lineage","numreads")], lineage, numreads, fill=0) #dims: 109 rows, 18 species, should be 128 rows
bug.wide <- merge(bug.wide, sk.final[,c("ScatID","scatseason")], by.x=c("sample"), by.y=c("ScatID"), all.y=T) #add in samples that don't have any verts
bug.wide[is.na(bug.wide)] <- 0 #replace NAs with 0
bug.wide[bug.wide > 0] <- 1

###############
#create and plot curves
data.inext <- list(
  vertebrate = t(vert.wide[,!names(vert.wide) %in% c("sample","season","scatseason"),]),
  invertebrate = t(bug.wide[,!names(bug.wide) %in% c("sample","season","scatseason"),]),
  plant = t(plant.wide[,!names(plant.wide) %in% c("sample","season","scatseason"),]))

out.inext <- iNEXT(x=data.inext, datatype="incidence_raw", endpoint=150) #data needs to be rows = species, columns = samples
p <- ggiNEXT(out.inext, type=1) + 
  xlab("Number of scat samples") + ylab("Taxonomic richness") +
  theme_bw(base_size=15) + theme(legend.position=c(0.15,0.75))
p
# ggsave(p, filename="Figures/Rarefaction_Spilogale_gracilis.tiff", height=6, width=8, units="in", dpi=300, compression="lzw")