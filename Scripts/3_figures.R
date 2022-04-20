########################
#plot metacoder figures#
########################

require(dplyr)
require(reshape2)
require(metacoder)
require(ggplot2)
require(ggpubr)
require(lubridate)
require(plyr)

#create functions for metacoder
plots <- function(season, meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.001, 0.01),
                  node_color = temp$data$carnivore_occ[[season]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.015),
                  tree_color_range = tree_color_range,
                  layout="davidson-harel",
                  node_color_axis_label="season focc",
                  node_size_axis_label="total focc")
  return(ht)
}

key_plot <- function(meta.obj, phylum, tree_color_range, seednum, nobs = 0)
{
  set.seed(seednum)
  temp <- meta.obj %>%
    filter_taxa(n_obs > nobs, reassign_obs = F)
  ht <- heat_tree(temp,
                  node_label = taxon_names, 
                  node_label_size_range = c(0.02, 0.04),
                  #node_size=n_obs,
                  node_size = temp$data$carnivore_occ[["total"]],
                  node_size_range = c(0.001, 0.015),
                  node_color = temp$data$carnivore_occ[["total"]],
                  node_color_range = tree_color_range,
                  edge_size_range = c(0.005, 0.01),
                  tree_color_range = tree_color_range,
                  layout="davidson-harel",
                  node_color_axis_label="freq occ")
  return(ht)
}

############
#load data

#all data
data.all <- read.table("Data_Raw/dataall_meta_man_2022-04-18.txt", sep=",") #data with metabarcoding and manual sorting
sk.final <- read.table(file="Data_Raw/spilogale_gracilis_metabarcoding_manual_byphylum_2022-01-25.txt", sep=",", header=T) #NA values in rows for MT904 and F37-193, samples only had clupea
sk.final <- sk.final[!sk.final$ScatID %in% c("MT904","F37-193"),] #remove MT904 and F37-193 from sk.final, samples only had clupea

#give data.all new season information
data.all <- merge(data.all, sk.final[,c("ScatID","CollectionType","CollectionDate","scatseason","logged","e","e.m","YrsSinceDi")], by.x="sample", by.y="ScatID", all.x=T)

data.all$lineage <- paste(data.all$p, data.all$c, data.all$o, data.all$f, data.all$g, data.all$s, sep=",") #make lineages line up across metabarcoding data and manual sort data
data.all$lineage <-  gsub(data.all$lineage, pattern=",{2,}", replacement="")

#vertebrate data
vertdata <- data.all[data.all$carnivore == "Spilogale_gracilis" & data.all$p == "Chordata",]
#remove MT906 and MT907 because probably clupea bones via manual sort
vertdata <- vertdata[!vertdata$sample %in% c("MT906","MT907"),]

#summarize number of verts per scat
vertcount <- ddply(vertdata, .(sample), nrow)
#
skunksamples <- unique(vertdata[,c("sample","scatseason")])
table(skunksamples$scatseason)
#
# F  S
#45 30
#77 samples total

sk <- vertdata %>% select(sample, lineage, numreads)
sk <- melt(sk, id=c("lineage","sample"))
sk <- dcast(data=sk, lineage~sample, fun.aggregate=length)

obj <- parse_tax_data(sk, class_cols="lineage", class_sep=",")
obj$data$carnivore_abund <- calc_taxon_abund(obj, 'tax_data', cols=unique(vertdata$sample))
obj$data$carnivore_occ <- calc_n_samples(obj, 'carnivore_abund', cols=unique(vertdata$sample), groups=unique(vertdata[,c("sample","scatseason")])$scatseason)
obj$data$carnivore_occ$total <- obj$data$carnivore_occ$wet + obj$data$carnivore_occ$dry

v_color <- c("azure3" , "sienna3" , "darkcyan")

f <- plots("wet", meta.obj=obj, phylum="spilogale", tree_color_range=v_color, seednum=2)
s <- plots("dry", meta.obj=obj, phylum="spilogale", tree_color_range=v_color, seednum=2)
vkey <- key_plot(meta.obj=obj, phylum="spilogale", tree_color_range=v_color, seednum=2)

# ggsave(vkey, filename=paste("Figures/spilogale_all_bestmatch_vertebrates_", Sys.Date(), ".tiff", sep=""), height=8, width=8, units="in", compression="lzw", dpi=300)
# ggsave(ggarrange(s,f, nrow=1), filename=paste("Figures/spilogale_all_bestmatch_vertebrates_byseason_", Sys.Date(), ".tiff", sep=""), height=8, width=16, units="in", compression="lzw", dpi=300)

######
#plant data
plantdata <- data.all[data.all$carnivore == "Spilogale_gracilis" & data.all$p == "Streptophyta",]
plantdata$lineage <- paste(plantdata$p, plantdata$c, plantdata$o, plantdata$f, plantdata$g, sep=",") #only go up to genus, since so many congeners
plantdata$lineage <-  gsub(plantdata$lineage, pattern=",{2,}", replacement="")

plantsamples <- unique(plantdata[,c("sample","scatseason")])
table(plantsamples$scatseason)

# F  S 
#19 18 
#37 samples

pl <- plantdata %>% select(sample, lineage, numreads)
pl <- melt(pl, id=c("lineage","sample"))
pl <- dcast(data=pl, lineage~sample, fun.aggregate=length)
pl2 <- pl
pl2$focc <- rowSums(pl2[,-1])
pl2 <- pl2[rev(order(pl2$focc)),]

#being more conservative for plant species since could be contamination
pl <- pl[rowSums(pl[,-1]) > 1,] #remove species only present in one sample

obj.pl <- parse_tax_data(pl, class_cols="lineage", class_sep=",")
obj.pl$data$carnivore_abund <- calc_taxon_abund(obj.pl, 'tax_data', cols=unique(plantdata$sample))
obj.pl$data$carnivore_occ <- calc_n_samples(obj.pl, 'carnivore_abund', cols=unique(plantdata$sample), groups=unique(plantdata[,c("sample","scatseason")])$scatseason)

obj.pl$data$carnivore_occ$total <- obj.pl$data$carnivore_occ$wet + obj.pl$data$carnivore_occ$dry

nobs <- 0
seednum <- 4 #for plants
f <- plots("wet", meta.obj=obj.pl, phylum="spilogale", tree_color_range=v_color, seednum=seednum, nobs = nobs)
s <- plots("dry", meta.obj=obj.pl, phylum="spilogale", tree_color_range=v_color, seednum=seednum, nobs = nobs)
pkey <- key_plot(meta.obj=obj.pl, phylum="spilogale", tree_color_range=v_color, seednum=seednum, nobs = nobs)

# ggsave(pkey, filename=paste("Figures/spilogale_all_bestmatch_plants_", Sys.Date(), ".tiff", sep=""), height=8, width=8, units="in", compression="lzw", dpi=300)
# ggsave(ggarrange(pkey, s, f, nrow=1, labels="AUTO", font.label=list(size=24)),
#        filename=paste("Figures/spilogale_season_bestmatch_plants_", Sys.Date(), ".tiff", sep=""), height=8, width=20, units="in", compression="lzw", dpi=300)

##########
#invert data
bugdata <- data.all[data.all$carnivore == "Spilogale_gracilis" & (data.all$p == "Arthropoda" | data.all$p == "Mollusca"),]
bugdata$lineage <- paste(bugdata$superk, bugdata$p, bugdata$c, bugdata$o, bugdata$f, bugdata$g, sep=",") #include superk so gastropods and arthropods on same figure
bugdata$lineage <-  gsub(bugdata$lineage, pattern=",{2,}", replacement="")

bugsamples <- unique(bugdata[,c("sample","scatseason")])
table(bugsamples$scatseason)

#wet dry 
# 55  54
#109 samples

bugs <- bugdata %>% select(sample, lineage, numreads)
bugs <- melt(bugs, id=c("lineage","sample"))
bugs <- dcast(data=bugs, lineage~sample, fun.aggregate=length)

obj.bugs <- parse_tax_data(bugs, class_cols="lineage", class_sep=",")
obj.bugs$data$carnivore_abund <- calc_taxon_abund(obj.bugs, 'tax_data', cols=unique(bugdata$sample))
obj.bugs$data$carnivore_occ <- calc_n_samples(obj.bugs, 'carnivore_abund', cols=unique(bugdata$sample), groups=unique(bugdata[,c("sample","scatseason")])$scatseason)

obj.bugs$data$carnivore_occ$total <- obj.bugs$data$carnivore_occ$wet + obj.bugs$data$carnivore_occ$dry

seednum <- 2
f <- plots("wet", meta.obj=obj.bugs, phylum="spilogale", tree_color_range=v_color, seednum=seednum)
s <- plots("dry", meta.obj=obj.bugs, phylum="spilogale", tree_color_range=v_color, seednum=seednum)
ikey <- key_plot(meta.obj=obj.bugs, phylum="spilogale", tree_color_range=v_color, seednum=seednum)

# ggsave(ikey, filename=paste("Figures/spilogale_all_bestmatch_bugs_", Sys.Date(), ".tiff", sep=""), height=8, width=8, units="in", compression="lzw", dpi=300)

############
#figure with 3 panels: verts, inverts, plants
ggarrange(vkey, ikey, pkey, nrow=1, labels="AUTO", align="hv")

############
#code each scat by scat contents (phylum)

skunksamples #75 samples
plantsamples #37 samples
bugsamples #109 samples

#calculate counts by taxonomic class and season
# make plants all the same c to tabulate below
# data.all$c <- as.character(data.all$c)
# data.all[data.all$p == "Streptophyta",]$c <- "plant"

data.long <- ddply(data.all, .(sample, superk, k, p, c), nrow)
data.long[data.long$c == "",]$c <- data.long[data.long$c == "",]$p #some rows missing taxonomic class information
data.long <- data.long[!duplicated(data.long[,c("sample","c")]),]
data.wide <- tidyr::spread(data.long[,c("sample","c","V1")], c, V1, fill=0)
data.wide$plant <- data.wide$Magnoliopsida + data.wide$Pinopsida + data.wide$Polypodiopsida + data.wide$Streptophyta
data.wide$myriapoda <- data.wide$Chilopoda + data.wide$Diplopoda
data.wide <- data.wide[,c("sample","Amphibia","Aves","Mammalia","Reptilia","Gastropoda","Arachnida","myriapoda","Insecta","plant")]

classes <- data.wide[,-1]
classes[classes > 0] <- 1
colSums(classes)

# Amphibia         Aves     Mammalia     Reptilia   Gastropoda    Arachnida    myriapoda      Insecta        plant 
#       17           18           60            3            9            3           45           92           37

table(sk.final$scatseason)
#57 dry, 71 wet

sk.vert <- data.frame(sample=skunksamples$sample, season=skunksamples$scatseason)
sk.vert$vert <- 1
nrow(sk.vert) #75 skunk scats with vertebrates

sk.plant <- data.frame(sample=plantsamples$sample, season=plantsamples$scatseason)
sk.plant$plant <- 1
nrow(sk.plant) #37 skunk scats with plants

sk.bug <- data.frame(sample=bugsamples$sample, season=bugsamples$scatseason)
sk.bug$bug <- 1
nrow(sk.bug) #109 skunk scats with plants

###########
#bar charts
#tabulate number of occurences and mean rra
###########

#calculate overall frequency of occurrence and relative read abundance
require(plyr)
vert.colors <- c(#"#6DC0D5",
  "#5A716A","#EAC435","#64403E","#E0D3DE","#FF6201")
legend.params <- theme(legend.position = c(0.95, 0.025),
                       legend.justification = c("right", "bottom"),
                       legend.box.just = "right",
                       legend.margin = margin(6, 6, 6, 6))

s.tab <- ddply(data.all[data.all$p == "Chordata",], .(superk, k, p, c, o, f, g, s), summarize, freq=length(rra), focc=length(rra)/nrow(sk.vert), finalrra_all=sum(rra)/nrow(sk.final), finalrra_vert=sum(rra)/nrow(sk.vert))
s.tab$c <- factor(s.tab$c)
s.tab$sciname <- as.character(s.tab$s)
s.tab[s.tab$sciname == "",]$sciname <- paste(s.tab[s.tab$sciname == "",]$g, "sp", sep=" ")
s.tab[s.tab$sciname == " sp",]$sciname <- paste(s.tab[s.tab$sciname == " sp",]$o, s.tab[s.tab$sciname == " sp",]$sciname, sep="")
s.tab[s.tab$sciname == " sp",]$sciname <- paste(s.tab[s.tab$sciname == " sp",]$c, s.tab[s.tab$sciname == " sp",]$sciname, sep="")
s.tab$sciname <- factor(s.tab$sciname, levels=s.tab[order(s.tab$finalrra_vert, decreasing = F),]$sciname) #order by vert RRA
#s.tab$sciname <- factor(s.tab$sciname, levels=s.tab[order(s.tab$c, decreasing = F),]$sciname) #order by class
#s.tab$sciname <- factor(s.tab$sciname, levels=s.tab[order(s.tab$c, s.tab$finalrra_vert),]$sciname) #order by class, then by vert RRA
#s.tab$sciname <- factor(s.tab$sciname, levels=s.tab[order(s.tab$finalrra_all, decreasing = F),]$sciname) #order by total RRA

rra <- ggplot(s.tab[!is.na(s.tab$finalrra_all),], aes(x=sciname, y=finalrra_vert, fill=c)) + 
  #rra <- ggplot(s.tab, aes(x=sciname, y=finalrra_all, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="vertebrate class", values=vert.colors) + 
  legend.params +
  ylab("relative read abundance") +
  xlab("species name")
rra
# ggsave(rra, filename="Figures/Spilogale_gracilis_rra_vertebrate.tiff", height=8, width=6, units="in", dpi=300, compression="lzw")

s.tab$sciname <- factor(s.tab$sciname, levels=s.tab[order(s.tab$focc, decreasing = F),]$sciname)
#s.tab$sciname <- factor(s.tab$sciname, levels=s.tab[order(s.tab$c, s.tab$focc),]$sciname) #order by class, then by focc
focc <- ggplot(s.tab[!s.tab$o == "",], aes(x=sciname, y=focc, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="prey class", values=vert.colors) + 
  legend.params +
  ylab("conditional FOO") +
  xlab("species name") +
  theme(axis.text.y = element_text(size=10))
focc
# ggsave(focc, filename="Figures/Spilogale_gracilis_focc_vertebrate.tiff", height=4, width=6.5, units="in", dpi=300, compression="lzw")

fr <- ggplot(s.tab[!s.tab$o == "",], aes(x=sciname, y=freq, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="prey class", values=vert.colors) + 
  legend.params +
  ylab("number of occurrences") +
  xlab("species name")
fr
# ggsave(fr, filename="Figures/Spilogale_gracilis_freq_vertebrate.tiff", height=8, width=6, units="in", dpi=300, compression="lzw")

###
#plant.colors <- c("#0B3948","#7A8E90","#D0F0C0")
plant.colors <- c("#5B9279","#444054","#8FCB9B")
p.tab <- ddply(data.all[data.all$p == "Streptophyta",], .(superk, k, p, subp, c, o, f, g), summarize, freq=length(rra), focc=length(rra)/nrow(sk.plant), finalrra_all=sum(rra)/nrow(sk.final), finalrra_plant=sum(rra)/nrow(sk.plant))
p.tab <- p.tab[!p.tab$f == "",]
p.tab$gname <- factor(p.tab$g, levels=p.tab[order(p.tab$finalrra_plant, decreasing=F),]$g)
rra2 <- ggplot(p.tab[p.tab$freq > 1,], aes(x=gname, y=finalrra_plant, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="vegetation class", values=plant.colors) + 
  legend.params +
  ylab("relative read abundance") +
  xlab("genus name")
rra2
# ggsave(rra2, filename="Figures/Spilogale_gracilis_rra_plant.tiff", height=8, width=6, units="in", dpi=300, compression="lzw") 

p.tab$g <- factor(p.tab$g, levels=p.tab[order(p.tab$focc, decreasing = F),]$g)
focc2 <- ggplot(p.tab[p.tab$freq > 1,], aes(x=g, y=focc, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="prey class", values=plant.colors) + 
  legend.params +
  ylab("conditional FOO") +
  xlab("genus name")
focc2
# ggsave(focc2, filename="Figures/Spilogale_gracilis_focc_plant.tiff", height=4, width=5.5, units="in", dpi=300, compression="lzw")

####
bugcolors <- c("#943F0A","#51344D","#FF934F","#AB8D7B")
bugdata <- ddply(bugdata, .(sample, locus, superk, k, p, c, o, f, g, s, carnivore, lineage), summarize, rra=mean(rra), numreads=sum(numreads), totalreads=mean(totalreads))
b.tab <- ddply(bugdata, .(superk, k, p, c, o, f, g), summarize, freq=length(rra), focc=length(rra)/nrow(sk.bug), finalrra_all=sum(rra)/130, finalrra_bug=sum(rra)/nrow(sk.bug))
b.tab$gname <- as.character(b.tab$g)
b.tab[b.tab$gname == "",]$gname <- as.character(b.tab[b.tab$g == "",]$f) #fill in gname
b.tab[b.tab$f == "",]$gname <- as.character(b.tab[b.tab$f == "",]$o) #fill in gname
b.tab[b.tab$o == "",]$gname <- as.character(b.tab[b.tab$o == "",]$c) #fill in gname

b.tab$gname <- factor(b.tab$gname, levels=b.tab[order(b.tab$finalrra_bug, decreasing=F),]$gname)
rra3 <- ggplot(b.tab[!is.na(b.tab$finalrra_all),], aes(x=gname, y=finalrra_bug, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="invertebrate class", values=bugcolors) + 
  legend.params +
  ylab("relative read abundance") +
  xlab("genus name")
rra3
# ggsave(rra3, filename="Figures/Spilogale_gracilis_rra_bug.tiff", height=8, width=6, units="in", dpi=300, compression="lzw") 

bugcolors <- c("#943F0A","#51344D","#FCBB92","#FF934F","#AB8D7B")
b.tab$gname <- factor(b.tab$gname, levels=b.tab[order(b.tab$focc, decreasing = F),]$gname)
focc3 <- ggplot(b.tab, aes(x=gname, y=focc, fill=c)) + 
  geom_bar(stat="identity") + theme_bw(base_size=15) + coord_flip() +
  scale_fill_manual(name="prey class", values=bugcolors) + 
  legend.params +
  ylab("conditional FOO") +
  xlab("genus name")
focc3
# ggsave(focc3, filename="Figures/Spilogale_gracilis_focc_bug.tiff", height=4, width=4.5, units="in", dpi=300, compression="lzw")
