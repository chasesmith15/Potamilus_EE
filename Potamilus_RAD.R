<<<<<<< HEAD
#######R script for Potamilus RAD filtering and analysis#########
library(dartR)
library(adegenet)
library(magrittr)

####Bring in DArT Data
#2 SNP data
PotDArT <- gl.read.dart(filename="Report_DPota20-4942_SNP_2.csv", covfilename = "Pot_cov.csv") %>%
  #Remove failed sample
  gl.recode.ind(ind.recode="rm_PohiBra036.csv")  

#Clean all
glpot_clean <- 
  #Quality check
  gl.filter.reproducibility(PotDArT, t=1) %>% 
  #1 SNP per locus and filter loci that are called at percentage e.g., 0.9 = 90%
  gl.filter.callrate(method = "loc", threshold = 0.70) %>%   
  #Remove individuals with >10% missing data
  gl.filter.callrate(method = "ind", threshold = 0.70) %>%
  #Remove maf
  gl.filter.maf(threshold = 0.05) %>%
  #1 SNP per locus
  gl.filter.secondaries(method = "best") %>%

#Linkage
library(SNPRelate)
gl2gds(glpot_clean, outfile="Pot_all.gds", outpath = "C:/Users/Chase Smith/Desktop/Chase Files/Completed Papers/Pamp_str_RAD/Report-DPota20-4942")
gds <- snpgdsOpen("Pot_all.gds")
snpgdsLDpruning(gds, maf=0.05, missing.rate = 0.3)
glpot_clean1 <- gl.drop.loc(glpot_clean, loc.list = c("55366436-39-A/T","55366461-55-T/C","55366479-61-A/T","55366505-66-C/T","55366506-65-G/T",
                            "55366509-61-G/T","55366555-27-G/A","55366742-24-C/A","55366882-22-A/C","55367206-52-A/G",
                            "55367271-67-C/T","55367279-35-G/A","55367384-48-T/G","55367448-40-T/A","55367803-28-A/C",
                            "55368009-11-A/G","55368132-58-T/C","55368420-49-T/C","55370687-26-A/T","55370880-35-G/A",
                            "55371153-45-A/T","55371219-28-G/T","55371377-31-T/A","55371472-17-C/T","55371496-29-T/A",
                            "55371880-35-T/A","55371982-47-G/A","55372246-20-G/A","55372538-22-C/T","55372868-48-T/C",
                            "55372928-33-G/T","55373799-29-T/G","55374237-15-T/G","55374652-47-T/A","55375089-39-A/T",
                            "55375102-65-G/A","55377855-18-A/G","55379046-33-C/T","55379481-58-C/A","55380191-22-T/C",
                            "55380772-68-A/G","55381995-63-C/A","55382372-12-G/A","55383620-12-T/C","55383672-35-C/G",
                            "55384378-16-T/G","55388665-33-G/A","55390481-36-G/A","55391579-25-A/G","55391693-7-G/A",
                            "55393740-21-A/G","55394789-67-C/T","55396291-30-C/T","55399115-22-C/A","55401475-42-G/A",
                            "55401651-49-A/G","55403922-18-G/A", "55739788-10-G/A"))

#Create frame only for amphichaenus
glamp_clean <- gl.recode.ind(glpot_clean1, ind.recode="amponly_ind_assignments.csv")
gl2gds(glamp_clean, outfile="Pot_amp.gds", outpath = "C:/Users/Chase Smith/Desktop/Chase Files/Completed Papers/Pamp_str_RAD/Report-DPota20-4942")
gds1 <- snpgdsOpen("Pot_amp.gds")
snpgdsLDpruning(gds1, maf=0.05, missing.rate = 0.3)
glamp_clean1 <- gl.drop.loc(glamp_clean, loc.list = c("55371782-16-A/G","55379526-34-T/A","55384337-10-G/A","55384619-36-G/A","55384626-24-G/A",
                                                      "55384633-48-C/A","55384644-38-C/A","55384646-46-C/G","55384664-26-C/T","55384665-25-T/C",
                                                      "55384677-37-T/G","55384688-8-G/T","55384725-9-A/G","55384750-41-G/T","55384763-60-A/T",
                                                      "55384778-57-C/T","55384879-14-A/T","55384892-38-G/C","55384976-20-T/A","55384993-9-A/T",
                                                      "55384994-39-C/A","55385097-55-T/A","55385106-54-T/A","55385110-43-T/A","55385263-35-C/T",
                                                      "55385559-18-C/T","55386219-30-C/A","55386287-42-C/T","55386300-34-G/A","55386495-25-A/G",
                                                      "55386695-66-C/T","55386948-68-C/A","55387083-52-T/A","55387633-10-G/A","55388396-67-T/A",
                                                      "55388827-15-T/C","55390330-48-T/G","55401319-17-G/A","55401834-17-C/A")) %>%
  gl.filter.monomorphs()

#Create frame only for streckersoni
glstr_clean <- gl.recode.ind(glpot_clean1, ind.recode="stronly_ind_assignments.csv")
gl2gds(glstr_clean, outfile="Pot_str.gds", outpath = "C:/Users/Chase Smith/Desktop/Chase Files/Completed Papers/Pamp_str_RAD/Report-DPota20-4942")
gds2 <- snpgdsOpen("Pot_str.gds")
snpgdsLDpruning(gds2, maf=0.05, missing.rate = 0.3)
glstr_clean1 <- gl.drop.loc(glstr_clean, loc.list = c("55371372-58-C/T","55371384-44-A/C","55371402-52-C/T","55371422-30-A/T","55371441-22-T/C",
                                                      "55371515-39-A/T","55371588-20-A/T","55371718-24-G/A","55371743-57-A/G","55371747-62-C/G",
                                                      "55371753-11-C/G","55372561-37-A/G","55373070-34-T/G","55373257-64-C/T","55373965-64-T/A",
                                                      "55376199-21-T/C","55376970-40-A/G","55384477-47-G/A","55386392-67-T/C","55399305-18-T/C",
                                                      "55399928-16-C/T"))%>%
  gl.filter.monomorphs()
#Convert to genind and genepop
potgi <- gl2gi(glpot_clean1)

#Basic stats
gl.report.pa(glpot_clean1)
gl.fixed.diff(glpot_clean1)

library(diveRsity)
gl2structure(glpot_clean1, outfile="Pot_all.str", outpath = "C:/Users/Chase Smith/Desktop")
#Converted to genepop
basicStats(infile = 'Pot_Gp.txt', outfile = 'Diversity', fis_ci = T,
           ar_ci = T, fis_boots = 999, ar_boots = 999, 
           mc_reps = 9999, rarefaction = TRUE, ar_alpha = 0.05, 
           fis_alpha = 0.05)

#####FSTATS
library(StAMPP)
#Fst
potFst <-stamppFst(glpot_clean1, nboots=999, percent=95, nclusters=8)
potFst

#####Visualization
#PCoA
pc <- gl.pcoa(glpot_clean1, n.cores = 8)
gl.pcoa.scree(pc)
plot(pc$eig)
library(ggplot2)
gl.pcoa.plot(pc, glpot_clean1, labels = "pop", xaxis=1, yaxis=2, ellipse = T)+
  labs(color="Taxon")+
  theme_classic()
library(reshape2)
library(plotly)
ggplotly()
ggsave("Fig.1_PCoA.pdf", dpi = 600, device = "pdf")


#####DAPC
###All Pot
#K Means cluster and DAPC
grp <- find.clusters.genind(potgi, n.pca = 2, method = "kmeans")

#Groups
grp$size
table(pop(potgi), grp$grp)
table.value(table(pop(potgi), grp$grp), col.lab=paste("inf", 1:3), row.lab=paste("ori", 1:4))

#dapc
dapc1 <- dapc(potgi, grp$grp, n.pca = 2, n.da = 2)
scatter(dapc1, scree.da = F)
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, solid=.4, cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:4))
summary(dapc1)

#####FastStructure###########
gl2faststructure(glpot_clean1, outfile="Pot_all_FS.str", outpath = "C:/Users/Chase Smith/Desktop/")
gl2faststructure(glamp_clean1, outfile="Pot_amp_FS.str", outpath = "C:/Users/Chase Smith/Desktop/")
gl2faststructure(glstr_clean1, outfile="Pot_str_FS.str", outpath = "C:/Users/Chase Smith/Desktop/")
gl2fasta(glpot_clean1, outfile="Pot_all.fas", outpath = "C:/Users/Chase Smith/Desktop")

#####TESS3##########
#All taxa
geno2 <- glpot_clean1
geno <- as.matrix(geno2)
sample <- row.names(geno)
pop.names <- pop(geno2)
ploidy <- ploidy(geno2)
geno = geno
geno[is.na(geno)] = NaN
format <- vector(length = length(geno[, 1]))
format[1:length(geno[, 1])] = "genlight"
pops <- unique(pop.names)
pop.num <- vector(length = length(geno[, 1]))
for (i in 1:length(geno[, 1])) {
  pop.num[i] = which(pop.names[i] == pops)
}
genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                               ploidy, format))
geno <- cbind(genoLHS, geno)
geno[, 2] = as.character(pop.names)
geno[, 4] = as.numeric(as.character(geno[, 4]))
row.names(geno) = NULL

coordinates <- potgi$other$ind.metrics[,c(3,4)]
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

library(tess3r)
# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:10, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 2
q.matrix <- qmatrix(tess3.obj, K = 2)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)


###Map the results for K = 4####
library(ggplot2)
library(rworldmap)
library(maps)
library(tidyverse)
library(maptools)
library(rgdal)
library(MazamaSpatialUtils)
#DEFINE URL
# - this is the location of the file
url.river_data <- url("http://sharpsightlabs.com/wp-content/datasets/usa_rivers.RData")
# LOAD DATA
load(url.river_data)
# INSPECT
summary(lines.rivers)
lines.rivers@data %>% glimpse()
levels(lines.rivers$FEATURE)
table(lines.rivers$FEATURE)

# REMOVE MISC FEATURES
lines.rivers <- subset(lines.rivers, !(FEATURE %in% c("Shoreline"
                                                      ,"Shoreline Intermittent"
                                                      ,"Null"
                                                      ,"Closure Line"
                                                      ,"Apparent Limit")))


lines.rivers <- subset(lines.rivers, (STATE %in% c('TX')))
df.usa_rivers <- fortify(lines.rivers)

map.polygon <- getMap(resolution = "low")

#HUCS
ogrListLayers("NHD_H_Texas_State_GDB.gdb")
WBD8 <- readOGR("NHD_H_Texas_State_GDB.gdb", layer = "WBDHU8", stringsAsFactors = F)

#Lakes
library(sf)
library(ggplot2)
lakes <- st_read("C:/Users/Chase Smith/Desktop/Pamp_str_RAD/Distribution/ne_10m_lakes.shp")
lakes_add <- st_read("C:/Users/Chase Smith/Desktop/Pamp_str_RAD/Distribution/ne_10m_lakes_north_america.shp")
abundance <- read.csv(file = "C:/Users/Chase Smith/Desktop/Pamp_str_RAD/Distribution/Test_abundance.csv", header = T)
pl <- ggtess3Q(q.matrix, coordinates, map.polygon = map.polygon)
pl +
  geom_path(data = df.usa_rivers, aes(x = long, y = lat, group = group), color = "blue") +
  geom_polygon(data = WBD8, aes(x = long, y = lat), color="black", fill=NA) +
  geom_sf(data = lakes, color = "blue", fill = "NA") +
  geom_sf(data = lakes_add, color = "blue", fill = "NA") +
  geom_point(data = abundance, aes(x = Longitude, y = Latitude, size = CPUE))+
  xlim(-97.5, -93.7) + 
  ylim(29.4, 32.8) + 
  coord_sf() + 
  #geom_point(data = as.data.frame(coordinates), aes(x = long, y = lat), size = 3, shape = 21, color = "Pink") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#str only
strgi <- gl2gi(glstr_clean1)
geno2 <- glstr_clean1
geno <- as.matrix(geno2)
sample <- row.names(geno)
pop.names <- pop(geno2)
ploidy <- ploidy(geno2)
geno = geno
geno[is.na(geno)] = NaN
format <- vector(length = length(geno[, 1]))
format[1:length(geno[, 1])] = "genlight"
pops <- unique(pop.names)
pop.num <- vector(length = length(geno[, 1]))
for (i in 1:length(geno[, 1])) {
  pop.num[i] = which(pop.names[i] == pops)
}
genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                               ploidy, format))
geno <- cbind(genoLHS, geno)
geno[, 2] = as.character(pop.names)
geno[, 4] = as.numeric(as.character(geno[, 4]))
row.names(geno) = NULL

coordinates <- strgi$other$ind.metrics[,c(3,4)]
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:10, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 1
q.matrix <- qmatrix(tess3.obj, K = 1)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

#amp only
ampgi <- gl2gi(glamp_clean1)
geno2 <- glamp_clean1
geno <- as.matrix(geno2)
sample <- row.names(geno)
pop.names <- pop(geno2)
ploidy <- ploidy(geno2)
geno = geno
geno[is.na(geno)] = NaN
format <- vector(length = length(geno[, 1]))
format[1:length(geno[, 1])] = "genlight"
pops <- unique(pop.names)
pop.num <- vector(length = length(geno[, 1]))
for (i in 1:length(geno[, 1])) {
  pop.num[i] = which(pop.names[i] == pops)
}
genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                               ploidy, format))
geno <- cbind(genoLHS, geno)
geno[, 2] = as.character(pop.names)
geno[, 4] = as.numeric(as.character(geno[, 4]))
row.names(geno) = NULL

coordinates <- ampgi$other$ind.metrics[,c(3,4)]
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:10, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 2
q.matrix <- qmatrix(tess3.obj, K = 2)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)
=======
#######R script for Potamilus RAD filtering and analysis#########
library(dartR)
library(adegenet)
library(magrittr)

####Bring in DArT Data
#2 SNP data
PotDArT <- gl.read.dart(filename="Report_DPota20-4942_SNP_2.csv", covfilename = "Pot_cov.csv") %>%
  #Remove failed sample
  gl.recode.ind(ind.recode="rm_PohiBra036.csv")  

#Clean all
glpot_clean <- 
  #Quality check
  gl.filter.reproducibility(PotDArT, t=1) %>% 
  #1 SNP per locus and filter loci that are called at percentage e.g., 0.9 = 90%
  gl.filter.callrate(method = "loc", threshold = 0.70) %>%   
  #Remove individuals with >10% missing data
  gl.filter.callrate(method = "ind", threshold = 0.70) %>%
  #Remove maf
  gl.filter.maf(threshold = 0.05) %>%
  #1 SNP per locus
  gl.filter.secondaries(method = "best") %>%

#Linkage
library(SNPRelate)
gl2gds(glpot_clean, outfile="Pot_all.gds", outpath = "C:/Users/Chase Smith/Desktop/Chase Files/Completed Papers/Pamp_str_RAD/Report-DPota20-4942")
gds <- snpgdsOpen("Pot_all.gds")
snpgdsLDpruning(gds, maf=0.05, missing.rate = 0.3)
glpot_clean1 <- gl.drop.loc(glpot_clean, loc.list = c("55366436-39-A/T","55366461-55-T/C","55366479-61-A/T","55366505-66-C/T","55366506-65-G/T",
                            "55366509-61-G/T","55366555-27-G/A","55366742-24-C/A","55366882-22-A/C","55367206-52-A/G",
                            "55367271-67-C/T","55367279-35-G/A","55367384-48-T/G","55367448-40-T/A","55367803-28-A/C",
                            "55368009-11-A/G","55368132-58-T/C","55368420-49-T/C","55370687-26-A/T","55370880-35-G/A",
                            "55371153-45-A/T","55371219-28-G/T","55371377-31-T/A","55371472-17-C/T","55371496-29-T/A",
                            "55371880-35-T/A","55371982-47-G/A","55372246-20-G/A","55372538-22-C/T","55372868-48-T/C",
                            "55372928-33-G/T","55373799-29-T/G","55374237-15-T/G","55374652-47-T/A","55375089-39-A/T",
                            "55375102-65-G/A","55377855-18-A/G","55379046-33-C/T","55379481-58-C/A","55380191-22-T/C",
                            "55380772-68-A/G","55381995-63-C/A","55382372-12-G/A","55383620-12-T/C","55383672-35-C/G",
                            "55384378-16-T/G","55388665-33-G/A","55390481-36-G/A","55391579-25-A/G","55391693-7-G/A",
                            "55393740-21-A/G","55394789-67-C/T","55396291-30-C/T","55399115-22-C/A","55401475-42-G/A",
                            "55401651-49-A/G","55403922-18-G/A", "55739788-10-G/A"))

#Create frame only for amphichaenus
glamp_clean <- gl.recode.ind(glpot_clean1, ind.recode="amponly_ind_assignments.csv")
gl2gds(glamp_clean, outfile="Pot_amp.gds", outpath = "C:/Users/Chase Smith/Desktop/Chase Files/Completed Papers/Pamp_str_RAD/Report-DPota20-4942")
gds1 <- snpgdsOpen("Pot_amp.gds")
snpgdsLDpruning(gds1, maf=0.05, missing.rate = 0.3)
glamp_clean1 <- gl.drop.loc(glamp_clean, loc.list = c("55371782-16-A/G","55379526-34-T/A","55384337-10-G/A","55384619-36-G/A","55384626-24-G/A",
                                                      "55384633-48-C/A","55384644-38-C/A","55384646-46-C/G","55384664-26-C/T","55384665-25-T/C",
                                                      "55384677-37-T/G","55384688-8-G/T","55384725-9-A/G","55384750-41-G/T","55384763-60-A/T",
                                                      "55384778-57-C/T","55384879-14-A/T","55384892-38-G/C","55384976-20-T/A","55384993-9-A/T",
                                                      "55384994-39-C/A","55385097-55-T/A","55385106-54-T/A","55385110-43-T/A","55385263-35-C/T",
                                                      "55385559-18-C/T","55386219-30-C/A","55386287-42-C/T","55386300-34-G/A","55386495-25-A/G",
                                                      "55386695-66-C/T","55386948-68-C/A","55387083-52-T/A","55387633-10-G/A","55388396-67-T/A",
                                                      "55388827-15-T/C","55390330-48-T/G","55401319-17-G/A","55401834-17-C/A")) %>%
  gl.filter.monomorphs()

#Create frame only for streckersoni
glstr_clean <- gl.recode.ind(glpot_clean1, ind.recode="stronly_ind_assignments.csv")
gl2gds(glstr_clean, outfile="Pot_str.gds", outpath = "C:/Users/Chase Smith/Desktop/Chase Files/Completed Papers/Pamp_str_RAD/Report-DPota20-4942")
gds2 <- snpgdsOpen("Pot_str.gds")
snpgdsLDpruning(gds2, maf=0.05, missing.rate = 0.3)
glstr_clean1 <- gl.drop.loc(glstr_clean, loc.list = c("55371372-58-C/T","55371384-44-A/C","55371402-52-C/T","55371422-30-A/T","55371441-22-T/C",
                                                      "55371515-39-A/T","55371588-20-A/T","55371718-24-G/A","55371743-57-A/G","55371747-62-C/G",
                                                      "55371753-11-C/G","55372561-37-A/G","55373070-34-T/G","55373257-64-C/T","55373965-64-T/A",
                                                      "55376199-21-T/C","55376970-40-A/G","55384477-47-G/A","55386392-67-T/C","55399305-18-T/C",
                                                      "55399928-16-C/T"))%>%
  gl.filter.monomorphs()
#Convert to genind and genepop
potgi <- gl2gi(glpot_clean1)

#Basic stats
gl.report.pa(glpot_clean1)
gl.fixed.diff(glpot_clean1)

library(diveRsity)
gl2structure(glpot_clean1, outfile="Pot_all.str", outpath = "C:/Users/Chase Smith/Desktop")
#Converted to genepop
basicStats(infile = 'Pot_Gp.txt', outfile = 'Diversity', fis_ci = T,
           ar_ci = T, fis_boots = 999, ar_boots = 999, 
           mc_reps = 9999, rarefaction = TRUE, ar_alpha = 0.05, 
           fis_alpha = 0.05)

#####FSTATS
library(StAMPP)
#Fst
potFst <-stamppFst(glpot_clean1, nboots=999, percent=95, nclusters=8)
potFst

#####Visualization
#PCoA
pc <- gl.pcoa(glpot_clean1, n.cores = 8)
gl.pcoa.scree(pc)
plot(pc$eig)
library(ggplot2)
gl.pcoa.plot(pc, glpot_clean1, labels = "pop", xaxis=1, yaxis=2, ellipse = T)+
  labs(color="Taxon")+
  theme_classic()
library(reshape2)
library(plotly)
ggplotly()
ggsave("Fig.1_PCoA.pdf", dpi = 600, device = "pdf")


#####DAPC
###All Pot
#K Means cluster and DAPC
grp <- find.clusters.genind(potgi, n.pca = 2, method = "kmeans")

#Groups
grp$size
table(pop(potgi), grp$grp)
table.value(table(pop(potgi), grp$grp), col.lab=paste("inf", 1:3), row.lab=paste("ori", 1:4))

#dapc
dapc1 <- dapc(potgi, grp$grp, n.pca = 2, n.da = 2)
scatter(dapc1, scree.da = F)
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, solid=.4, cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:4))
summary(dapc1)

#####FastStructure###########
gl2faststructure(glpot_clean1, outfile="Pot_all_FS.str", outpath = "C:/Users/Chase Smith/Desktop/")
gl2faststructure(glamp_clean1, outfile="Pot_amp_FS.str", outpath = "C:/Users/Chase Smith/Desktop/")
gl2faststructure(glstr_clean1, outfile="Pot_str_FS.str", outpath = "C:/Users/Chase Smith/Desktop/")
gl2fasta(glpot_clean1, outfile="Pot_all.fas", outpath = "C:/Users/Chase Smith/Desktop")

#####TESS3##########
#All taxa
geno2 <- glpot_clean1
geno <- as.matrix(geno2)
sample <- row.names(geno)
pop.names <- pop(geno2)
ploidy <- ploidy(geno2)
geno = geno
geno[is.na(geno)] = NaN
format <- vector(length = length(geno[, 1]))
format[1:length(geno[, 1])] = "genlight"
pops <- unique(pop.names)
pop.num <- vector(length = length(geno[, 1]))
for (i in 1:length(geno[, 1])) {
  pop.num[i] = which(pop.names[i] == pops)
}
genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                               ploidy, format))
geno <- cbind(genoLHS, geno)
geno[, 2] = as.character(pop.names)
geno[, 4] = as.numeric(as.character(geno[, 4]))
row.names(geno) = NULL

coordinates <- potgi$other$ind.metrics[,c(3,4)]
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

library(tess3r)
# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:10, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 2
q.matrix <- qmatrix(tess3.obj, K = 2)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)


###Map the results for K = 4####
library(ggplot2)
library(rworldmap)
library(maps)
library(tidyverse)
library(maptools)
library(rgdal)
library(MazamaSpatialUtils)
#DEFINE URL
# - this is the location of the file
url.river_data <- url("http://sharpsightlabs.com/wp-content/datasets/usa_rivers.RData")
# LOAD DATA
load(url.river_data)
# INSPECT
summary(lines.rivers)
lines.rivers@data %>% glimpse()
levels(lines.rivers$FEATURE)
table(lines.rivers$FEATURE)

# REMOVE MISC FEATURES
lines.rivers <- subset(lines.rivers, !(FEATURE %in% c("Shoreline"
                                                      ,"Shoreline Intermittent"
                                                      ,"Null"
                                                      ,"Closure Line"
                                                      ,"Apparent Limit")))


lines.rivers <- subset(lines.rivers, (STATE %in% c('TX')))
df.usa_rivers <- fortify(lines.rivers)

map.polygon <- getMap(resolution = "low")

#HUCS
ogrListLayers("NHD_H_Texas_State_GDB.gdb")
WBD8 <- readOGR("NHD_H_Texas_State_GDB.gdb", layer = "WBDHU8", stringsAsFactors = F)

#Lakes
library(sf)
library(ggplot2)
lakes <- st_read("C:/Users/Chase Smith/Desktop/Pamp_str_RAD/Distribution/ne_10m_lakes.shp")
lakes_add <- st_read("C:/Users/Chase Smith/Desktop/Pamp_str_RAD/Distribution/ne_10m_lakes_north_america.shp")
abundance <- read.csv(file = "C:/Users/Chase Smith/Desktop/Pamp_str_RAD/Distribution/Test_abundance.csv", header = T)
pl <- ggtess3Q(q.matrix, coordinates, map.polygon = map.polygon)
pl +
  geom_path(data = df.usa_rivers, aes(x = long, y = lat, group = group), color = "blue") +
  geom_polygon(data = WBD8, aes(x = long, y = lat), color="black", fill=NA) +
  geom_sf(data = lakes, color = "blue", fill = "NA") +
  geom_sf(data = lakes_add, color = "blue", fill = "NA") +
  geom_point(data = abundance, aes(x = Longitude, y = Latitude, size = CPUE))+
  xlim(-97.5, -93.7) + 
  ylim(29.4, 32.8) + 
  coord_sf() + 
  #geom_point(data = as.data.frame(coordinates), aes(x = long, y = lat), size = 3, shape = 21, color = "Pink") + 
  xlab("Longitude") +
  ylab("Latitude") + 
  theme_bw()

#str only
strgi <- gl2gi(glstr_clean1)
geno2 <- glstr_clean1
geno <- as.matrix(geno2)
sample <- row.names(geno)
pop.names <- pop(geno2)
ploidy <- ploidy(geno2)
geno = geno
geno[is.na(geno)] = NaN
format <- vector(length = length(geno[, 1]))
format[1:length(geno[, 1])] = "genlight"
pops <- unique(pop.names)
pop.num <- vector(length = length(geno[, 1]))
for (i in 1:length(geno[, 1])) {
  pop.num[i] = which(pop.names[i] == pops)
}
genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                               ploidy, format))
geno <- cbind(genoLHS, geno)
geno[, 2] = as.character(pop.names)
geno[, 4] = as.numeric(as.character(geno[, 4]))
row.names(geno) = NULL

coordinates <- strgi$other$ind.metrics[,c(3,4)]
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:10, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 1
q.matrix <- qmatrix(tess3.obj, K = 1)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

#amp only
ampgi <- gl2gi(glamp_clean1)
geno2 <- glamp_clean1
geno <- as.matrix(geno2)
sample <- row.names(geno)
pop.names <- pop(geno2)
ploidy <- ploidy(geno2)
geno = geno
geno[is.na(geno)] = NaN
format <- vector(length = length(geno[, 1]))
format[1:length(geno[, 1])] = "genlight"
pops <- unique(pop.names)
pop.num <- vector(length = length(geno[, 1]))
for (i in 1:length(geno[, 1])) {
  pop.num[i] = which(pop.names[i] == pops)
}
genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                               ploidy, format))
geno <- cbind(genoLHS, geno)
geno[, 2] = as.character(pop.names)
geno[, 4] = as.numeric(as.character(geno[, 4]))
row.names(geno) = NULL

coordinates <- ampgi$other$ind.metrics[,c(3,4)]
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:10, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 2
q.matrix <- qmatrix(tess3.obj, K = 2)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix") -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)
>>>>>>> 954d019b32cf57b727365f2adf3d339e2049e022
