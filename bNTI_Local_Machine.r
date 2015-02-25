
## change to the directory on your computer that contains the OTU table and phylogeny
## note that the 'slash' needs to be changed to a forward slash like this /
setwd("C:/Users/steg815/Desktop/Stegen_PNNL/")

## load this library
## if not already installed, use install.packages('picante')
library(picante)

## read in OTU table

otu = read.csv("bacteria-abundance of OTU.csv",header=T,row.names=1);
dim(otu); # this gives the dimensions
otu[1:5,1:5]; # this gives a look at the first 5 rows and columns

## read in the phylogeny

phylo = read.tree("bacteria-phylogeny.txt");
phylo; # a summary of the phylogeny
plot.phylo(phylo,typ="fan"); # a quick plot

## make sure the names on the phylogeny are ordered the same as the names in otu table

match.phylo.otu = match.phylo.data(phylo, t(otu));
str(match.phylo.otu);

## calculate empirical betaMNTD

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
write.csv(beta.mntd.weighted,'betaMNTD_weighted.csv',quote=F);

identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE

# calculate randomized betaMNTD

beta.reps = 999; # number of randomizations

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
write.csv(weighted.bNTI,"weighted_bNTI.csv",quote=F);

pdf("weighted_bNTI_Histogram.pdf")
  hist(weighted.bNTI)
dev.off()

