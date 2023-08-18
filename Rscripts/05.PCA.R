
# Run PCA on gamma measures 

CreatePCAdimensions <- function(inputDF) {
sys.source("H:/Projects/7TBrainMech/scripts/eeg/Shane/Rscripts/01.PrepData.R", envir = knitr::knit_global(), chdir = TRUE)


gammavars <- grep("Gamma",names(inputDF),value=TRUE)
allGamma <- inputDF[,c("Subject", "age",gammavars)]

this <- allGamma[c(3:10)]
completeidx <- which(complete.cases(this))
this <- this[completeidx,]
s.pca <- prcomp(scale(this))

summary(s.pca)
print(s.pca$rotation)
print(ggplot(data=reshape2::melt(s.pca$rotation[,1:5])) + geom_bar(aes(x=Var2, y=value, fill=Var1), stat="identity", position='dodge'))

allGamma$eeg_pc1 <- NA
allGamma[completeidx,]$eeg_pc1 <- unname(unlist(s.pca$x[,1]))

allGamma$eeg_pc2 <- NA
allGamma[completeidx,]$eeg_pc2 <- unname(unlist(s.pca$x[,2]))

allGamma$eeg_pc3 <- NA
allGamma[completeidx,]$eeg_pc3 <- unname(unlist(s.pca$x[,3]))

allGamma$eeg_pc4 <- NA
allGamma[completeidx,]$eeg_pc4 <- unname(unlist(s.pca$x[,4]))

allGamma$eeg_pc5 <- NA
allGamma[completeidx,]$eeg_pc5 <- unname(unlist(s.pca$x[,5]))

allGamma$eeg_pc6 <- NA
allGamma[completeidx,]$eeg_pc6 <- unname(unlist(s.pca$x[,6]))

return(allGamma)

}