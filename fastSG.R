FastSG<-function (Y, X, outTopN = 15L, numThreads = 1L)
{
 #Modified S test (Fast S) for Genome-Wide Association
 Mark<-Y
 Geno<-X
 zMark=0;k=1;
 num_sample <- dim(Mark)[1]
 nMark<-dim(Mark)[2]
 nGeno<-dim(Geno)[2]
 markerNames <- colnames(Mark)
 SNPNames <- colnames(Geno)
 for (i in 1:nMark)
 for (j in 1:nGeno)
 {
 veri <- data.frame(Geno[,j],Mark[,i])
 colnames(veri)<-c("A","B")
 newdata <- veri[order(veri$A,veri$B),]
 newdata$A<-as.factor(newdata$A)
 zMark[k]<-fastS(B~A,newdata,verbose = FALSE)$Z
 k<-k+1
 }
 MarkZ<-NULL
 for (i in 1:nMark)
 {
 L<-(i-1)*nGeno+1
 U<-i*nGeno
 #MarkZ<-cbind(MarkZ,as.matrix(zMark[L:U],1,1000))
 MarkZ<-cbind(MarkZ,as.matrix(zMark[L:U]))
 }
 rownames(MarkZ)<-SNPNames
 BestJMark<-NULL
 BestSNP<-NULL
 for (i in 1:nMark)
 {
 BestSNP <- cbind(BestSNP, names(MarkZ[order(-abs(MarkZ[,i])),][1:outTopN,i]))
 BestJMark <- cbind(BestJMark, MarkZ[order(-abs(MarkZ[,i])),][1:outTopN,i])
 }
 rownames(BestJMark)<-NULL
 class(res) <- "JTGenome"
 res<-list()
 res$Mark<-BestJMark
 res$SNP<-BestSNP
 attr(res, "class") <- "owt"
 invisible(res)
}
