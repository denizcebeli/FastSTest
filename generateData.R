# Gerekli paketler
#library(fastJT)
#library(MASS)
#library(tidyverse)
#library(GGally)
# Girdi parametreleri
SampleSize=20 # <-- Burayı çizelgeye göre değiştir.
nSNP=100
nMArk=2
REALSNP<-c("SNP:98","SNP:99","SNP:100")
FastJPower<-0
FastSPower<-0
itNum<-10 # 5000 olacak
############################################################
# Yapay verinin üretilmesi 3 düzeyli korelasyonlu 3 SNP için
############################################################
# create the variance covariance matrix and the mean vector
sigma<-rbind(c(1,0.8,0.7), c(0.8,1, 0.9), c(0.7,0.9,1))
mu<-c(10, 5, 2)
# generate the multivariate normal distribution
SNPDataCont<-as.data.frame(mvrnorm(SampleSize, mu=mu, Sigma=sigma))
#################### SIMULATION ###############################################
for (it in 1:itNum)
{
# create categorical dependent vectors
SNPData<-SNPDataCont%>%transmute(SNP1= case_when(V1<quantile(V1,0.33)~1,
 V1<quantile(V1,0.67)~2,
TRUE~3),
 SNP2= case_when(V2<quantile(V2,0.33)~1,
 V2<quantile(V2,0.67)~2,
TRUE~3),
 SNP3= case_when(V3<quantile(V3,0.33)~1,
 V3<quantile(V3,0.67)~2,
TRUE~3))
# Create random marker data for each SNP
# MARK1
SNPData<-SNPData%>%mutate(Mrk1=case_when(SNP1==1~rnorm(SampleSize,0,1), # Buralarda
değişiklik olacak
 SNP1==2~rnorm(SampleSize,0.2,1),
SNP1==3~rnorm(SampleSize,0.4,1)))
# MARK2
SNPData<-SNPData%>%mutate(Mrk2=case_when(SNP1==1~rnorm(SampleSize,0,1),
 SNP1==2~rnorm(SampleSize,0.2,1),
SNP1==3~rnorm(SampleSize,0.4,1)))
# For chi square distribution --> rchisq(SampleSize, 1, ncp = 0) + sabit
# For exponential distribution --> rexp(SampleSize, rate = 1) + sabit
# For t distribution --> rt(SampleSize, 3, ncp) + sabit
# For normal distribution --> rnorm(SampleSize,0,1) + sabit
#boxplot(SNPData$Mrk1~SNPData$SNP1)
# Create 100-3 SNP
y<-rep(sample( 1:3, SampleSize, replace=TRUE, prob=c(1/3, 1/3, 1/3) ),times=nSNP-3)
X<-matrix(y,nrow=SampleSize)
MrkData<-cbind(SNPData$Mrk1,SNPData$Mrk2)
MrkD<-as.matrix(MrkData)
SNPData<-cbind(X,SNPData$SNP1,SNPData$SNP2,SNPData$SNP3)
SNPD<-as.matrix(SNPData)
colnames(MrkD) <- paste0("Mrk:",1:nMArk)
colnames(SNPD) <- paste0("SNP:",1:nSNP)
# Running fastJT and fastSG tests
res1 <- fastJT(Y=MrkD, X=SNPD, outTopN=3)
res2 <- FastSG(Y=MrkD, X=SNPD, outTopN=3)
# Detection of correct SNPs for Mark1
Marker1FastJ<-tibble(Z=cbind(res1$J[,1],res1$XIDs[,1]))
FastJSNP1<-filter(Marker1FastJ,Z[,1]>1.64)
Marker1FastS<-tibble(Z=cbind(res2$Mark[,1],res2$SNP[,1]))
FastSSNP1<-filter(Marker1FastS,Z[,1]>1.64)
FastJPower<-FastJPower+length(intersect(REALSNP, FastJSNP1$Z[,2]))
FastSPower<-FastSPower+length(intersect(REALSNP, FastSSNP1$Z[,2]))
# Detection of correct SNPs for Mark2
Marker2FastJ<-tibble(Z=cbind(res1$J[,2],res1$XIDs[,2]))
FastJSNP2<-filter(Marker2FastJ,Z[,2]>1.64)
Marker2FastS<-tibble(Z=cbind(res2$Mark[,2],res2$SNP[,2]))
FastSSNP2<-filter(Marker2FastS,Z[,2]>1.64)
FastJPower<-FastJPower+length(intersect(REALSNP, FastJSNP2$Z[,2]))
FastSPower<-FastSPower+length(intersect(REALSNP, FastSSNP2$Z[,2]))
}
FastJPower/(dim(res1$J)[1]*dim(res1$J)[2]*itNum)
FastSPower/(dim(res1$J)[1]*dim(res1$J)[2]*itNum)
