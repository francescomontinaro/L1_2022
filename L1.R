library(ggplot2)
library(plotly)
#This is the first Practical for the Course of Molecular Antropology and Genetics of Dependencies

# Here we will analyse a dataset of more 2300+ individuals from many populations.
# Originally the dataset was genotyped for 2.5 million SNPs using the Illumina Omni 2.5

#As we saw, It's extremely important to perform a strong Quality Control analysis

#First, we will remove Snps wit missing data more than 1%

system("./plink --bfile 1KGPApril2016NewMapAutoChr2 --geno .01 --make-bed --out subset1KGPApril2016NewMapAutoChr2Geno01")

# Then we will remove individuals with more than .1 missingness

system("./plink --bfile subset1KGPApril2016NewMapAutoChr2Geno01 --mind .01 --make-bed --out subset1KGPApril2016NewMapAutoChr2Geno01Missing01")


# Second, we will remove SNPs that are not in Hardy-Weinberg Equilibrium

system("./plink --bfile subset1KGPApril2016NewMapAutoChr2Geno01Missing01 --hwe 1e-6 --make-bed --out subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWE")

#
system("./plink --bfile subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWE --indep-pairwise 250 25 .4 --out subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWELD")

system("./plink --bfile subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWE --exclude subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWELD.prune.out --make-bed --out subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWELD")


system("./plink --bfile subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWELD --pca 3 --out 1KGP_PCA")


pca=read.table("1KGP_PCA.eigenvec")

head(pca)


colnames(pca) <- c("Pop","Ind","PC1","PC2","PC3")

p1=ggplot(pca)+
  geom_point(aes(x=PC1,y=PC2,col=Pop),pch=21)+
  theme_light()


uniquepop <- unique(pca$Pop)

continents=c("Europe","Europe","Asia","Americas","Asia","Americas","Europe","Americas",
             "Asia","Americas","Europe","Africa","Asia","Asia","Asia","Africa",
             "Americas","Americas","Europe","Asia","Africa")

popContinentDF <- data.frame(Pop=uniquepop,Continent=continents)

pca <- merge(pca,popContinentDF,by="Pop")

p2=ggplot(pca)+
  geom_point(aes(x=PC1,y=PC2,fill=Continent,col=Pop),lwd=3,pch=21,stroke=1)+
  theme_light()

p2

ggplotly(p2+
  theme(legend.position = "none"))



#FST

#fst(Ht-Hs)/Ht

system("./plink --bfile subset1KGPApril2016NewMapAutoChr2Geno01Missing01HWE --freq --family --out 1KGP_PCA")

freq=read.table("1KGP_PCA.frq.strat",header=T)

head(freq)

freqGBR=freq[freq$CLST=="GBR",]

freqCHB=freq[freq$CLST=="TSI",]

#First, we will estimate heterozigosisy in single population

freqGBR$H=freqGBR$MAF*(1-freqGBR$MAF)*2

freqCHB$H=freqCHB$MAF*(1-freqCHB$MAF)*2

Hs=(freqGBR$H+freqGBR$H)/2

AF_tot=(freqCHB$MAF+freqGBR$MAF)/2

Ht=AF_tot*(1-AF_tot)*2

fst <- (Ht-Hs)/Ht

fstDF <- data.frame(CHR=freqGBR$CHR,SNP=freqGBR$SNP,fst=fst)

pairfst= function(p1,p2){
# p1=.3
# p2=.4

Hs=((p1*(1-p1)*2)+(p2*(1-p2)*2))/2
ptot=(p1+p2)/2
Ht=(ptot*(1-ptot)*2)

fst=round((Ht-Hs)/Ht,4)
return (fst)
}


fstGBR_CHB <- pairfst(freqGBR$MAF,freqCHB$MAF)


fstAnnotated <- data.frame(pos=freqGBR$SNP,fst=fstGBR_CHB)

fstAnnotated[order(fstAnnotated$fst,decreasing = T),]

fstAnnotated$pos=gsub(pattern = "2_",replacement = "",fstAnnotated$pos)

#plot(fstAnnotated[fstAnnotated$pos %in% 130e6:140e6,2],type="l")

plot(fstAnnotated$pos,fstAnnotated$fst,pch=19)

fstAnnotated[order(fstAnnotated$fst,decreasing = T),][1:20,]
