library(sommer)

pheno = read.csv("pheno.csv")
pheno$Family = factor(pheno$Family)
pheno$Year = factor(pheno$Year)

A = as.matrix(read.csv("A.csv",row.names = 1,check.names = F))

#Fit the linear mixed model
modelfit = mmer(Status~Year,
                random=~vs(Family,Gu=A[match(pheno$Family,rownames(A)),match(pheno$Family,rownames(A))]),
                rcov=~units,
                data=pheno)

#Estimate and SE of heritability on the observed scale
vpredict(modelfit,h_o~V1/(V1+V2))

#Observed to underlying liability scale factor
scalefactor = mean(pheno$Status)*(1-mean(pheno$Status))/dnorm(qnorm(1-mean(pheno$Status)))^2

#Estimate and SE of heritability on the underlying liability scale
vpredict(modelfit,h_u~V1/(V1+V2)*scalefactor)

#Estimated breeding values of families
EBV = modelfit$U$`u:Family`$Status

#Genetic trend (%)
mean(EBV[grep("27.",names(EBV))])*100
mean(EBV[grep("28.",names(EBV))])*100
mean(EBV[grep("29.",names(EBV))])*100

#Inbreeding coefficients of cohorts
mean(diag(A)[grep("27.",rownames(A))])-1
mean(diag(A)[grep("28.",rownames(A))])-1
mean(diag(A)[grep("29.",rownames(A))])-1
