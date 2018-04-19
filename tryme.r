
source("relatedness.r") ## calls a bunch of functions, but you care only about one: relat

load("mydat.RData") ## an example data set
## Notes:
## mydat is a data.frame
## each line is an individual
## each column is a gene, ie two consecutive columns stock a genotype
## the first column is the individuals ID (not compulsory)

## The arguments of relat:
## data= like mydat
## id= 0 if there is no ID is in the 1st column, 1 if there is
## li= T/F (calculate or not Li's estimator)
## lynch= T/F (calculate or not Lynch's estimator)
## queller= T/F (calculate or not Queller's estimator)
## ritland= T/F (calculate or not Ritland's estimator)
## wang= T/F (calculate or not Wang's estimator)
  
r <- relat(mydat) ## call relat with its defaults:
head(r)
##                wang       lynch
## 10L-10P  0.04333763  0.06472117
## 10L-10S  0.13513891  0.10101831
## 10L-10Y -0.08301690 -0.10591922
## 10L-10Z -0.11552752 -0.13236066
## 10L-11F -0.16790429 -0.04995266
## 10L-11I  0.48369941  0.35700626

r <- relat(mydat, li=TRUE, lynch=F, queller=T, ritland=T)
head(r)
##                  li        wang     queller     ritland
## 10L-10P -0.02564872  0.04333763  0.05669502  0.09887280
## 10L-10S  0.14283713  0.13513891  0.19223713  0.13297426
## 10L-10Y -0.05962285 -0.08301690 -0.05246626 -0.14238056
## 10L-10Z -0.11758963 -0.11552752 -0.10984387 -0.14364186
## 10L-11F -0.17485357 -0.16790429 -0.15663328 -0.07818908
## 10L-11I  0.45017932  0.48369941  0.44623673  0.35139520

plot(r) ## pretty good correlation between the different estimators
