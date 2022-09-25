library(mclust)

Fitbit <- read.csv("C:/Users/zxu11/OneDrive/Desktop/JHU Paper/sensi_analysis.csv", header=T,row.names = 1)


# BIC results
BIC <- mclustBIC(Fitbit)
#write.csv(BIC,"C:/Users/zxu11/OneDrive/Desktop/JHU Paper/LPA/BIC_criteria.csv")

## BIC values sensitivity

plot(BIC)
summary(BIC)

# ICL results
ICL <- mclustICL(Fitbit)

plot(ICL)
summary(ICL)

## model fit
mod1 <- Mclust(Fitbit,G=3) # G=2,3,4
## export the raw data of patients in each group
summary(mod1,parameters=TRUE)

### dimension reduction
drmod <- MclustDR(mod1, lambda = 1)
summary(drmod)

## export the coordinates for scatterplot, which subject in which cluster
grp_label <- mod1[["classification"]]
val_dir <- drmod[["dir"]]
#write.csv(val_dir,"C:/Users/zxu11/OneDrive/Desktop/JHU Paper/LPA/value_dir_cluster2.csv")
#write.csv(grp_label,"C:/Users/zxu11/OneDrive/Desktop/JHU Paper/LPA/grp_label_cluster4.csv")

## scatter plot of patients in each direction
plot(drmod,what = 'scatterplot')
legend("topleft",
       legend = c("Cluster1", "Cluster2",'Cluster3'),
       cex=0.7
       )
#plot(drmod, what = "evalues")
