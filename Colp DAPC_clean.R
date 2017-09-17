setwd("C:\\Users\\Tharsis\\Dropbox\\JacobsLab Files\\Clive folder\\Colpichthys\\Colp R")

####DAPC on Colpichthys MSAT data####


library(adegenet)

obj1 <- read.structure("/Users/Tharsis/Dropbox/Colp R/Colp_6-loci_locale.str")
obj1


x <- scaleGen(obj1, NA.method ="mean")
head(x)
dim(x)
grp <- find.clusters(x, max.n.clust = 10)
grp
dapc1 <- dapc(x, grp$grp)
temp <- optim.a.score(dapc1) #optimal: retain 1 PC
temp <- xvalDapc(x, grp$grp) #optimal: retain 20 PC <- this one is better
names(temp)
dapc1
scatter.dapc(dapc1)
posterior <- round(dapc1$posterior,3)
DF1 <- dapc1$ind.coord
hybrids <- DF1[ which(DF1[, 1]  <= 1.5 & DF1[, 1] >= -0.5), ]



####DAPC on Colpichthys morpho data####

obj2 <- read.csv("/Users/Tharsis/Dropbox/JacobsLab Files/Clive folder/Colpichthys/Colp R/non-MSAT_dapc.csv")
head(obj2)
rownames(obj2) <- obj2[, 1] ## set rownames
obj2 <- obj2[, -1]
obj2 <- subset(obj2[ , 3:6]) #subset only the morpho data
obj2 <- na.omit(obj2)
grp2 <- find.clusters(obj2, max.n.clust = 10)
dapc2 <- dapc(obj2, grp = grp2$grp)
temp <- optim.a.score(dapc2)
temp <- xvalDapc(obj2, grp2$grp) # optimal: retain 1 PC
dapc2
scatter.dapc(dapc2, leg=TRUE, txt.leg=paste(c("C. hubbsi", "C. regis")))
title(main="DAPC Microsatellite")
DF1 <- dapc2$ind.coord
hybrids <- DF1[ which(DF1[, 1]  <= 1 & DF1[, 1] >= -1.75), ]

include <- c("CohPRIc1", "CorPRIc1", "CorPRIc8", "CohSGU1", "CohCHY14", "CohSFGb4")
RAG1.het <- subset(DF1, rownames(DF1) %in% include)
y_val <- rep(-.025, 6)
RAG1.het <- RAG1.het[, 1]
RAG1.het <- cbind(RAG1.het, y_val)
scatter.dapc(dapc2, leg=TRUE, txt.leg=paste(c("C. hubbsi", "C. regis")))
points(x = RAG1.het[, 1] , y = RAG1.het[ , 2], pch = 17)
legend(x = 3.25, y = 1.05, legend = c("RAG1 heterozygotes"), pch = 17, cex = 0.7)
title(main = "DAPC Morphology")

##regress the discriminant function of morpho on that of MSAT
DF_MSAT <- dapc1$ind.coord
colnames(DF_MSAT) <- "MSAT"
DF_Morpho <- dapc2$ind.coord
colnames(DF_Morpho) <- "Morpho"
DF.regress <- merge(DF_MSAT, DF_Morpho, by = 0)
regress <- lm( MSAT ~ Morpho, data=DF.regress, x=TRUE )
summary(regress)

plot(DF.regress$MSAT~DF.regress$Morpho, ylab = "microsatellite discriminant function 1", xlab ="morphology discriminant function 1")
abline(regress)

#where are the RAG het?
edge.het <- RAG1.het[1:4, ]
delta.het <- RAG1.het[5:6, ]
plot(DF.regress$MSAT~DF.regress$Morpho, 
     ylab = "discriminant function 1 (microsatellite)", 
     xlab ="discriminant function 1 (morphology)",
     cex = 1.5,
     pch = ifelse(DF.regress$Row.names %in% row.names(RAG1.het), 17, 1),
     col = ifelse(DF.regress$Row.names %in% row.names(edge.het), "#CC79A7",
                  ifelse(DF.regress$Row.names %in% row.names(delta.het), "#0072B2", "black"))
)
abline(regress)
legend(x = 1.35, y = 4, legend = c("RAG1 heterozygotes"), pch = 17)

# color code by locale
location <- substring(DF.regress$Row.names, 4,6)
DF.regress <- cbind(DF.regress, location)
# color code by region -include heterozygotes
locale <- unique(location)
region <- cbind(locale, c("Delta", "Delta", "Delta Edge", "Delta","Delta Edge","Delta", "Sonora Coast", "Baja Coast", "Baja Coast", "Sonora Coast", "Sonora Coast", "Sonora Coast"))
library(data.table)
dt1 <- data.table(DF.regress, key = "location")
dt2 <- data.table(region, key = "locale")
joined <- dt1[dt2]
DF.regress <- as.data.frame(joined)
colnames(DF.regress) <- c("Specimens", "MSAT", "Morpho", "location", "region")
cols_region <- as.factor(DF.regress$region)
palette(c("red", "skyblue", "#E69F00", "#008000"))
par(bg = "white")
plot(DF.regress$MSAT~DF.regress$Morpho, col = cols_region,
     pch = ifelse(DF.regress$Specimens %in% row.names(RAG1.het), 17, 8),
     ylim = c(-6, 3), ylab = "microsatellite discriminant function 1",
     xlab ="morphology discriminant function 1",
     cex = 1.5)
legend("topright", levels(factor(DF.regress$region)) , pch = 15, col = c("red", "skyblue", "#E69F00", "#008000"), cex = 1)
abline(regress)

# using ggplot2
library(ggplot2)
library(ggExtra)
library(ggrepel)
regression.plot <- ggplot(DF.regress, aes(x = Morpho, y = MSAT, color = region, shape = RAG1)) +
  geom_abline(slope = -0.8016, intercept = 0.2946) +
  geom_point(size = 3, stroke = 2) + 
  scale_shape_manual(values = c(2, 16)) +
  scale_color_manual(values=c("#33a02c", "#fb9a99", "#1f78b4", "#DAA520")) +
  labs(x = "Phenotype Discriminant Function 1",
  y = "Microsatellite Discriminant Function 1") #+
  # geom_label_repel(data = subset(DF.regress, 
  #                               location == "TRC" | 
  #                                 location == "SGU" |
  #                                 location == "PRI"), 
  #                 aes(label = location), label.size = 0.05,
  #                 box.padding = unit(0.5, "lines"))
  

pp <- regression.plot + theme(legend.position="left", legend.text=element_text(size = 13),
                              plot.margin = unit(c(1, 1, 2, 1), "cm"))
ggMarginal(pp, type = "histogram", bins = 30)

# Same plot as above, but with location labels instead of points
regression.plot <- ggplot(DF.regress, aes(x = Morpho, y = MSAT, color = region)) +
  geom_abline(slope = -0.8016, intercept = 0.2946) +
  scale_color_manual(values=c("#33a02c", "#fb9a99", "#1f78b4", "#DAA520")) +
  labs(x = "Phenotype Discriminant Function 1",
       y = "Microsatellite Discriminant Function 1") +
  geom_text(data = DF.regress, aes(Morpho,label = location)) 

pp <- regression.plot + theme(legend.position="left", legend.text=element_text(size = 13),
                              plot.margin = unit(c(1, 1, 2, 1), "cm"))
ggMarginal(pp, type = "histogram", bins = 30)



## data frame with RAG1 labels
DF.regress <- read.csv("C:\\Users\\Tharsis\\Dropbox\\JacobsLab Files\\Clive folder\\Colpichthys\\Colp R\\msat-morpho regress.csv")
head(DF.regress)

#############################################################
# Geographic cline across the southwestern Delta border

geog_cline <- subset(DF.regress, location == "TRC" | 
                         location == "SGU" |
                         location == "PRI")

par(mfrow = c(1, 2))
plot(geog_cline$Morpho~droplevels(geog_cline)$location,
     xlab = "Hybrid zone locales", ylab = "Phenotype Discriminant Function 1")
plot(geog_cline$MSAT~droplevels(geog_cline)$location,
     xlab = "Hybrid zone locales", ylab = "Microsatellite Discriminant Function 1")
