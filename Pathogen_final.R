library(data.table)
library(inlabru)
library(INLA)
library(mgcv)
library(tidyverse)
library(ggplot2)

setwd('~/Documents/Postdoc/Disease/')

#### Import Data ####
human.reads <- setDT(read.table("./diseases.krakenuniq_class.summary.tsv", header = T))
DT <- fread("./diseases.all.csv")
DT <- merge(DT, human.reads, by = 'sampleId')
DT <- DT[!is.na(longitude)]

PDT <- melt(DT, measure.vars = c("Yersinia_pestis", "Streptococcus_mutans", "Clostridium_septicum", 
                                 "Borrelia_recurrentis", "Leptospira_interrogans", "sapronotic", "anthroponotic", 
                                 "other", "opportunistic", "zoonotic"), 
            variable.name = c("Taxa"),
            value.name = c("Presence"))
PDT <- PDT[order(sampleId)]
PDT[, Presence_fctr := factor(Presence, levels = c(0, 1))]
colnames(PDT)[which(colnames(PDT) %in% c("ANCE2", "ANCE4", "ANCE6", "ANCE8"))] <- c("EHG", "WHG", "LVN", "CHG")
colnames(PDT) <- sub("ov_", "mobility_", colnames(PDT))

## Plot
countries <- sf::st_as_sf(terra::vect(rnaturalearthdata::countries50))
cols <- c('gold', 'midnightblue')
names(cols) <- as.character(c(1,0))

ggplot() +
  geom_point(data = PDT, aes(x = ageAverage/1000, y = Presence, fill = Presence_fctr), shape = 21) +
  facet_wrap(~Taxa, nrow = 5, ncol = 2) +
  scale_x_reverse(n.breaks = 8) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(legend.position = 'none') +
  xlab('Time (kyrs BP)')

ggplot() +
  geom_density(data = PDT[Presence == 0 & ageAverage < 14000], aes(x = ageAverage/1000, 
                                                                   fill = Presence_fctr), alpha = 0.5) +
  geom_density(data = PDT[Presence == 1 & ageAverage < 14000], aes(x = ageAverage/1000, 
                                                                   fill = Presence_fctr), alpha = 0.5) +
  facet_wrap(~Taxa, nrow = 5, ncol = 2) +
  scale_fill_manual(values = cols) +
  scale_x_reverse(n.breaks = 8) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  xlab('Time (kyrs BP)')

ggplot() +
  geom_sf(data = countries) +
  geom_point(data = PDT[order(Presence)], aes(x = longitude, y = latitude, fill = Presence_fctr), shape = 21, position = 'jitter') +
  scale_fill_manual('Presence', values = cols) +
  facet_wrap(~Taxa, nrow = 5, ncol = 2) +
  coord_sf(xlim = range(PDT$longitude, na.rm = T), ylim = range(PDT$latitude, na.rm = T)) +
  theme_bw() +
  theme(aspect.ratio = 0.2) +
  xlab('Longitude') +
  ylab('Latitude')

#### Process Data ####
## Filter for ≤ 14 ka BP
PDT <- PDT[ageHigh <= 14000]

## Log-transform human reads
PDT$nClassHuman_log <- log(PDT$nClassHuman)
head(PDT, 3)

## Transform (scale) data and split into list with one data table for each taxon
dt.list <- lapply(unique(PDT$Taxa), function(x){
  d <- PDT[Taxa == x]
  d <- cbind(d[, .(sampleId, region_id, Taxa, Presence)],
             apply(d[, !c("sampleId", "region_id", "Taxa", "Presence", 
                          "ClimateData", "Presence_fctr")], 
                   2, scale))
  return(d)
})
names(dt.list) <- unique(PDT$Taxa)

#### Models ####
## Models list
models.list <- c(m1 = ~ ageAverage + latitude + longitude + nClassHuman_log + 
                   SS_mobility_mds2 + EHG + WHG + LVN + CHG + 
                   bio01 + bio12 + Intercept(1),
                 m2 = ~ ageAverage + latitude + longitude + nClassHuman_log + 
                   SS_mobility_mds2 + bio01 + bio12 + 
                   Intercept(1),
                 m3 = ~ ageAverage + latitude + longitude + nClassHuman_log + 
                   EHG + WHG + LVN + CHG + bio01 + bio12 + 
                   Intercept(1),
                 m4 = ~ ageAverage + latitude + longitude + nClassHuman_log + 
                   bio01 + bio12 + Intercept(1),
                 m5 = ~ ageAverage + latitude + longitude + nClassHuman_log +
                   SS_mobility_mds2 + EHG + WHG + LVN + CHG + 
                   Intercept(1),
                 m6 = ~ ageAverage + latitude + longitude + nClassHuman_log +
                   EHG + WHG + LVN + CHG + Intercept(1),
                 m7 = ~ ageAverage + latitude + longitude + nClassHuman_log +
                   SS_mobility_mds2 + Intercept(1),
                 m8 = ~ ageAverage + latitude + longitude + nClassHuman_log +
                   Intercept(1))

## Run models in parallel
c <- parallel::detectCores() - 2
cls <- parallel::makeCluster(c)
parallel::clusterExport(cls, c("dt.list", "models.list"))

# go through list and run the 8 models for each taxon
# this results in a nested list, where each element (taxon) contains the results of the 8 models 
taxa.models <- pbapply::pblapply(dt.list, function(d){
  require(inlabru)
  lapply(models.list, function(m){
    bru(components = m, like(Presence ~ ., data = d, family = 'binomial'))
  })
}, cl = cls)
stopCluster(cls); gc()

# example
names(taxa.models)
summary(taxa.models$Yersinia_pestis$m1)
plot(taxa.models$Yersinia_pestis$m1, "CHG")
pd <- taxa.models$Yersinia_pestis$m1$marginals.fixed # data of the posterior densities
plot(pd$CHG)

#### Results ####
## Delta WAIC
waic <- rbindlist(pbapply::pblapply(1:length(taxa.models), function(x){
  d <- data.table(Taxa = unique(dt.list[[x]]$Taxa),
             Model = as.factor(names(models.list)),
             WAIC = sapply(taxa.models[[x]], function(s){s$waic$waic})) # extract WAIC scores
  d[, DeltaWAIC := WAIC - min(WAIC)] # calculate delta
  return(d)
}))

labels <- c("All", "Climate + Mobility", "Climate + Ancestry", "Only Climate", "Ancestry + Mobility",
            "Only Ancestry", "Only Mobility", "None")

ggplot(data = waic[Taxa %in% c("Yersinia_pestis", "Borrelia_recurrentis", "zoonotic")],
       aes(Model, DeltaWAIC, color = Taxa)) +
  scale_color_manual(values = c("#332288", "#117733", "#AA4499")) +
  geom_point() +
  geom_line(aes(group = 1), linetype = 2, linewidth = 0.1) +
  facet_wrap(~Taxa) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(Delta*"WAIC")) +
  scale_x_discrete(labels = labels)

ggplot(data = waic,
       aes(Model, DeltaWAIC, color = Taxa)) +
  scale_color_manual(values = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", 
                                "#661100", "#CC6677", "#882255", "#AA4499")) +
  geom_point() +
  geom_line(aes(group = 1), linetype = 2, linewidth = 0.1) +
  facet_wrap(~Taxa) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(Delta*"WAIC")) +
  scale_x_discrete(labels = labels)

ggplot(data = waic[Taxa != "Leptospira_interrogans"],
       aes(Model, DeltaWAIC, color = Taxa)) +
  scale_color_manual(values = c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", 
                                "#661100", "#CC6677", "#882255", "#AA4499")) +
  geom_point() +
  geom_line(aes(group = 1), linetype = 2, linewidth = 0.1) +
  facet_wrap(~Taxa) +
  theme_bw() +
  theme(plot.title = element_text(face = 'bold', hjust = 0.5),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab(expression(Delta*"WAIC")) +
  scale_x_discrete(labels = labels)

## Variable importance matrix
varnames <- c("ageAverage", "latitude", "longitude", "nClassHuman_log", "SS_mobility_mds2",
              "EHG", "WHG", "LVN", "CHG", "bio01", "bio12")
varlabels <-  c("Age", "Latitude", "Longitude", "Num. Human Reads", "Mobility Dist.", "EHG",
                "WHG", "LVN", "CHG", "Annual Temp.", "Annual Prec.")
taxa <- gsub(names(taxa.models), pattern = "_", replacement = " ")

# Build table
Disease <- data.table()
Disease[, Taxa := unlist(lapply(names(taxa.models), rep, length(varnames)))]
Disease[, Variable := rep(varnames, length(taxa))]
Disease <- merge(Disease, waic[DeltaWAIC == 0], by = "Taxa", all.x = T)
Disease[, Selected := lapply(.SD, function(x){
  m <- unique(Model)
  t <- unique(Taxa)
  c <- rownames(summary(taxa.models[[t]][[m]])$inla$fixed) # extract summary of model
  ix <- ifelse(x %in% c, "Yes", "No") # record whether variable is in selected model or not
  return(ix)
}), .SDcols = "Variable", by = "Taxa"]
Disease[, Relationship := unlist(pbapply::pblapply(unique(Disease$Taxa), function(t){
  d <- Disease[Taxa == t]
  m <- unique(d$Model)
  s <- summary(taxa.models[[t]][[m]])$inla$fixed
  # for the variables in the selected model, check whether the relationship is positive or negative
  effect <- sapply(d$Variable, function(x){
    if(x %in% rownames(s)){
      bool <- between(0, s[x, "0.025quant"], s[x, "0.975quant"])
      if(bool == T){
        eff <- "None"
      } else {
        eff <- ifelse(0 < s[x, "0.025quant"], "Positive", "Negative")
      }
      return(eff)
    } else {
      return(NA)
    }
  })
  return(effect)
}))]
Disease$Taxa <- factor(Disease$Taxa, 
                       levels = unique(Disease$Taxa), 
                       labels = gsub(pattern = "_", replacement = " ", x = unique(Disease$Taxa)))
Disease$Variable <- factor(Disease$Variable, 
                       levels = varnames, 
                       labels = varlabels)

## Plot
# color1 <- c("#8E1327", "#DDDDDD", "#F89241")
color2 <- c("#F89241", "#DDDDDD", "cyan3")
names(color2) <- c("Negative", "None", "Positive")

p1 <- ggplot(data = Disease[Taxa %in% taxa[1:5]]) +
  geom_point(aes(x = Taxa, y = Variable, color = as.factor(Relationship)), size = 4) +
  scale_color_manual("Relationship", values = color2, na.value = "transparent",
                     breaks =  c("Negative", "Positive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

p2 <- ggplot(data = Disease[Taxa %in% taxa[6:10]]) +
  geom_point(aes(x = Taxa, y = Variable, color = as.factor(Relationship)), size = 4) +
  scale_color_manual("Relationship", values = color2, na.value = "transparent",
                     breaks =  c("Negative", "Positive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

p3 <- ggplot(data = Disease) +
  geom_point(aes(x = Taxa, y = Variable, color = as.factor(Relationship)), size = 4) +
  scale_color_manual("Relationship", values = color2, na.value = "transparent",
                     breaks =  c("Negative", "Positive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

png("./PathogenCovariatesRelationship.png", width = 12, height = 10, units = 'in', res = 330)
p3
dev.off()

#### Different climatic variables ####
## Add Isothermality and Precipitation Seasonality
models.list <- lapply(models.list, function(x){
  formula(gsub(pattern = "bio12", replacement = "bio12 + bio03 + bio15", x = x))
})

## Run models
c <- parallel::detectCores() - 2
cls <- parallel::makeCluster(c)
parallel::clusterExport(cls, c("dt.list", "models.list"))

taxa.models2 <- pbapply::pblapply(dt.list, function(d){
  require(inlabru)
  lapply(models.list, function(m){
    bru(components = m, like(Presence ~ ., data = d, family = 'binomial'))
  })
}, cl = cls)
stopCluster(cls); gc()

## Calculate WAIC
waic2 <- rbindlist(pbapply::pblapply(1:length(taxa.models2), function(x){
  d <- data.table(Taxa = unique(dt.list[[x]]$Taxa),
                  Model = as.factor(names(models.list)),
                  WAIC = sapply(taxa.models2[[x]], function(s){s$waic$waic}))
  d[, DeltaWAIC := WAIC - min(WAIC)]
  return(d)
}))

## Matrix
varnames <- c(varnames, "bio03", "bio15")
varlabels <- c(varlabels, "Isothermality", "Prec. Seasonality")
Clim <- data.table()
Clim[, Taxa := unlist(lapply(names(taxa.models2), rep, length(varnames)))]
Clim[, Variable := rep(varnames, length(taxa))]
Clim <- merge(Clim, waic2[DeltaWAIC == 0], by = "Taxa", all.x = T)
Clim[, Selected := lapply(.SD, function(x){
  m <- unique(Model)
  t <- unique(Taxa)
  c <- rownames(summary(taxa.models2[[t]][[m]])$inla$fixed)
  ix <- ifelse(x %in% c, "Yes", "No")
  return(ix)
}), .SDcols = "Variable", by = "Taxa"]
Clim[, Relationship := unlist(pbapply::pblapply(unique(Clim$Taxa), function(t){
  d <- Clim[Taxa == t]
  m <- unique(d$Model)
  s <- summary(taxa.models2[[t]][[m]])$inla$fixed
  effect <- sapply(d$Variable, function(x){
    if(x %in% rownames(s)){
      bool <- between(0, s[x, "0.025quant"], s[x, "0.975quant"])
      if(bool == T){
        eff <- "None"
      } else {
        eff <- ifelse(0 < s[x, "0.025quant"], "Positive", "Negative")
      }
      return(eff)
    } else {
      return(NA)
    }
  })
  return(effect)
}))]
Clim$Taxa <- factor(Clim$Taxa, 
                       levels = unique(Clim$Taxa), 
                       labels = gsub(pattern = "_", replacement = " ", x = unique(Clim$Taxa)))
Clim$Variable <- factor(Clim$Variable, 
                           levels = varnames, 
                           labels = varlabels)

p1.2 <- ggplot(data = Clim[Taxa %in% taxa[1:5]]) +
  geom_point(aes(x = Taxa, y = Variable, color = as.factor(Relationship)), size = 4) +
  scale_color_manual("Relationship", values = color2, na.value = "transparent",
                     breaks =  c("Negative", "Positive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

p2.2 <- ggplot(data = Clim[Taxa %in% taxa[6:10]]) +
  geom_point(aes(x = Taxa, y = Variable, color = as.factor(Relationship)), size = 4) +
  scale_color_manual("Relationship", values = color2, na.value = "transparent",
                     breaks =  c("Negative", "Positive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

p3.2 <- ggplot(data = Clim) +
  geom_point(aes(x = Taxa, y = Variable, color = as.factor(Relationship)), size = 4) +
  scale_color_manual("Relationship", values = color2, na.value = "transparent",
                     breaks =  c("Negative", "Positive")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

png("./PathogenCovariatesRelationship_Clim.png", width = 12, height = 10, units = 'in', res = 330)
p3.2
dev.off()