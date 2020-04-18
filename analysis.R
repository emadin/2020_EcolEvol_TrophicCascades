#### Analyses for Madin et al. (2020) Ecology and Evolution

### Load libaries
library(nlme)
library(plyr)
library(plm)
library(lme4)
library(ggplot2)
library(reshape)
library(stringr)
library(maps)
library(mapdata)
library(oz)
library(Hmisc)


### LOAD DATA

store <- read.csv("data/store.csv", header = TRUE)

### CREATE FUNCTIONS

# Back-transform z to correlation coef r
r <- function(z) {
  (exp(2*z)-1) / (exp(2*z)+1)
}


### SUBSET DATA 

## For alternating trophic trend analysis 
# Note: uses standard "carn" (carnivore) category to integrate different 
# predators of interest in temperate vs tropical datasets

# Keep subset of store with only significant carnivore changes
store_pt <- store[store$p_sig < 0.05,]
store_pt$PART <- factor(store_pt$PART)


### GENERATE MAP OF STUDY LOCATIONS
# where red points = significant, 3-level alternating trophic trends 

pdf("output/figs/Fig_map.pdf", width = 6, height = 10)

  map('worldHires','australia', xlim = c(142,160), ylim = c(-45, -10), fill = TRUE, col = "grey85", border=NA)
  points(store$LONGITUDE[store$ts==0], store$LATITUDE[store$ts==0])
  points(store$LONGITUDE[store$ts==1], store$LATITUDE[store$ts==1], pch = 16, col="red")

  text(store$LONGITUDE[store$LOCATION_DESCRIPTION == "Bicheno External"], store$LATITUDE[store$LOCATION_DESCRIPTION == "Bicheno External"] + 0, "Bicheno External (b)", cex=0.8, pos=4)
  text(store$LONGITUDE[store$LOCATION_DESCRIPTION == "Gannett Cay Reef"], store$LATITUDE[store$LOCATION_DESCRIPTION == "Gannett Cay Reef"], "Gannett Cay Reef (a)", cex=0.8, pos=4)
  text(store$LONGITUDE[store$LOCATION_DESCRIPTION == "Maria Island Reserve"], store$LATITUDE[store$LOCATION_DESCRIPTION == "Maria Island Reserve"] - 0.1, "Maria Island Reserve (c)", cex=0.8, pos=4)
  text(store$LONGITUDE[store$LOCATION_DESCRIPTION == "Tinderbox Reserve"], store$LATITUDE[store$LOCATION_DESCRIPTION == "Tinderbox Reserve"] - 0.4, "Tinderbox Reserve (d)", cex=0.8, pos=4)

  map.axes()
  title(xlab = "Longitude", ylab = "Latitude")
  oz(sections = c(13,15), add = TRUE)
  text(x = c(146.5,146.5,146.5,146.5), y = c(-24.27816, -32.5, -37.3, -42), c("QLD", "NSW", "VIC", "TAS"), cex = 0.8)
  par(fig = c(0.64, 0.795, 0.75, 0.90), mar=c(0,0,0,0), new = TRUE)
  map('worldHires','australia', add = TRUE, fill = TRUE, col="grey85", border=NA)
  # box(which = "figure")
  polygon(x=c(142,154.5,154.5,142), y=c(-9.5,-9.5,-45,-45), border="black")

dev.off()


### MODEL SELECTION

## Cascade occurrence
lat_status_ts <- glm(ts ~ LATITUDE + PART, family=binomial, data=store_pt)
drp_cascade <- drop1(lat_status_ts, test="Chisq")

## Individual trophic group trends
lat_status_pred <- lm(p_z ~ LATITUDE + PART, weights=1/Vz, data=store)
drp_pred <- drop1(lat_status_pred, test="F")

lat_status_herb <- lm(h_z ~ LATITUDE + PART, weights=1/Vz, data=store)
drp_herb <- drop1(lat_status_herb, test="F")

lat_status_alga <- lm(a_z ~ LATITUDE + PART, weights=1/Vz, data=store)
drp_alga <- drop1(lat_status_alga, test="F")

## Output model selection stats
stat_store_cascade <- cbind(stat="lat_status_ts", data.frame(drp_cascade))
stat_store_cascade[,2:5] <- round(stat_store_cascade[,2:5],2) # rounding stats 
stat_store_cascade[,6] <- round(stat_store_cascade[,6],3)
stat_store_trophgrps <- data.frame()
stat_store_trophgrps <- rbind(stat_store_trophgrps, cbind(stat="lat_status_pred", data.frame(drp_pred)))
stat_store_trophgrps <- rbind(stat_store_trophgrps, cbind(stat="lat_status_herb", data.frame(drp_herb)))
stat_store_trophgrps <- rbind(stat_store_trophgrps, cbind(stat="lat_status_alga", data.frame(drp_alga)))
stat_store_trophgrps[,2:6] <- round(stat_store_trophgrps[,2:6],2) # rounding stats 
stat_store_trophgrps[,7] <- round(stat_store_trophgrps[,7],3)

# Generate Table S2: model selection results
write.csv(stat_store_cascade, "output/tables/TableS2_stat_output_modselect_binom_cascade.csv", quote=FALSE, row.names=TRUE)
write.csv(stat_store_trophgrps, "output/tables/TableS2_stat_output_modselect_linear_trophgrps.csv", quote=FALSE, row.names=TRUE)


### LATITUDE ANALYSIS AND PLOTS

reef <- c("a", "b", "c", "d") 

## Generate 4-panel plot of latitude (x-axis) vs. trophic trends (y-axis)
pdf("output/figs/Fig_latitude.pdf", 6, 6)

par(mfrow=c(2,2), mar=c(2, 4, 1, 1), oma=c(3, 1, 1, 1))

# Generalized linear model
lat_ts <- glm(ts ~ LATITUDE, family=binomial, data=store_pt)
summary(lat_ts)
drop1(lat_ts, test="Chisq")

# (Unused) Linear mixed effects model
# lat_ts_mixed <- glmer(ts ~ LATITUDE + (1 | LOCATION_DESCRIPTION), family=binomial, data=store_pt)
# summary(lat_ts_mixed)

plot(jitter(ts, factor=0.1) ~ LATITUDE, data=store_pt[store_pt$ts == 0,], axes=FALSE, xlab="", ylab="Prevalence", ylim=c(0,1))
points(jitter(ts, factor=0.1) ~ LATITUDE, data=store_pt[store_pt$ts == 1,], pch=21, bg=rgb(1,0,0,0.3), col="red", cex=1.2)
mtext("A. Alternating trends", 3, adj=0, cex=1.1)
text(store_pt$LATITUDE[store_pt$LOCATION_DESCRIPTION == "Bicheno External"], 0.95, "b", cex=0.9)
text(store_pt$LATITUDE[store_pt$LOCATION_DESCRIPTION == "Gannett Cay Reef"], 0.95, "a", cex=0.9)
text(store_pt$LATITUDE[store_pt$LOCATION_DESCRIPTION == "Maria Island Reserve"], 0.90, "c", cex=0.9)
text(store_pt$LATITUDE[store_pt$LOCATION_DESCRIPTION == "Tinderbox Reserve"], 0.95, "d", cex=0.9)
axis(1, label=NA)
axis(2, las=2)
LATITUDE <- sort(store_pt$LATITUDE)
pred_ts <- predict(lat_ts, list(LATITUDE), type="response", se.fit=TRUE)
lines(LATITUDE, pred_ts$fit)
polygon(c(LATITUDE, rev(LATITUDE)), c(pred_ts$fit + pred_ts$se.fit, rev(pred_ts$fit - pred_ts$se.fit)), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

lat_pred <- lm(p_z ~ LATITUDE, weights=1/Vz, data=store)
summary(lat_pred)
drop1(lat_pred, test="Chisq")

plot(p_cor ~ LATITUDE, data=store[store$ts == 0,], axes=FALSE, ylab="Standardized trend", xlab="", ylim=c(-1, 1))
points(p_cor ~ LATITUDE, data=store[store$ts == 1,], pch=21, bg=rgb(1,0,0,0.3), col="red", cex=1.2)
mtext("B. Predators", 3, adj=0, cex=1.1)
abline(h=0, lty=2)
axis(1, label=NA)
axis(2, las=2)
LATITUDE <- sort(store$LATITUDE)
pred_lat_pred <- r(predict(lat_pred, list(LATITUDE), interval="confidence"))
lines(LATITUDE, pred_lat_pred[,1])
polygon(c(LATITUDE, rev(LATITUDE)), c(pred_lat_pred[,2], rev(pred_lat_pred[,3])), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

lat_herb <- lm(h_z ~ LATITUDE, weights=1/Vz, data=store)
summary(lat_herb)
drop1(lat_herb, test="Chisq")

plot(h_cor ~ LATITUDE, data=store[store$ts == 0,], axes=FALSE, ylab="Standardized trend", ylim=c(-1, 1))
points(h_cor ~ LATITUDE, data=store[store$ts == 1,], pch=21, bg=rgb(1,0,0,0.3), col="red", cex=1.2)
mtext("C. Herbivores", 3, adj=0, cex=1.1)
abline(h=0, lty=2)
axis(1)
axis(2, las=2)
LATITUDE <- sort(store$LATITUDE)
pred_lat_herb <- r(predict(lat_herb, list(LATITUDE), interval="confidence"))
lines(LATITUDE, pred_lat_herb[,1])
polygon(c(LATITUDE, rev(LATITUDE)), c(pred_lat_herb[,2], rev(pred_lat_herb[,3])), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

lat_alga <- lm(a_z ~ LATITUDE, weights=1/Vz, data=store)
summary(lat_alga)
drop1(lat_alga, test="Chisq")

plot(a_cor ~ LATITUDE, data=store[store$ts == 0,], axes=FALSE, ylab="Standardized trend", ylim=c(-1, 1))
points(a_cor ~ LATITUDE, data=store[store$ts == 1,], pch=21, bg=rgb(1,0,0,0.3), col="red", cex=1.2)
mtext("D. Algae", 3, adj=0, cex=1.1)
abline(h=0, lty=2)
axis(1)
axis(2, las=2)
LATITUDE <- sort(store$LATITUDE)
pred_lat_alga <- r(predict(lat_alga, list(LATITUDE), interval="confidence"))
lines(LATITUDE, pred_lat_alga[,1])
polygon(c(LATITUDE, rev(LATITUDE)), c(pred_lat_alga[,2], rev(pred_lat_alga[,3])), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

mtext("Latitude", 1, outer=TRUE, line=1)

dev.off()


### MARINE RESERVE ANALYSIS AND PLOTS

## Generate 4-panel plot of marine reserve status (x-axis) vs. trophic trends (y-axis)
pdf("output/figs/Fig_status.pdf", 6, 6)

par(mfrow=c(2,2), mar=c(2, 4, 1, 1), oma=c(5, 1, 1, 1))

# cascades and status
ns <- paste0("(", c(summary(store_pt$PART), 0), ")")

sta_ts <- glm(ts ~ PART, family=binomial, data=store_pt)
summary(sta_ts)

pred_ts <- predict(sta_ts, list(PART=c("always fished", "early", "late")), type="response", se.fit=TRUE)

bp <- barplot(c(pred_ts$fit, 0), ylim=c(0, 1), names.arg="", col=c("white", "white", "white"), ylab="Prevalence")
mtext("A. Alternating trends", 3, adj=0, cex=1.1)
# arrows(bp, c(ci[,1], 0), bp, c(ci[,2], 0), code=3, angle=90, length=0.1)
arrows(bp, c(pred_ts$fit + pred_ts$se.fit, 0), bp, c(pred_ts$fit - pred_ts$se.fit, 0), code=3, angle=90, length=0.1)
mtext(ns, 1, at=bp, line=0.0, cex=0.6)

points(bp[1], 0.27, pch=21, bg=rgb(1,0,0,0.5), col="red", cex=1.2)
points(bp[2]-0.1, 0.4, pch=21, bg=rgb(1,0,0,0.5), col="red", cex=1.2)
points(bp[2]+0.1, 0.4, pch=21, bg=rgb(1,0,0,0.5), col="red", cex=1.2)
points(bp[3], 0.55, pch=21, bg=rgb(1,0,0,0.5), col="red", cex=1.2)

text(bp[1], 0.32, "b", cex=0.9)
text(bp[2]-0.1, 0.46, "c", cex=0.9)
text(bp[2]+0.1, 0.46, "d", cex=0.9)
text(bp[3], 0.6, "a", cex=0.9)

nss <- c("Always fished", "Early reserve", "Late reserve", "Newly fished")
ns <- paste0("(", summary(store$PART), ")")

sta_pred <- lm(p_z ~ PART, weights=1/Vz, data=store)
# pred1_mixed <- lmer(p_cor ~ PART + (1 | LOCATION_DESCRIPTION), data=store)
summary(sta_pred)

pred_sta_pred <- r(predict(sta_pred, list(PART=c("always fished", "early", "late", "newly fished")), interval="confidence"))

bp <- barplot(pred_sta_pred[,1], ylim=c(-0.7,0.7), names.arg="", col=c("white", "white", "white", "white"), ylab="Standardized trend")
abline(h=0, lty=1)
arrows(bp, pred_sta_pred[,2], bp, pred_sta_pred[,3], code=3, angle=90, length=0.1)
mtext("B. Predators", 3, adj=0, cex=1.1)
mtext(ns, 1, at=bp, line=0.0, cex=0.6)

sta_herb <- lm(h_z ~ PART, weights=1/Vz, data=store)
summary(sta_herb)
pred_sta_herb <- r(predict(sta_herb, list(PART=c("always fished", "early", "late", "newly fished")), interval="confidence"))
barplot(pred_sta_herb[,1], ylim=c(-0.7,0.7), names.arg="", col=c("white", "white", "white", "white"), ylab="Standardized trend")
abline(h=0, lty=1)
arrows(bp, pred_sta_herb[,2], bp, pred_sta_herb[,3], code=3, angle=90, length=0.1)
mtext("C. Herbivores", 3, adj=0, cex=1.1)
mtext(ns, 1, at=bp, line=0.0, cex=0.6)
mtext(nss, 1, at=bp, line=1.5, las=2, cex=0.8)

sta_alga <- lm(a_z ~ PART, weights=1/Vz, data=store)
summary(sta_alga)
pred_sta_alga <- r(predict(sta_alga, list(PART=c("always fished", "early", "late", "newly fished")), interval="confidence"))
barplot(pred_sta_alga[,1], ylim=c(-0.7,0.7), names.arg="", col=c("white", "white", "white", "white"), ylab="Standardized trend")
abline(h=0, lty=1)
arrows(bp, pred_sta_alga[,2], bp, pred_sta_alga[,3], code=3, angle=90, length=0.1)
mtext("D. Algae", 3, adj=0, cex=1.1)
mtext(ns, 1, at=bp, line=0.0, cex=0.6)
mtext(nss, 1, at=bp, line=1.5, las=2, cex=0.8)

dev.off()


### ALTERNATING TROPHIC TRENDS PLOT

## Generate single-panel plot of trophic groups (x-axis) vs. trophic trends (y-axis)
pdf("output/figs/Fig_cascade.pdf", 5, 4.5)
par(mar=c(3, 4, 1, 1))

plot(0,0, type="n", ylim=c(-1, 1), xlim=c(0.5, 3.5), axes=FALSE, ylab="Standardized trend", xlab="")
axis(2, las=2)
axis(1, at=1:3, labels=c("Predators", "Herbivores", "Algae"))
abline(h=0, lty=2)
for (i in 1:nrow(store_pt)) {
  if (store_pt$ts[i] == 1) {
    points(1:3, store_pt[i,][c("p_cor", "h_cor", "a_cor")], pch=21, bg=rgb(1,0,0,0.5), col="red", cex=1.2)
    lines(1:3, store_pt[i,][c("p_cor", "h_cor", "a_cor")], col="red")
  } else {
    points(1:3, store_pt[i,][c("p_cor", "h_cor", "a_cor")], pch=16, col=rgb(0,0,0,0.3))
    lines(1:3, store_pt[i,][c("p_cor", "h_cor", "a_cor")], col=rgb(0,0,0,0.3))
  }
}
text(0.9, store_pt$p_cor[store_pt$LOCATION_DESCRIPTION == "Bicheno External"], "b", cex=0.9)
text(0.9, store_pt$p_cor[store_pt$LOCATION_DESCRIPTION == "Gannett Cay Reef"], "a", cex=0.9)
text(0.92, store_pt$p_cor[store_pt$LOCATION_DESCRIPTION == "Maria Island Reserve"], "c", cex=0.9)
text(0.85, store_pt$p_cor[store_pt$LOCATION_DESCRIPTION == "Tinderbox Reserve"], "d", cex=0.9)

dev.off()


### TROPHIC FORCING ANALYSIS AND PLOTS

## Generate 2-panel plot of latitude (x-axis) vs. pairwise, between-trophic group trends (y-axis)
pdf("output/figs/Fig_trophicforcing.pdf", width = 4, height = 8) # allows plot size to be adjusted but suppresses plot in RStudio
par(mfrow=c(2,1), mar=c(4, 5, 2, 2))

# Trophic groups' correlations vs. latitude: carnivores-herbivores
mod <- lm(ch_cor ~ LATITUDE, store)
summary(mod) # non-significant, positive relationship; matches expected trend direction
plot(store$LATITUDE, store$ch_cor, axes=FALSE, ylim=c(-1,1), xlab="", ylab="Carnivore-herbivore trend coefficient", col="grey70")
axis(1)
axis(2, las=2)
mtext("A. Predators and herbivores", 3, adj=0, cex=1.1)
abline(mod)
abline(h = 0, lty = 2)
# put SEs around model's fitted values
lat <- sort(store$LATITUDE)
predict_ch <- predict(mod, list(LATITUDE=lat), se.fit=TRUE) 
lines(lat, predict_ch$fit)
polygon(c(lat, rev(lat)), c(predict_ch$fit + predict_ch$se.fit, rev(predict_ch$fit - predict_ch$se.fit)), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

# Trophic groups' correlations vs. latitude: herbivores-algae
mod <- lm(ha_cor ~ LATITUDE, store)
summary(mod) # non-significant, positive relationship; matches expected trend direction
plot(store$LATITUDE, store$ha_cor, axes=FALSE, ylim=c(-1,1), xlab="Latitude", ylab="Herbivore-algae trend coefficient", col="grey70")
axis(1)
axis(2, las=2)
mtext("B. Herbivores and algae", 3, adj=0, cex=1.1)
abline(mod)
abline(h = 0, lty = 2)
# put SEs around model's fitted values
lat <- sort(store$LATITUDE)
predict_ha <- predict(mod, list(LATITUDE=lat), se.fit=TRUE) 
lines(lat, predict_ha$fit)
polygon(c(lat, rev(lat)), c(predict_ha$fit + predict_ha$se.fit, rev(predict_ha$fit - predict_ha$se.fit)), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

dev.off()

### TOP DOWN VS. BOTTOM UP EFFECTS ANALYSIS AND PLOTS

## Address question: where carnivore / algal trends are stronger, are trophic cascades more likely? 
## (with the latter, circumstantially tests whether some envt'l (?) factor - the same factor causing algae to increase over time across most sites - may be responsible for TC occurrence rather than carnivores)

## Analysis of effect of strength of carnivore trend on TC occurrence
# Create subset of store with only significant carn changes
store_pt_carn <- store[store$p_sig < 0.05,]
lat_ts_carn <- glm(ts ~ p_cor, family=binomial, data=store_pt_carn)
summary(lat_ts_carn)
drop1(lat_ts_carn, test="Chisq")

plot(jitter(ts, factor=0.1) ~ p_cor, data=store_pt_carn[store_pt_carn$ts == 0,], axes=FALSE, xlab="Standardized predator trend", ylab="Prevalence", ylim=c(0,1), xlim=c(-1,1)) 
points(jitter(ts, factor=0.1) ~ p_cor, data=store_pt_carn[store_pt_carn$ts == 1,], pch=21, bg=rgb(1,0,0,0.3), col="red", cex=1.2)
#mtext("C. Predators", 3, adj=0, cex=1.1)
text(store_pt_carn$p_cor[store_pt_carn$LOCATION_DESCRIPTION == "Bicheno External"], 0.95, "b", cex=0.9)
text(store_pt_carn$p_cor[store_pt_carn$LOCATION_DESCRIPTION == "Gannett Cay Reef"], 0.95, "a", cex=0.9)
text(store_pt_carn$p_cor[store_pt_carn$LOCATION_DESCRIPTION == "Maria Island Reserve"], 0.90, "c", cex=0.9)
text(store_pt_carn$p_cor[store_pt_carn$LOCATION_DESCRIPTION == "Tinderbox Reserve"], 0.95, "d", cex=0.9)
axis(1)
axis(2, las=2)
PRED_TREND <- sort(store_pt_carn$p_cor)
predict_ts_carn <- predict(lat_ts_carn, list(p_cor=PRED_TREND), type="response", se.fit=TRUE)
lines(PRED_TREND, predict_ts_carn$fit)
polygon(c(PRED_TREND, rev(PRED_TREND)), c(predict_ts_carn$fit + predict_ts_carn$se.fit, rev(predict_ts_carn$fit - predict_ts_carn$se.fit)), border=rgb(0,0,0,0.2), col=rgb(0,0,0,0.2))

dev.copy(pdf,"output/figs/Fig_topdown.pdf")
dev.off()


### GENERATE TABLE S3: MODEL RESULTS (GLM; LMES)

stat_store <- data.frame()
stat_store <- rbind(stat_store, cbind(stat="lat_ts", data.frame(summary(lat_ts)$coefficients)))
stat_store <- rbind(stat_store, cbind(stat="sta_ts", data.frame(summary(sta_ts)$coefficients)))
stat_store[,2:5] <- round(stat_store[,2:5],2) # rounding stats 
stat_store[,5] <- round(stat_store[,5],3)
stat_store2 <- data.frame()
stat_store2 <- rbind(stat_store2, cbind(stat="lat_pred", data.frame(summary(lat_pred)$coefficients)))
stat_store2 <- rbind(stat_store2, cbind(stat="sta_pred", data.frame(summary(sta_pred)$coefficients)))
stat_store2 <- rbind(stat_store2, cbind(stat="lat_herb", data.frame(summary(lat_herb)$coefficients)))
stat_store2 <- rbind(stat_store2, cbind(stat="sta_herb", data.frame(summary(sta_herb)$coefficients)))
stat_store2 <- rbind(stat_store2, cbind(stat="lat_alga", data.frame(summary(lat_alga)$coefficients)))
stat_store2 <- rbind(stat_store2, cbind(stat="sta_alga", data.frame(summary(sta_alga)$coefficients)))
stat_store2[,2:4] <- round(stat_store2[,2:4],2) # rounding stats 
stat_store2[,5] <- round(stat_store2[,5],3)
# Generate Table S3: model results output 
write.csv(stat_store, "output/tables/TableS3_stat_output_binom.csv", quote=FALSE, row.names=TRUE)
write.csv(stat_store2, "output/tables/TableS3_stat_output_linear.csv", quote=FALSE, row.names=TRUE)


