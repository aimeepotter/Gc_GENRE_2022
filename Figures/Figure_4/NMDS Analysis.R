# Read in data
T1_130alone_fluxes <- read.delim('~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/T1_130alone_complete_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
T1_130withPMN_fluxes <- read.delim('~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/T1_130withPMN_complete_maxfit/flux_samples.tsv', sep='\t', header=TRUE, row.names=1)
# Identify context-specific reactions
T1_130alone_only <- setdiff(colnames(T1_130alone_fluxes), colnames(T1_130withPMN_fluxes))
T1_130withPMN_only <- setdiff(colnames(T1_130withPMN_fluxes), colnames(T1_130alone_fluxes))
# Get median absolute flux for each
T1_130alone_med_flux <- c()
for (x in T1_130alone_only) {T1_130alone_med_flux <- c(T1_130alone_med_flux, abs(median(T1_130alone_fluxes[,x])))}
T1_130alone_only_flux <- as.data.frame(cbind(T1_130alone_only, T1_130alone_med_flux, rep('T1_130alone',length(T1_130alone_only))))
colnames(T1_130alone_only_flux) <- c('reaction','abs_med_flux','context')
T1_130withPMN_med_flux <- c()
for (x in T1_130withPMN_only) {T1_130withPMN_med_flux <- c(T1_130withPMN_med_flux, abs(median(T1_130withPMN_fluxes[,x])))}
T1_130withPMN_only_flux <- as.data.frame(cbind(T1_130withPMN_only, T1_130withPMN_med_flux, rep('T1_130withPMN',length(T1_130withPMN_only))))
colnames(T1_130withPMN_only_flux) <- c('reaction','abs_med_flux','context')
unique_fluxes <- as.data.frame(rbind(T1_130alone_only_flux, T1_130withPMN_only_flux))
rm(T1_130alone_only, T1_130alone_med_flux, T1_130alone_only_flux, T1_130withPMN_only, T1_130withPMN_med_flux, T1_130withPMN_only_flux, x)
write.table(unique_fluxes, file='~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/PMN_unique_flux_complete_maxfit_median.tsv',
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
# Biomass flux
T1_130withPMN_biomass <- as.vector(T1_130withPMN_fluxes[,'biomass'])
T1_130alone_biomass <- as.vector(T1_130alone_fluxes[,'biomass'])
biomass_pval <- round(wilcox.test(T1_130withPMN_biomass, T1_130alone_biomass, exact=FALSE)$p.value, 10)
# Format data for intersection analysis
shared_rxns <- intersect(colnames(T1_130alone_fluxes), colnames(T1_130withPMN_fluxes))
T1_130alone_fluxes <- T1_130alone_fluxes[,shared_rxns]
T1_130withPMN_fluxes <- T1_130withPMN_fluxes[,shared_rxns]
rm(shared_rxns)
# Remove Biomass components
T1_130alone_fluxes[,c('biomass')] <- NULL
T1_130withPMN_fluxes[,c('biomass')] <- NULL
# Subsample data
sample_size <- min(c(nrow(T1_130alone_fluxes), nrow(T1_130withPMN_fluxes), 500))
sub_sample <- sample(1:min(c(nrow(T1_130alone_fluxes), nrow(T1_130withPMN_fluxes))), sample_size, replace=FALSE)
T1_130alone_fluxes <- T1_130alone_fluxes[sub_sample,]
T1_130withPMN_fluxes <- T1_130withPMN_fluxes[sub_sample,]
rm(sample_size, sub_sample)
# Add column and merge
condition_fluxes <- as.data.frame(rbind(T1_130alone_fluxes, T1_130withPMN_fluxes))
rownames(condition_fluxes) <- c(paste0('T1_130alone_', c(1:nrow(T1_130alone_fluxes))), paste0('T1_130withPMN_', c(1:nrow(T1_130withPMN_fluxes))))
T1_130alone_names <- paste0('T1_130alone_', c(1:nrow(T1_130alone_fluxes)))
T1_130withPMN_names <- paste0('T1_130withPMN_', c(1:nrow(T1_130withPMN_fluxes)))
rownames(condition_fluxes) <- c(T1_130alone_names, T1_130withPMN_names)
condition_fluxes <- condition_fluxes + abs(min(condition_fluxes))
# Create metadata
T1_130alone_metadata <- cbind(T1_130alone_names, rep('T1_130alone', length(T1_130alone_names)))
T1_130withPMN_metadata <- cbind(T1_130withPMN_names, rep('T1_130withPMN', length(T1_130withPMN_names)))
condition_metadata <- rbind(T1_130alone_metadata, T1_130withPMN_metadata)
colnames(condition_metadata) <- c('label', 'condition')
condition_metadata <- as.data.frame(condition_metadata)
rm(T1_130alone_metadata, T1_130withPMN_metadata)
# Calculate dissimilarity (Bray-Curtis)
library(vegan)
flux_dist <- vegdist(condition_fluxes, method='bray')
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=4, trymax=25, )$points)
# Center points
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_xlim <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_ylim <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01
write.table(flux_nmds, file='~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/PMN_flux_nmds_complete_maxfit.tsv', quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
flux_nmds <- read.delim('~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/PMN_flux_nmds_complete_maxfit.tsv', sep='\t',
                        header=TRUE, row.names=1)
# Subset axes
T1_130alone_names <- paste0('T1_130alone_', c(1:nrow(T1_130alone_fluxes)))
T1_130withPMN_names <- paste0('T1_130withPMN_', c(1:nrow(T1_130withPMN_fluxes)))
T1_130alone_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% T1_130alone_names)
T1_130withPMN_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% T1_130withPMN_names)
rm(T1_130alone_names, T1_130withPMN_names)
# Statistical testing (permANOVA)
test <- merge(x=condition_metadata, y=condition_fluxes, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
condition_pval <- adonis(flux_dist ~ condition, data=test, perm=999, method='bray')
condition_pval <- condition_pval$aov.tab[[6]][1]
rm(test)
# Supervised learning
library(randomForest)
flux_groups <- as.factor(c(rep('T1_130alone', nrow(T1_130alone_fluxes)), rep('T1_130withPMN', nrow(T1_130withPMN_fluxes))))
rf_obj <- randomForest(flux_groups ~ ., data=condition_fluxes, importance=TRUE, err.rate=TRUE, ntree=1500, mtry=20)
rf_obj <- importance(rf_obj, type=1, scale=TRUE)
rf_mda <- as.data.frame(subset(rf_obj, rf_obj > (abs(min(rf_obj)))))
rf_mda$reaction <- rownames(rf_mda)
rm(T1_130alone_fluxes, T1_130withPMN_fluxes, flux_groups)
# Subset for plotting
rf_mda <- rf_mda[order(-rf_mda$MeanDecreaseAccuracy),]
rf_mda <- rf_mda[c(1:20),]
rf_mda <- rf_mda[order(rf_mda$MeanDecreaseAccuracy),]
write.table(rf_mda, file='~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/PMN_flux_mda_complete_maxfit.tsv',
            quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)

# Generate figures
T1_130withPMN_col <- 'lightsteelblue2'
T1_130alone_col <- 'red3'
library(scales)
tiff(file='~/Metabolic_Modeling/Gc_GENRE_2022/Figures/Figure_4/PMN_NMDSFigure_complete_maxfit.tiff', units='in', res=1200,  width=5, height=3.5)
layout(matrix(c(1,2,2), nrow=1, ncol=3, byrow=TRUE))
par(mar=c(5,3,1,1), xpd=FALSE, las=1, mgp=c(2.0,0.8,0), lwd=2)
boxplot(T1_130alone_biomass, at=0.5, xlim=c(0,2), ylab='Simulated Biomass Flux',
        ylim=c(0,2), yaxt='n', col=T1_130alone_col,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1,  outline=TRUE, width = 0.5, notch = TRUE)
boxplot(T1_130withPMN_biomass, at=1.5, add=TRUE,  yaxt='n', col=T1_130withPMN_col, outline=TRUE, cex.lab=1.2,
        boxlwd=2, medlwd=2, staplelwd=2, whisklwd=2, whisklty=1, width = 0.5, notch=TRUE)
axis(side=2, at=seq(0,2,.5), cex.axis=0.9, lwd=1.7)
segments(x0=0.5, y0=1.9, x1=1.5, lwd=2)
text(x=1, y=2, '***', cex=1.5, font=2)
par(xpd=TRUE)
text(x=c(0.5,1.5), y=-.35, labels=c('Gc without \nPMNs','Gc with \nPMNs'), cex=1.1, srt=90)
text(x=-0.7, y=2.6,'A', cex=1.4, font=2)
par(xpd=FALSE)
par(mar=c(3.5,3.5,1,1), las=1, mgp=c(2.5,0.8,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim = c(-0.004, 0.004), ylim=c(-0.004, 0.004), lwd.tick=2,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.8, xaxs='i', yaxs='i')
points(x=T1_130alone_nmds_points$MDS1, y=T1_130alone_nmds_points$MDS2, bg=alpha(T1_130alone_col,0.8), pch=21, cex=1.7)
points(x=T1_130withPMN_nmds_points$MDS1, y=T1_130withPMN_nmds_points$MDS2, bg=alpha(T1_130withPMN_col,0.8), pch=21, cex=1.7)
legend('topleft', legend=c('Gc without PMNs','Gc with PMNs'), bg='white',
       pt.bg=c(T1_130alone_col, T1_130withPMN_col), pch=21, pt.cex=1.5, cex=1, box.lwd=2)
legend('bottomright', as.expression(bquote(paste(italic('p'),'-value = 0.001 ***'))), bty='n', pt.cex=0, cex=0.9)
par(xpd=TRUE)
text(x=-0.045, y=0.04, 'B', cex=1.4, font=2)
par(xpd=FALSE)
dev.off()
