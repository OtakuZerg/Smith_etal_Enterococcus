
solo_col <- 'black'
coculture_col <- '#FE2900'
shared_col <- 'white'
solo_fluxsamples_file <- '~/Desktop/repos/Smith_etal_Enterococcus/data/solo_bhi_new/flux_samples.tsv'
coculture_fluxsamples_file <- '~/Desktop/repos/Smith_etal_Enterococcus/data/coculture_bhi_new/flux_samples.tsv'
solo_transcription_file <- '/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/data/cdiff_only.tsv'
coculture_transcription_file <- '/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/data/coculture.tsv'

library(vegan)
library(ape)
library(vioplot)
library(AUCRF)
library(scales)
library(plotrix)

#-------------------------------------------------------------------------------------------------------------#

# Read in transcriptomic data
solo_transcription <- read.delim(solo_transcription_file, sep='\t', header=TRUE, row.names=1)
coculture_transcription <- read.delim(coculture_transcription_file, sep='\t', header=TRUE, row.names=1)

# Subsample data
sample_size <- round(min(c(colSums(solo_transcription), colSums(coculture_transcription))) * 0.9)
solo_transcription$norm_1 <- as.vector(rrarefy(solo_transcription$norm_1, sample=sample_size))
solo_transcription$norm_2 <- as.vector(rrarefy(solo_transcription$norm_2, sample=sample_size))
solo_transcription$norm_3 <- as.vector(rrarefy(solo_transcription$norm_3, sample=sample_size))
coculture_transcription$norm_1 <- as.vector(rrarefy(coculture_transcription$norm_1, sample=sample_size))
coculture_transcription$norm_2 <- as.vector(rrarefy(coculture_transcription$norm_2, sample=sample_size))
coculture_transcription$norm_3 <- as.vector(rrarefy(coculture_transcription$norm_3, sample=sample_size))
rm(sample_size)

# Calculate medians
solo_transcription$solo_median <- apply(solo_transcription, 1, median)
solo_transcription$norm_1 <- NULL
solo_transcription$norm_2 <- NULL
solo_transcription$norm_3 <- NULL
coculture_transcription$coculture_median <- apply(coculture_transcription, 1, median)
coculture_transcription$norm_1 <- NULL
coculture_transcription$norm_2 <- NULL
coculture_transcription$norm_3 <- NULL
transcription <- merge(solo_transcription, coculture_transcription, by='row.names')
rownames(transcription) <- transcription$Row.names
transcription$Row.names <- NULL

# Save data table
write.table(solo_transcription, file='/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/solo_norm.tsv', quote=FALSE, sep='\t', col.names=FALSE)
write.table(coculture_transcription, file='/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/coculture_norm.tsv', quote=FALSE, sep='\t', col.names=FALSE)
rm(solo_transcription, coculture_transcription)

# Correlate transcriptomes
cor.test(x=transcription$solo_median, y=transcription$coculture_median, method='spearman', exact=FALSE)
# R = 0.968, p << 0.001

# Transform data
transcription$diff <- transcription$solo_median - transcription$coculture_median
transcription$solo_median_log <- log2(transcription$solo_median + 1)
transcription$coculture_median_log <- log2(transcription$coculture_median + 1)

# Subset for enriched genes
solo_enriched <- subset(transcription, diff > 1)
solo_enriched$diff_log <- log2(solo_enriched$diff)
solo_enriched <- subset(solo_enriched, diff_log > 10)
coculture_enriched <- subset(transcription, diff < -1)
coculture_enriched$diff_log <- log2(abs(coculture_enriched$diff))
coculture_enriched <- subset(coculture_enriched, diff_log > 10)

# Generate figure
png(filename='/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/transcriptome_log.png', 
    units='in', width=5, height=4, res=300)
par(mar=c(3,3,0.5,0.5), las=1, mgp=c(1.9,0.75,0), lwd=1.5)
plot(x=transcription$solo_median_log, y=transcription$coculture_median_log, cex=1.1, pch=20,
     xlab='Solo Culture Transcript (Log2)', ylab='Co-culture Transcript (Log2)', 
     xlim=c(0,20), ylim=c(0,20), cex.axis=0.8)
abline(lm(transcription$coculture_median_log ~ transcription$solo_median_log), lwd=3, col='gray')
points(x=solo_enriched$solo_median_log, y=solo_enriched$coculture_median_log, pch=21, cex=1.3, bg=solo_col)
points(x=coculture_enriched$solo_median_log, y=coculture_enriched$coculture_median_log, pch=21, cex=1.3, bg=coculture_col)
text(x=0, y=18, labels='R = 0.968', pos=4)
text(x=0, y=16, labels='p', font=3, pos=4)
text(x=3.2, y=16, labels=' << 0.001***')
legend('bottomright', legend=c('Solo-enriched','Coculture-enriched'), pt.bg=c(solo_col, coculture_col), 
       pch=21, pt.cex=1.5, cex=0.9, box.lwd=2)
box(lwd=2)
dev.off()

rm(transcription, solo_enriched, coculture_enriched)

#-------------------------------------------------------------------------------------------------------------#

# Optimal flux results
#solo_fluxes <- '/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/solo_fluxes.tsv'
#solo_fluxes <- read.delim(solo_fluxes, sep='\t', header=FALSE)
#colnames(solo_fluxes) <- c('reaction','solo_flux')
#coculture_fluxes <- '/home/mjenior/Desktop/repos/Smith_etal_Enterococcus/coculture_fluxes.tsv'
#coculture_fluxes <- read.delim(coculture_fluxes, sep='\t', header=FALSE)
#colnames(coculture_fluxes) <- c('reaction','coculture_flux')
#fluxes <- merge(solo_fluxes, coculture_fluxes, by='reaction')
#rm(solo_fluxes, coculture_fluxes)

#fluxes$abs_diff <- abs(fluxes$solo_flux - fluxes$coculture_flux)
#fluxes <- subset(fluxes, abs_diff > 1.0)

#-------------------------------------------------------------------------------------------------------------#

# Model topology differences

# iCdR700 tology
# genes, rxns, cpds
base <- c(703, 1313, 1243)

# Context-specific models
# solo, coculture, shared
genes <- c(6, 4, 216)
rxns <- c(15, 18, 291)
cpds <- c(7, 10, 299)

# Calculate percentages
genes_perc <- (genes/sum(genes)) * 100.0
rxns_perc <- (rxns/sum(rxns)) * 100.0
cpds_perc <- (cpds/sum(cpds)) * 100.0

# Generate figures  
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/topology.png', 
    units='in', width=4, height=4, res=300)
layout(matrix(c(1,2,
                3,4), nrow=2, ncol=2, byrow=TRUE))
par(mar=c(0.5,0.5,1,0.5), xpd=FALSE, las=1, mgp=c(1.4,0.6,0), lwd=1.7)
pie(rev(genes_perc), labels=c('', '', ''), col=rev(c(solo_col, coculture_col, shared_col)), main='Genes')
pie(rev(rxns_perc), labels=c('', '', ''), col=rev(c(solo_col, coculture_col, shared_col)), main='Reactions')
pie(rev(cpds_perc), labels=c('', '', ''), col=rev(c(solo_col, coculture_col, shared_col)), main='Metabolites')
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
par(mar=c(0,0,0,0), lwd=1.5)
v1 <- expression(italic("C. difficile")*" alone")
v2 <- expression(italic("C. difficile")*" from co-culture")
#v2 <- expression(italic("C. difficile")*" + "*italic("E. faecalis"))
legend('center', legend=c(v1, v2, 'Conserved pathways'), 
       pt.bg=c(solo_col,coculture_col,shared_col), pch=22, pt.cex=1.7, box.lwd=1.5)
dev.off()

#-------------------------------------------------------------------------------------------------------------#

# Read in data
solo_samples <- read.delim(solo_fluxsamples_file, sep='\t', header=TRUE)
solo_samples$X <- NULL
coculture_samples <- read.delim(coculture_fluxsamples_file, sep='\t', header=TRUE)
coculture_samples$X <- NULL

# Format data
overlap <- intersect(colnames(coculture_samples), colnames(solo_samples))
solo_samples <- solo_samples[, overlap]
coculture_samples <- coculture_samples[, overlap]
rm(overlap)

# Subsample data
sample_size <- min(c(nrow(solo_samples), nrow(coculture_samples), 250))
sub_sample <- sample(1:min(c(nrow(solo_samples), nrow(coculture_samples))), sample_size, replace=FALSE)
solo_samples <- solo_samples[sub_sample,]
coculture_samples <- coculture_samples[sub_sample,]
rm(sample_size, sub_sample)

# Format row names
solo_names <- paste('solo_', 1:nrow(solo_samples), sep='')
rownames(solo_samples) <- solo_names
coculture_names <- paste('coculture_', 1:nrow(coculture_samples), sep='')
rownames(coculture_samples) <- coculture_names

# Create metadata
solo_metadata <- cbind(solo_names, rep('solo', length(solo_names)))
coculture_metadata <- cbind(coculture_names, rep('coculture', length(coculture_names)))
metadata <- rbind(solo_metadata, coculture_metadata)
colnames(metadata) <- c('label', 'group')
metadata <- as.data.frame(metadata)
rm(solo_metadata, coculture_metadata)

# Merge data and prep for unsupervised learning
all_samples <- rbind(solo_samples, coculture_samples)
all_samples <- all_samples + abs(min(all_samples))

# Calculate dissimilarity (Bray-Curtis)
flux_dist <- vegdist(all_samples, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
flux_groups <- as.factor(c(rep('solo',nrow(solo_samples)), rep('coculture',nrow(coculture_samples))))
meandist(flux_dist, grouping=flux_groups)
#            coculture  solo
# coculture 0.03563527 0.04101468
# solo      0.04101468 0.03274338

# Unsupervised learning
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=25)$points)
flux_pcoa <- pcoa(flux_dist)
flux_pcoa_axes <- as.data.frame(flux_pcoa$vectors[,c(1,2)])
flux_pcoa_values <- flux_pcoa$values$Relative_eig[1:2]
flux_pcoa_xaxis <- paste0('PC1 (',round(flux_pcoa_values[1],2),'%)')
flux_pcoa_yaxis <- paste0('PC1 (',round(flux_pcoa_values[2],2),'%)')
rm(flux_pcoa, flux_pcoa_values)

# Center points
# NMDS
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01
# PCoA
flux_x <- (abs(max(flux_pcoa_axes$Axis.1)) - abs(min(flux_pcoa_axes$Axis.1))) / 2
flux_y <- (abs(max(flux_pcoa_axes$Axis.2)) - abs(min(flux_pcoa_axes$Axis.2))) / 2
flux_pcoa_axes$Axis.1 <- flux_pcoa_axes$Axis.1 - flux_x
flux_pcoa_axes$Axis.2 <- flux_pcoa_axes$Axis.2 - flux_y
flux_x <- max(abs(max(flux_pcoa_axes$Axis.1)), abs(min(flux_pcoa_axes$Axis.1))) + 0.01
flux_y <- max(abs(max(flux_pcoa_axes$Axis.2)), abs(min(flux_pcoa_axes$Axis.2))) + 0.01

# Subset axes
solo_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% solo_names)
coculture_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% coculture_names)
solo_pcoa_points <- subset(flux_pcoa_axes, rownames(flux_pcoa_axes) %in% solo_names)
coculture_pcoa_points <- subset(flux_pcoa_axes, rownames(flux_pcoa_axes) %in% coculture_names)
rm(solo_names, coculture_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
pval <- adonis(flux_dist ~ group, data=test, perm=99, method='bray')
pval <- round(pval$aov.tab[[6]][1], 3)
if (pval > 0.05) {pval <- 'n.s.'} else {pval <- as.character(pval)}
rm(all_samples, test, flux_dist, metadata)

# Generate figure
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/shared_rxn_nmds.png', 
    units='in', width=5, height=4.5, res=300)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.04, 0.04), ylim=c(-0.04, 0.04),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=solo_nmds_points$MDS1, y=solo_nmds_points$MDS2, bg=alpha(solo_col,0.8), pch=21, cex=1.7)
points(x=coculture_nmds_points$MDS1, y=coculture_nmds_points$MDS2, bg=alpha(coculture_col,0.8), pch=21, cex=1.7)
v1 <- expression(italic("C. difficile")*" alone")
v2 <- expression(italic("C. difficile")*" from co-culture")
#v2 <- expression(italic("C. difficile")*" + "*italic("E. faecalis"))
legend('topright', legend=c(v1, v2), 
       pt.bg=c(alpha(solo_col,0.8),alpha(coculture_col,0.8)), pch=21, pt.cex=1.7, box.lwd=2)
pval <- expression(italic("p")*"-value = 0.01 **")
legend('bottomleft', legend=pval, pt.cex=0, bty='n', cex=0.8)
legend('bottomright', legend='Conserved pathways', pt.cex=0, bty='n')
box(lwd=2)
dev.off()

rm(flux_nmds, solo_nmds_points, coculture_nmds_points, pval, flux_x, flux_y)

#-------------------------------------------------------------------------------------------------------------#

# Limit to Biomass flux
solo_biomass <- solo_samples[,'biomass']
coculture_biomass <- coculture_samples[,'biomass']

# Test differences
pvals <- c()
for (x in c(1:2500)) {
  test_1 <- sample(solo_biomass, size=25)
  test_2 <- sample(coculture_biomass, size=25)
  pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)}
biomass_pval <- median(pvals)
rm(pvals, test_1, test_2)

# Convert to doubling time
solo_doubling <- (1 / solo_biomass) * 3600
coculture_doubling <- (1 / coculture_biomass) * 3600
rm(solo_biomass, coculture_biomass)

# Center on plot area
solo_doubling <- solo_doubling - 30
coculture_doubling <- coculture_doubling - 30

# Generate figure
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/doubling.png', 
    units='in', width=2, height=3, res=300)
par(mar=c(2,2.5,0.5,0.5), xpd=FALSE, las=1, mgp=c(1.4,0.6,0), lwd=1.7)
vioplot(solo_doubling, coculture_doubling, col=c(solo_col,coculture_col),
        ylim=c(0,30), ylab='Predicted Doubling Time (min)', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0,40,5), labels=c(0,seq(40,75,5)), cex.axis=0.6, lwd=1.7)
axis.break(2, 2.5, style='slash', brw=0.03)
segments(x0=1, y0=25, x1=2)
text(x=1.5, y=27, '**', font=2, cex=1.5)
pval <- expression(italic("p")*"-value = 0.006 **")
legend('bottomright', legend=pval, pt.cex=0, bty='n', cex=0.6)
par(xpd=TRUE)
text(x=c(1,2), y=-1.2, labels=c('C. difficile','C. difficile'), cex=0.6, font=3)
text(x=c(1,2), y=-2.6, labels=c('alone','from co-culture'), cex=0.6)
par(xpd=FALSE)
dev.off()

#-------------------------------------------------------------------------------------------------------------#

# Merge data for supervised learning
solo_fluxsamples_file <- '~/Desktop/repos/Smith_etal_Enterococcus/data/old/solo_bhi_reps/flux_samples.tsv'
solo_samples <- read.delim(solo_fluxsamples_file, sep='\t', header=TRUE)
solo_samples$X <- NULL
coculture_fluxsamples_file <- '~/Desktop/repos/Smith_etal_Enterococcus/data/old/coculture_bhi_reps2/flux_samples.tsv'
coculture_samples <- read.delim(coculture_fluxsamples_file, sep='\t', header=TRUE)
coculture_samples$X <- NULL

conserved <- intersect(colnames(solo_samples),colnames(coculture_samples))
keep <- c()
for (x in conserved) {if (!grepl('EX_', x, fixed=TRUE)) {keep <- c(keep, x)}}
solo_samples <- solo_samples[,keep]
coculture_samples <- coculture_samples[,keep]

solo_samples$condition <- 1
coculture_samples$condition <- 0
all_samples <- rbind(solo_samples, coculture_samples)
all_samples$condition <- as.factor(all_samples$condition)
solo_samples$condition <- NULL
coculture_samples$condition <- NULL
  
# Run AUCRF and obtain feature lists
set.seed(906801)
all_aucrf <- AUCRF(condition ~ ., data=all_samples, pdel=0, k0=10)
print(all_aucrf)
rm(all_samples)

# Assemble feature table
top_rxns_importance <- all_aucrf$ranking[1:all_aucrf$Kopt]
top_aucrf <- as.data.frame(cbind(labels(top_rxns_importance), as.vector(top_rxns_importance)))
colnames(top_aucrf) <- c('id','mda')
top_aucrf$mda <- as.numeric(as.character(top_aucrf$mda))
rm(all_aucrf, top_rxns_importance)
write.table(top_aucrf, file='~/Desktop/repos/Smith_etal_Enterococcus/data/all_aucrf_results.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

# Read in previous results
top_aucrf <- read.delim('~/Desktop/repos/Smith_etal_Enterococcus/data/all_aucrf_results.tsv', sep='\t', header=TRUE)
top_aucrf <- top_aucrf[order(top_aucrf$mda),] 
top_aucrf$name <- gsub('_',' ',top_aucrf$name)

# Generate figure
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/aucrf_shared_rxns.png', 
    units='in', width=3, height=4, res=300)
par(mar=c(3, 0.5, 2, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(top_aucrf$mda, bg='gray60', xlim=c(0,max(top_aucrf$mda)*1.1),  
         pch=21, lwd=1.7, pt.cex=1.5, cex=0.8, main='Conserved Pathways')
text(x=-1, y=seq(1.4,10.4,1), labels=top_aucrf$name, cex=0.6, pos=4)
mtext('Mean Decrease Accuracy (%)', side=1, padj=2.5)
dev.off()


#png(filename='~/Desktop/repos/Smith_etal_Enterococcus/.png', 
#    units='in', width=2.5, height=3, res=300)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(1.8,0.7,0), lwd=2)
for (x in c(1:ncol(solo_features))) {
  vioplot(solo_features[,x], coculture_features[,x], col=c(solo_col,coculture_col), main=top_aucrf$label[x], cex.main=1,
          ylim=c(-1000,1000), ylab='Predicted Efflux Rate', lwd=1.7, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=seq(-1000,1000,250), cex.axis=0.7, lwd=2)
  box(lwd=2)
  mtext(c('C. difficile','C. difficile +\nE. faecalis'), side=1, padj=c(0.5, 0.6), adj=c(0.2,0.8), cex=0.7, font=3)
}


# Focus on exchange reactions (All and amino acids alone)
iCdR700_exchs <- c('EX_cpd01030_e','EX_cpd00276_e','EX_cpd00082_e','EX_cpd00080_e','EX_cpd00001_e','EX_cpd00023_e','EX_cpd00971_e','EX_cpd00314_e','EX_cpd00797_e','EX_cpd00051_e','EX_cpd00064_e','EX_cpd10515_e','EX_cpd00179_e','EX_cpd00009_e','EX_cpd00060_e','EX_cpd00076_e','EX_cpd01171_e','EX_cpd00067_e','EX_cpd00154_e','EX_cpd00637_e','EX_cpd00122_e','EX_cpd00035_e','EX_cpd00309_e','EX_cpd00220_e','EX_cpd00588_e','EX_cpd00117_e','EX_cpd00307_e','EX_cpd00254_e','EX_cpd00205_e','EX_cpd00138_e','EX_cpd00305_e','EX_cpd01242_e','EX_cpd03048_e','EX_cpd00065_e','EX_cpd00149_e','EX_cpd00099_e','EX_cpd00162_e','EX_cpd00395_e','EX_cpd00027_e','EX_cpd00794_e','EX_cpd00105_e','EX_cpd00226_e','EX_cpd00104_e','EX_cpd00007_e','EX_cpd00092_e','EX_cpd00208_e','EX_cpd00158_e','EX_cpd00011_e','EX_cpd00063_e','EX_cpd01080_e','EX_cpd00156_e','EX_cpd00039_e','EX_cpd00393_e','EX_cpd00264_e','EX_cpd00322_e','EX_cpd00030_e','EX_cpd00516_e','EX_cpd00084_e','EX_cpd00644_e','EX_cpd00263_e','EX_cpd00232_e','EX_cpd00492_e','EX_cpd00751_e','EX_cpd05161_e','EX_cpd00107_e','EX_cpd00069_e','EX_cpd00066_e','EX_cpd00129_e','EX_cpd00567_e','EX_cpd29317_e','EX_cpd00033_e','EX_cpd00119_e','EX_cpd00132_e','EX_cpd00054_e','EX_cpd00161_e','EX_cpd00053_e','EX_cpd00041_e','EX_cpd05178_e','EX_cpd01711_e','EX_cpd19585_e','EX_cpd00489_e','EX_cpd01042_e','EX_cpd00430_e','EX_cpd00339_e','EX_C21399_e','EX_cpd03343_e','EX_cpd03170_e','EX_cpd00036_e','EX_cpd00211_e','EX_cpd00029_e','EX_cpd00159_e','EX_cpd00221_e','EX_cpd00141_e','EX_cpd00013_e','EX_cpd00242_e','EX_cpd00133_e','EX_cpd00443_e','EX_cpd00210_e','EX_cpd11592_e','EX_cpd11590_e','EX_cpd11588_e','EX_cpd11586_e','EX_cpd11582_e','EX_cpd11581_e','EX_cpd15551_e','EX_cpd15603_e','EX_cpd15604_e','EX_cpd15605_e','EX_cpd15606_e','EX_cpd29695_e','EX_cpd29694_e','EX_cpd29697_e','EX_cpd29698_e','EX_cpd29696_e','EX_cpd29699_e','EX_cpd29693_e','EX_cpd29690_e','EX_cpd29691_e','EX_cpd29700_e','EX_cpd00320_e','EX_cpd02711_e','EX_cpd01553_e','EX_cpd00589_e','EX_cpd27607_e','EX_cpd00731_e','EX_cpd00072_e','EX_cpd00089_e','EX_cpd00079_e','EX_cpd00108_e','EX_cpd00136_e','EX_cpd03561_e','EX_cpd00106_e','EX_cpd00142_e','EX_cpd00383_e','EX_cpd00047_e','EX_cpd00139_e','EX_cpd00281_e','EX_cpd00165_e','EX_cpd03279_e','EX_cpd00246_e','EX_cpd00311_e','EX_cpd00367_e','EX_cpd00277_e','EX_cpd00182_e','EX_cpd00654_e','EX_cpd00412_e','EX_cpd00438_e','EX_cpd00274_e','EX_cpd00186_e','EX_cpd00424_e','EX_cpd00504_e','EX_cpd00155_e','EX_cpd00042_e','EX_cpd00085_e','EX_cpd00152_e','EX_cpd00100_e','EX_cpd00098_e','EX_cpd00207_e','EX_cpd11574_e','EX_cpd04097_e','EX_cpd01012_e','EX_cpd00531_e','EX_cpd00058_e')
aa_exchs <- c('EX_cpd00023_e','EX_cpd00051_e','EX_cpd00064_e','EX_cpd00060_e','EX_cpd00637_e','EX_cpd00035_e','EX_cpd00117_e','EX_cpd00065_e','EX_cpd00156_e','EX_cpd00039_e','EX_cpd00322_e','EX_cpd00084_e','EX_cpd00107_e','EX_cpd00069_e','EX_cpd00066_e','EX_cpd00129_e','EX_cpd00567_e','EX_cpd29317_e','EX_cpd00033_e','EX_cpd00119_e','EX_cpd00132_e','EX_cpd00054_e','EX_cpd00161_e','EX_cpd00053_e','EX_cpd00041_e','EX_cpd00210_e','EX_cpd11592_e','EX_cpd11590_e','EX_cpd11588_e','EX_cpd11586_e','EX_cpd11582_e','EX_cpd11581_e','EX_cpd15551_e','EX_cpd15603_e','EX_cpd15604_e','EX_cpd15605_e','EX_cpd15606_e','EX_cpd29695_e','EX_cpd29694_e','EX_cpd29697_e','EX_cpd29698_e','EX_cpd29696_e','EX_cpd29699_e','EX_cpd29693_e','EX_cpd29690_e','EX_cpd29691_e','EX_cpd29700_e','EX_cpd00320_e','EX_cpd00731_e','EX_cpd00186_e','EX_cpd00085_e')
overlap <- intersect(colnames(coculture_samples), colnames(solo_samples))
shared_exchs <- intersect(overlap, iCdR700_exchs)
solo_exchs <- solo_samples[, shared_exchs]
coculture_exchs <- coculture_samples[, shared_exchs]
rm(overlap, shared_exchs)
solo_aas <- solo_samples[, intersect(colnames(solo_samples), aa_exchs)]
coculture_aas <- coculture_samples[, intersect(colnames(coculture_samples), aa_exchs)]
solo_shared_aas <- solo_samples[, intersect(colnames(solo_aas), colnames(coculture_aas))]
coculture_shared_aas <- coculture_aas[, intersect(colnames(coculture_aas), colnames(solo_aas))]

# Merge data for supervised learning
solo_exchs$condition <- 1
coculture_exchs$condition <- 0
all_exchs <- rbind(solo_exchs, coculture_exchs)
all_exchs$condition <- as.factor(all_exchs$condition)
rm(solo_exchs, coculture_exchs)
solo_shared_aas$condition <- 1
coculture_shared_aas$condition <- 0
all_aas <- rbind(solo_shared_aas, coculture_shared_aas)
all_aas$condition <- as.factor(all_aas$condition)
rm(solo_shared_aas, coculture_shared_aas)

# Run AUCRF and obtain feature lists
exchs_aucrf <- AUCRF(condition ~ ., data=all_exchs, pdel=0, k0=10)
print(exchs_aucrf)
rm(all_exchs)
top_mda <- exchs_aucrf$ranking[1:exchs_aucrf$Kopt]
exchs_aucrf <- as.data.frame(cbind(labels(top_mda), as.vector(top_mda)))
colnames(exchs_aucrf) <- c('id','mda')
exchs_aucrf$mda <- as.numeric(as.character(exchs_aucrf$mda))
aa_aucrf <- AUCRF(condition ~ ., data=all_aas, pdel=0, k0=10)
print(aa_aucrf)
rm(all_aas)
top_mda <- aa_aucrf$ranking[1:aa_aucrf$Kopt]
aa_aucrf <- as.data.frame(cbind(labels(top_mda), as.vector(top_mda)))
colnames(aa_aucrf) <- c('id','mda')
aa_aucrf$mda <- as.numeric(as.character(aa_aucrf$mda))


exchs_aucrf <- AUCRF(condition ~ ., data=all_exchs, pdel=0, k0=10)
print(exchs_aucrf)
rm(all_exchs)
top_mda <- exchs_aucrf$ranking[1:exchs_aucrf$Kopt]
exchs_aucrf <- as.data.frame(cbind(labels(top_mda), as.vector(top_mda)))
colnames(exchs_aucrf) <- c('id','mda')
exchs_aucrf$mda <- as.numeric(as.character(exchs_aucrf$mda))



rm(top_mda)

# Save results
write.table(exchs_aucrf, file='~/Desktop/repos/Smith_etal_Enterococcus/data/exch_aucrf_results.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
write.table(aa_aucrf, file='~/Desktop/repos/Smith_etal_Enterococcus/data/aa_aucrf_results.tsv', 
            quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)

# Read in previous results
exchs_aucrf <- read.delim('~/Desktop/repos/Smith_etal_Enterococcus/data/exch_aucrf_results.tsv', sep='\t', header=TRUE)
exchs_aucrf <- exchs_aucrf[order(exchs_aucrf$mda),] 
aa_aucrf <- read.delim('~/Desktop/repos/Smith_etal_Enterococcus/data/aa_aucrf_results.tsv', sep='\t', header=TRUE)
aa_aucrf <- aa_aucrf[order(aa_aucrf$mda),] 

# Generate figures
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/aucrf_exchs.png', 
    units='in', width=2.5, height=4, res=300)
par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(exchs_aucrf$mda, bg='gray60', xlim=c(0,max(exchs_aucrf$mda)*1.1),  
         pch=21, lwd=1.7, pt.cex=1.5, cex=0.8)
text(x=-1, y=seq(1.4,10.4,1), labels=exchs_aucrf$substrate, cex=0.7, pos=4)
mtext('Mean Decrease Accuracy (%)', side=1, padj=2.5)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/aucrf_aas.png', 
    units='in', width=2.5, height=4, res=300)
par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=1.7)
dotchart(aa_aucrf$mda, bg='gray60', xlim=c(0,max(aa_aucrf$mda)*1.1),  
         pch=21, lwd=1.7, pt.cex=1.5, cex=0.8)
text(x=-1, y=seq(1.4,10.4,1), labels=aa_aucrf$substrate, cex=0.7, pos=4)
mtext('Mean Decrease Accuracy (%)', side=1, padj=2.5)
dev.off()

#-------------------------------------------------------------------------------------------------------------#

plot_flux <- function(exchange, substrate='test', ymax=1000) {
  solo_flux <- abs(subset(solo_samples[,exchange], solo_samples[,exchange] <= 0))
  coculture_flux <- abs(subset(coculture_samples[,exchange], coculture_samples[,exchange] <= 0))
  
  pvals <- c()
  for (x in c(1:2500)) {
    test_1 <- sample(solo_flux, size=25)
    test_2 <- sample(coculture_flux, size=25)
    pvals[x] <- round(wilcox.test(test_1, test_2, exact=FALSE)$p.value, 3)}
  print(median(pvals))
  
  par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1,
          ylim=c(0, ymax), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=3)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1, add=TRUE,
          ylim=c(0, ymax), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  abline(h=0, lwd=2, lty=2, col='gray25')
  box(lwd=3)
  mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
  mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
}


# Read in data
solo_fluxsamples_file <- '~/Desktop/repos/Smith_etal_Enterococcus/data/old/solo_bhi_reps/flux_samples.tsv'
solo_samples <- read.delim(solo_fluxsamples_file, sep='\t', header=TRUE)
solo_samples$X <- NULL
coculture_fluxsamples_file <- '~/Desktop/repos/Smith_etal_Enterococcus/data/old/coculture_bhi_reps2/flux_samples.tsv'
coculture_samples <- read.delim(coculture_fluxsamples_file, sep='\t', header=TRUE)
coculture_samples$X <- NULL

# Amino acids
aa_exchs <- read.delim('~/Desktop/repos/Smith_etal_Enterococcus/data/iCdR703_aa.tsv', sep='\t', header=TRUE)
solo_aas <- solo_samples[, intersect(colnames(solo_samples), aa_exchs$exch_id)]
coculture_aas <- coculture_samples[, intersect(colnames(coculture_samples), aa_exchs$exch_id)]
solo_shared_aas <- solo_samples[, intersect(colnames(solo_aas), colnames(coculture_aas))]
coculture_shared_aas <- coculture_aas[, intersect(colnames(coculture_aas), colnames(solo_aas))]
shared_aa_exchs <- subset(aa_exchs, exch_id %in% colnames(solo_shared_aas))
solo_only_aas <- setdiff(colnames(solo_aas), colnames(coculture_aas))
coculture_only_aas <- setdiff(colnames(coculture_aas), colnames(solo_aas))
rm(solo_aas, coculture_aas)


# Shared AA exchanges
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Serine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00054_e', 'L-Serine', ymax=480)
segments(x0=1, y0=400, x1=2, lwd=3)
text(x=1.5, y=420, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Histidine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00119_e', 'L-Histidine', ymax=4)
segments(x0=1, y0=3, x1=2, lwd=3)
text(x=1.5, y=3.25, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Phenylalanine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00066_e', 'L-Phenylalanine', ymax=12)
segments(x0=1, y0=10, x1=2, lwd=3)
text(x=1.5, y=11, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Tyrosine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00069_e', 'L-Tyrosine', ymax=12)
segments(x0=1, y0=10, x1=2, lwd=3)
text(x=1.5, y=11, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Leucine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00107_e', 'L-Leucine', ymax=24)
segments(x0=1, y0=21, x1=2, lwd=3)
text(x=1.5, y=22, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Cysteine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00084_e', 'L-Cysteine', ymax=4)
segments(x0=1, y0=3, x1=2, lwd=3)
text(x=1.5, y=3.25, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Isoleucine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00322_e', 'L-Isoleucine', ymax=1000)
segments(x0=1, y0=850, x1=2, lwd=3)
text(x=1.5, y=890, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Valine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00156_e', 'L-Valine', ymax=16)
segments(x0=1, y0=14, x1=2, lwd=3)
text(x=1.5, y=15, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Tryptophan.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00065_e', 'L-Tryptophan', ymax=2)
segments(x0=1, y0=1.5, x1=2, lwd=3)
text(x=1.5, y=1.6, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/Ornithine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00064_e', 'Ornithine', ymax=700)
segments(x0=1, y0=625, x1=2, lwd=3)
text(x=1.5, y=660, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Arginine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00051_e', 'L-Arginine', ymax=8)
segments(x0=1, y0=7, x1=2, lwd=3)
text(x=1.5, y=7.5, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Aspartate.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00041_e', 'L-Aspartate', ymax=580)
segments(x0=1, y0=500, x1=2, lwd=3)
text(x=1.5, y=530, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Methionine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00060_e', 'L-Methionine', ymax=8)
segments(x0=1, y0=6, x1=2, lwd=3)
text(x=1.5, y=6.5, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Alanine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00035_e', 'L-Alanine', ymax=1200)
segments(x0=1, y0=1050, x1=2, lwd=3)
text(x=1.5, y=1120, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/D-Alanine.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00117_e', 'D-Alanine', ymax=120)
segments(x0=1, y0=110, x1=2, lwd=3)
text(x=1.5, y=115, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/D-Glutamate.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00186_e', 'D-Glutamate', ymax=80)
segments(x0=1, y0=70, x1=2, lwd=3)
text(x=1.5, y=75, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Glutamate.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00023_e', 'L-Glutamate', ymax=24)
segments(x0=1, y0=21, x1=2, lwd=3)
text(x=1.5, y=22, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/D-Fructose.png', 
    units='in', width=2.5, height=3, res=300)
plot_flux('EX_cpd00082_e', 'D-Fructose', ymax=600)
segments(x0=1, y0=500, x1=2, lwd=3)
text(x=1.5, y=530, '***', font=2, cex=1.3)
dev.off()

# Discordant AA exchanges
coculture_only_flux <- abs(subset(coculture_samples[,'EX_cpd00132_e'], coculture_samples[,'EX_cpd00132_e'] <= 0))
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Asparagine.png',
    units='in', width=2.5, height=3, res=300)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='L-Asparagine', cex.main=1,
        ylim=c(0, 20), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0, 20, 5), cex.axis=0.7, lwd=3)
mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=2, 'inactive', cex=0.9)
dev.off()

coculture_only_flux <- abs(subset(coculture_samples[,'EX_cpd00039_e'], coculture_samples[,'EX_cpd00039_e'] <= 0))
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Lysine.png',
    units='in', width=2.5, height=3, res=300)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='L-Lysine', cex.main=1,
        ylim=c(0, 60), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0, 60, 15), cex.axis=0.7, lwd=3)
mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=3, 'inactive', cex=0.9)
dev.off()

coculture_only_flux <- abs(subset(coculture_samples[,'EX_cpd00161_e'], coculture_samples[,'EX_cpd00161_e'] <= 0))
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/L-Threonine.png',
    units='in', width=2.5, height=3, res=300)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='L-Threonine', cex.main=1,
        ylim=c(0, 15), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0, 15, 3), cex.axis=0.7, lwd=3)
mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=1, 'inactive', cex=0.9)
dev.off()


#------------------------------------------#

# Efflux
plot_efflux <- function(exchange, substrate='test', ymax=1000) {

  solo_flux <- subset(solo_samples[,exchange], solo_samples[,exchange] >= 0)
  coculture_flux <- subset(coculture_samples[,exchange], coculture_samples[,exchange] >= 0)

  pval <- round(wilcox.test(solo_flux, coculture_flux, exact=FALSE)$p.value, 3)
  print(pval)
  
  par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1,
          ylim=c(0, ymax), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=3)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1, add=TRUE,
          ylim=c(0, ymax), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  abline(h=0, lwd=2, lty=2, col='gray25')
  box(lwd=3)
  mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
  mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
}

# Shared
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/N-Acetyl-D-glucosamine.png', 
    units='in', width=2.5, height=3, res=300)
plot_efflux('EX_cpd00122_e', 'N-Acetyl-D-glucosamine', ymax=1200)
segments(x0=1, y0=1050, x1=2, lwd=3)
text(x=1.5, y=1100, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/Acetate.png', 
    units='in', width=2.5, height=3, res=300)
plot_efflux('EX_cpd00029_e', 'Acetate', ymax=1200)
segments(x0=1, y0=1050, x1=2, lwd=3)
text(x=1.5, y=1120, '***', font=2, cex=1.3)
dev.off()

png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/Formate.png', 
    units='in', width=2.5, height=3, res=300)
plot_efflux('EX_cpd00047_e', 'Formate', ymax=600)
segments(x0=1, y0=530, x1=2, lwd=3)
text(x=1.5, y=560, '**', font=2, cex=1.3)
dev.off()

# Discordant
coculture_flux <- subset(coculture_samples[,'EX_cpd19585_e'], coculture_samples[,'EX_cpd19585_e'] >= 0)
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/2-Methylbutyrate.png',
    units='in', width=2.5, height=3, res=300)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='2-Methylbutyrate', cex.main=1,
        ylim=c(0, 15), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0, 15, 3), cex.axis=0.7, lwd=3)
mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=1, 'inactive', cex=0.9)
dev.off()

coculture_only_flux <- subset(coculture_samples[,'EX_cpd00339_e'], coculture_samples[,'EX_cpd00339_e'] >= 0)
png(filename='~/Desktop/repos/Smith_etal_Enterococcus/results/5-Aminovalerate.png',
    units='in', width=2.5, height=3, res=300)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=3)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='5-Aminovalerate', cex.main=1,
        ylim=c(0,1000), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0,1000,250), cex.axis=0.7, lwd=3)
mtext(c('C. difficile','C. difficile'), side=1, padj=0.1, adj=c(0.16,0.86), cex=0.8, font=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=70, 'inactive', cex=0.9)
dev.off()


