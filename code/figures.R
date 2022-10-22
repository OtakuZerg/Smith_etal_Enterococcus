
# Define flux sampling files
solo_fluxsamples_file <- '../data/solo/flux_samples.tsv'
coculture_fluxsamples_file <- '../data/coculture/flux_samples.tsv'

# Read in data
solo_samples <- read.delim(solo_fluxsamples_file, sep='\t', header=TRUE)
coculture_samples <- read.delim(coculture_fluxsamples_file, sep='\t', header=TRUE)
rm(solo_fluxsamples_file, coculture_fluxsamples_file)

#----------------------------------------------------------------------#

# Biomass flux shift
solo_biomass <- solo_samples[,'biomass']
solo_samples[,'biomass'] <- NULL
coculture_biomass <- coculture_samples[,'biomass']
coculture_biomass[,'biomass'] <- NULL

# Test differences
pval <- wilcox.test(solo_biomass, coculture_biomass, exact=FALSE)$p.value

# Generate figure
library(vioplot)
pdf(file='../results/biomass.pdf', width=2, height=3)
par(mar=c(3.2,2.5,1,1), xpd=FALSE, las=1, mgp=c(1.4,0.6,0), lwd=2)
vioplot(solo_biomass, coculture_biomass, col=c(solo_col,coculture_col), cex.axis=0.8, xaxt='n', yaxt='n',
        ylim=c(0,120), ylab='Sampled Biomass Flux', lwd=2, drawRect=FALSE, yaxs='i')
axis(side=2, at=seq(0,120,20), cex.axis=0.7, lwd=2)
segments(x0=1, y0=110, x1=2)
text(x=1.5, y=114, '***', font=2, cex=1.2)
par(xpd=TRUE)
text(x=c(0.8,1.8), y=-15, labels=c('C. difficile','C. difficile'), cex=0.8, font=3, srt=45)
text(x=c(1,2), y=-19, labels=c('alone','in co-culture'), cex=0.8, srt=45)
par(xpd=FALSE)
dev.off()

#----------------------------------------------------------------------#

# Whole distribution analysis
overlap <- intersect(colnames(coculture_samples), colnames(solo_samples))
solo_only <- setdiff(colnames(solo_samples), colnames(coculture_samples))
coculture_only <- setdiff(colnames(coculture_samples), colnames(solo_samples))
solo_only <- solo_samples[, solo_only]
solo_samples <- solo_samples[, overlap]
coculture_only <- coculture_samples[, coculture_only]
coculture_samples <- coculture_samples[, overlap]
rm(overlap)

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

# Calulate dissimilarity
library(vegan)
flux_dist <- vegdist(all_samples, method='bray') # Bray-Curtis

# Mean wwithin-group dissimilarity
flux_groups <- as.factor(c(rep('solo',nrow(solo_samples)), rep('coculture',nrow(coculture_samples))))
meandist(flux_dist, grouping=flux_groups)
#            coculture  solo
# coculture 0.01419866 0.02881653
# solo      0.02881653 0.02075598

# coculture 0.03181850 0.04214376
# solo      0.04214376 0.03410532

# Unsupervised learning
flux_nmds <- as.data.frame(metaMDS(flux_dist, k=2, trymax=50)$points)

# NMDS
flux_x <- (abs(max(flux_nmds$MDS1)) - abs(min(flux_nmds$MDS1))) / 2
flux_y <- (abs(max(flux_nmds$MDS2)) - abs(min(flux_nmds$MDS2))) / 2
flux_nmds$MDS1 <- flux_nmds$MDS1 - flux_x
flux_nmds$MDS2 <- flux_nmds$MDS2 - flux_y
flux_x <- max(abs(max(flux_nmds$MDS1)), abs(min(flux_nmds$MDS1))) + 0.01
flux_y <- max(abs(max(flux_nmds$MDS2)), abs(min(flux_nmds$MDS2))) + 0.01

# Subset axes
solo_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% solo_names)
coculture_nmds_points <- subset(flux_nmds, rownames(flux_nmds) %in% coculture_names)

# Statistical testing (permANOVA)
test <- merge(x=metadata, y=all_samples, by.x='label', by.y='row.names')
rownames(test) <- test$label
test$label <- NULL
pval <- adonis(flux_dist ~ group, data=test, perm=999, method='bray')
pval <- round(pval$aov.tab[[6]][1], 4)
if (pval > 0.05) {pval <- 'n.s.'} else {pval <- as.character(pval)}
rm(all_samples, test, flux_dist, metadata)

library(scales)
solo_col <- 'black'
coculture_col <- '#FE2900'
shared_col <- 'white'
pdf(file='../results/shared_rxn_nmds.pdf', width=5, height=4.5)
par(mar=c(3.5,3.5,0.5,0.5), las=1, mgp=c(2.2,0.7,0), lwd=2)
plot(x=flux_nmds$MDS1, y=flux_nmds$MDS2, xlim=c(-0.05, 0.05), ylim=c(-0.05, 0.05),
     xlab='NMDS Axis 1', ylab='NMDS Axis 2', pch=19, cex.lab=1.1, cex=0, cex.axis=0.9)
points(x=solo_nmds_points$MDS1, y=solo_nmds_points$MDS2, bg=alpha(solo_col,0.8), pch=21, cex=1.7)
points(x=coculture_nmds_points$MDS1, y=coculture_nmds_points$MDS2, bg=alpha(coculture_col,0.8), pch=21, cex=1.7)
v1 <- expression(italic("C. difficile")*" alone")
v2 <- expression(italic("C. difficile")*" from co-culture")
legend('topright', legend=c(v1, v2), bg='white',
       pt.bg=c(alpha(solo_col,0.8),alpha(coculture_col,0.8)), pch=21, pt.cex=1.7, box.lwd=2)
pval <- expression(italic("p")*"-value = 0.001 ***")
legend('bottomleft', legend=pval, pt.cex=0, bty='n', cex=0.8)
legend('bottomright', legend='Conserved pathways', pt.cex=0, bty='n')
box(lwd=2)
dev.off()

# Save tables for submission
create_table <- function(data1, data2, cols, outFile) {
  data <- as.data.frame(cbind(data1, data2))
  colnames(data) <- cols
  outname <- paste0('../results/tables/',outFile,'.xlsx')
  write.table(data, outname, quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)}

write.table(all_samples, '../results/all_samples.xlsx', quote=FALSE, sep='\t', row.names=TRUE, col.names=TRUE)
create_table(solo_biomass, coculture_biomass, c('R20291_solo_biomass','R20291_coculture_biomass'), 'fig2b')
create_table(flux_nmds[,1], flux_nmds[,2], c('NMDS_x_coord','NMDS_y_coord'), 'fig2c')

exchs_aucrf <- read.delim('../results/exch_aucrf_results.tsv', sep='\t', header=TRUE)
aa_aucrf <- read.delim('../results/aa_aucrf_results.tsv', sep='\t', header=TRUE)
top_aucrf <- read.delim('../results/aucrf_results.tsv', sep='\t', header=TRUE)

reaction <- c('L-Ornithine ammonia-lyase (L-proline forming)','D-Glucose 6-phosphate ketol-isomerase',
              'NA+/H antiporter','D-Glucose 1-phosphate 1,6-phosphomutase','Malonyl-CoA:pyruvate carboxytransferase',
              'L-Aspartate:2-oxoglutarate aminotransferase','L-Alanine:2-oxoglutarate aminotransferase',
              'D-Fructose transport via PEP','Adenosine transport','Aspartate:Na+ symporter')
mda <- c(38,37,29.5,29,25,23,22,21,20,19)
top_aucrf <- as.data.frame(cbind(reaction, mda))
top_aucrf$mda <- as.numeric(top_aucrf$mda)
top_aucrf <- top_aucrf[order(top_aucrf$mda),]

pdf(file='../results/aucrf_shared.pdf', width=3, height=4)
par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1.4, 0.5, 0), xpd=FALSE, lwd=2.5)
dotchart(top_aucrf$mda, bg='gray60', xlim=c(0,40),  
         pch=21, pt.cex=1.6, cex=0.9)
text(x=-1, y=seq(1.4,10.4,1), labels=top_aucrf$reaction, cex=0.7, pos=4)
mtext('Mean Decrease Accuracy (%)', side=1, padj=2.7, cex=1.1)
dev.off()

#----------------------------------------------------------------------#

# Exchange analysis
solo_fluxsamples_file <- '../data/solo/flux_samples.tsv'
coculture_fluxsamples_file <- '../data/coculture/flux_samples.tsv'

# Read in data
solo_samples <- read.delim(solo_fluxsamples_file, sep='\t', header=TRUE)
coculture_samples <- read.delim(coculture_fluxsamples_file, sep='\t', header=TRUE)
rm(solo_fluxsamples_file, coculture_fluxsamples_file)

solo_col <- 'black'
coculture_col <- '#FE2900'
shared_col <- 'white'

plot_influx <- function(exchange, substrate='test', ymax=0, ignore_solo=FALSE, ignore_coculture=FALSE) {
  
  if (ignore_solo == FALSE) {
    solo_flux <- abs(subset(solo_samples[,exchange], solo_samples[,exchange] <= 0))
  } else {
    coculture_flux <- abs(subset(coculture_only[,exchange], coculture_only[,exchange] <= 0))
    solo_flux <- rep(0, length(coculture_flux))
    }
  
  if (ignore_coculture == FALSE) {
    coculture_flux <- abs(subset(coculture_samples[,exchange], coculture_samples[,exchange] <= 0))
  } else {
    solo_flux <- abs(subset(solo_only[,exchange], solo_only[,exchange] <= 0))
    coculture_flux <- rep(0, length(solo_flux))
    }
  
  pval <- wilcox.test(solo_flux, coculture_flux, exact=FALSE)$p.value
  if (ymax == 0) { ymax <- round(max(c(max(solo_flux),max(coculture_flux))) * 1.2)}
  
  par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1,
          ylim=c(0, ymax), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=3)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1, add=TRUE,
          ylim=c(0, ymax), ylab='Predicted Uptake Flux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  abline(h=0, lwd=2.5, lty=2, col='gray25')
  box(lwd=3)
  mtext(c('alone','co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
  
  if (pval <= 0.05) {segments(x0=1, x1=2, y0=ymax*0.9, lwd=2.5)}
  if (pval <= 0.001) {
    text(x=1.5, y=ymax*0.95, '***', font=2, cex=1.2)
    pval <- 10}
  if (pval <= 0.01) {
    text(x=1.5, y=ymax*0.95, '**', font=2, cex=1.2)
    pval <- 10}
  if (pval <= 0.05) {
    text(x=1.5, y=ymax*0.95, '*', font=2, cex=1.2)} 
}

# Shared Uptake exchanges
pdf(file='../results/D-Fructose.pdf', width=2, height=3)
plot_influx('EX_cpd00082_e', 'D-Fructose', ymax=500)
dev.off()

pdf(file='../results/Ornithine.pdf', width=2, height=3)
plot_influx('EX_cpd00064_e', 'Ornithine', ymax=700)
dev.off()

pdf(file='../results/L-Aspartate.pdf', width=2, height=3)
plot_influx('EX_cpd00041_e', 'L-Aspartate', ymax=550)
dev.off()

pdf(file='../results/L-Serine.pdf', width=2, height=3)
plot_influx('EX_cpd00054_e', 'L-Serine', ymax=450)
dev.off()

pdf(file='../results/L-Phenylalanine.pdf', width=2, height=3)
plot_influx('EX_cpd00066_e', 'L-Phenylalanine', ymax=10)
dev.off()

pdf(file='../results/L-Tyrosine.pdf', width=2, height=3)
plot_influx('EX_cpd00069_e', 'L-Tyrosine', ymax=10)
dev.off()

pdf(file='../results/L-Leucine.pdf', width=2, height=3)
plot_influx('EX_cpd00107_e', 'L-Leucine', ymax=25)
dev.off()

pdf(file='../results/L-Cysteine.pdf', width=2, height=3)
plot_influx('EX_cpd00084_e', 'L-Cysteine', ymax=3)
dev.off()

pdf(file='../results/L-Isoleucine.pdf', width=2, height=3)
plot_influx('EX_cpd00322_e', 'L-Isoleucine', ymax=1000)
dev.off()

pdf(file='../results/L-Valine.pdf', width=2, height=3)
plot_influx('EX_cpd00156_e', 'L-Valine', ymax=15)
dev.off()

pdf(file='../results/L-Tryptophan.pdf', width=2, height=3)
plot_influx('EX_cpd00065_e', 'L-Tryptophan', ymax=2)
dev.off()

pdf(file='../results/L-Arginine.pdf', width=2, height=3)
plot_influx('EX_cpd00051_e', 'L-Arginine', ymax=10)
dev.off()

pdf(file='../results/L-Methionine.pdf', width=2, height=3)
plot_influx('EX_cpd00060_e', 'L-Methionine')
dev.off()

pdf(file='../results/L-Glutamate.pdf', width=2, height=3)
plot_influx('EX_cpd00023_e', 'L-Glutamate')
dev.off()

pdf(file='../results/D-Glutamate.pdf', width=2, height=3)
plot_influx('EX_cpd00186_e', 'D-Glutamate', ymax=70)
dev.off()

pdf(file='../results/L-Alanine.pdf', width=2, height=3)
plot_influx('EX_cpd00035_e', 'L-Alanine')
dev.off()

# Discordant
pdf(file='../results/L-Asparagine.pdf', width=2, height=3)
plot_influx('EX_cpd00132_e', 'L-Asparagine', ignore_solo=TRUE)
text(x=1, y=(12)*0.083, 'inactive')
dev.off()

pdf(file='../results/L-Lysine.pdf', width=2, height=3)
plot_influx('EX_cpd00039_e', 'L-Lysine', ignore_solo=TRUE)
text(x=1, y=(57)*0.083, 'inactive')
dev.off()

pdf(file='../results/L-Threonine.pdf', width=2, height=3)
plot_influx('EX_cpd00161_e', 'L-Threonine', ignore_solo=TRUE)
text(x=1, y=(12)*0.083, 'inactive')
dev.off()


# Shared Efflux exchanges
plot_efflux <- function(exchange, substrate='test', ymax=1000) {
  solo_flux <- subset(solo_samples[,exchange], solo_samples[,exchange] >= 0)
  coculture_flux <- subset(coculture_samples[,exchange], coculture_samples[,exchange] >= 0)
  
  par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1,
          ylim=c(0, ymax), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  axis(side=2, at=seq(0, ymax, ymax/4), cex.axis=0.7, lwd=3)
  vioplot(solo_flux, coculture_flux, col=c(solo_col,coculture_col), main=substrate, cex.main=1, add=TRUE,
          ylim=c(0, ymax), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
  abline(h=0, lwd=2, lty=2, col='gray25')
  box(lwd=3)
  mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
}

pdf(file='../results/N-Acetyl-D-glucosamine.pdf', width=2, height=3)
plot_efflux('EX_cpd00122_e', 'N-Acetyl-D-glucosamine', ymax=1200)
segments(x0=1, y0=1050, x1=2, lwd=2.5)
text(x=1.5, y=1100, '*', font=2, cex=1.3)
dev.off()

pdf(file='../results/Acetate.pdf', width=2, height=3)
plot_efflux('EX_cpd00029_e', 'Acetate', ymax=1200)
segments(x0=1, y0=1050, x1=2, lwd=2.5)
text(x=1.5, y=1120, '***', font=2, cex=1.3)
dev.off()

# Discordant
coculture_only_flux <- abs(subset(coculture_only[,'EX_cpd00339_e'], coculture_only[,'EX_cpd00339_e'] >= 0))
pdf(file='../results/5-Aminovalerate.pdf', width=2, height=3)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='5-Aminovalerate', cex.main=1,
        ylim=c(0, 1200), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0, 1000, 250), cex.axis=0.7, lwd=3)
box(lwd=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=50, 'inactive', cex=0.7)
dev.off()

coculture_only_flux <- abs(subset(coculture_only[,'EX_cpd19585_e'], coculture_only[,'EX_cpd19585_e'] >= 0))
pdf(file='../results/2-Methylbutyrate.pdf', width=2, height=3)
par(mar=c(2,3,1.5,0.5), xpd=FALSE, las=1, mgp=c(2.1,0.7,0), lwd=2)
vioplot(0, coculture_only_flux, col=c(solo_col,coculture_col), main='2-Methylbutyrate', cex.main=1,
        ylim=c(0, 1000), ylab='Predicted Efflux', lwd=2.5, drawRect=FALSE, yaxt='n', yaxs='i')
axis(side=2, at=seq(0, 1000, 250), cex.axis=0.7, lwd=3)
box(lwd=3)
mtext(c('alone',' from co-culture'), side=1, padj=1.5, adj=c(0.2,0.97), cex=0.8)
text(x=1, y=50, 'inactive', cex=0.7)
dev.off()

