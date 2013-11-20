rm(list = ls())

library(OrgMassSpecR)

dn <- IsotopicDistributionN("NECFLQHK", incorp = 0.00368)
dn$t <- "incorp = 0.368%"
d25 <- IsotopicDistributionN("NECFLQHK", incorp = 0.25)
d25$t <- "incorp = 25%"
d50 <- IsotopicDistributionN("NECFLQHK", incorp = 0.5)
d50$t <- "incorp = 50%"
d75 <- IsotopicDistributionN("NECFLQHK", incorp = 0.75)
d75$t <- "incorp = 75%"
d99 <- IsotopicDistributionN("NECFLQHK", incorp = 0.99)
d99$t <- "incorp = 99%"
x <- rbind(dn, d25, d50, d75, d99)
library(lattice)

png("dist.png", width = 7, height = 5, units = "in", res = 300, bg = "transparent")
print(xyplot(percent ~ mz | t, data = x,
             type = "h", 
             xlab = "m/z", 
             ylab = "intensity (%)",
             main = "Peptide NECFLQHK"))
graphics.off()

png("chrom.png", width = 7, height = 5, units = "in", res = 300, bg = "transparent")
DrawChromatogram(example.chromatogram.multiple$time, 
                 example.chromatogram.multiple$intensity / 1000, 
                 range = list(start = c(21.5, 21.925, 23.1, 25.5, 27.35), 
                              stop = c(21.925, 22.4, 23.6, 26.2, 28.0)),
                 color = hcl(h = seq(0, 360, length = 6)[-6], c = 20, l = 90, fixup = FALSE), 
                 xlab = "retention time (min)", 
                 ylab = "intensity x 1000 (cps)", 
                 main = "Compound Mixture")
graphics.off()

t <- FragmentPeptide("NIDALSGMEGR")   # generate theoretical fragment ions
#dev.new(width = 10, height = 5.5)
png("spec.png", width = 7, height = 5, units = "in", res = 300, bg = "transparent")
PeptideSpectrum(example.spectrum.peptide, t, label = "CID", xlim = c(100, 1200))
title("Fragmentation of Peptide NIDALSGMEGR")
graphics.off()

png("similarity.png", width = 7, height = 5, units = "in", res = 300, bg = "transparent")
SpectrumSimilarity(example.spectrum.unknown, example.spectrum.authentic,
                   top.label = "unknown", 
                   bottom.label = "derivatized alanine",
                   xlim = c(25, 350))
## label peaks
plot.window(xlim = c(25,350), ylim = c(-125, 125))
text(c(73, 147, 158, 232, 260), c(100, 23, 44, 22, 15) + 10,
     c(73, 147, 158, 232, 260), cex = 0.75)
text(c(73, 147, 158, 232, 260), -c(100, 47, 74, 33, 20) - 10, 
     c(73, 147, 158, 232, 260), cex = 0.75)
title("Spectrum Similarity, GC/MS EI Mode")
graphics.off()

