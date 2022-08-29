# author: "Elsa Leitao"

library(karyoploteR)


### Scheme chrX
############################################################################################################################

# "karyoploteR" used to plot chrX with its karyotype bands
# image saved as "results/scheme_chrX/chrX.png"
# image used as base for a scheme genes with "confirmed" association with monogenic disorders with 
# "Neurologic" Clinical Synopsis plus SOX3 that has no written Clinical Synopsis for "Mental retardation, X-linked, with isolated growth hormone deficiency"

plotKaryotype(genome="hg38", plot.type=2, chromosomes=c("chrX"))
