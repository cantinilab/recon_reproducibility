library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
hummus_object_f <- args[6]
organism <- args[7]
extend <- as.numeric(args[8])
n_cores <- args[9]
n_cores <- if (n_cores == 'NULL') 1 else as.numeric(n_cores)
path_out <- args[10]
multilayer_f <- args[11]
print(extend)
print(multilayer_f)
print("n_cores of type")
print(class(n_cores))

# Set genome
if (organism == 'hg38'){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    specie <- 'human'
} else if (organism == 'mm10'){
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    specie <- 'mouse'
}

print("Calculating peak-to-gene scores")
# placeholder
enhancers = data.frame(gene=character(),
                       peak=character(),
                       score=numeric(),
                       stringsAsFactors=FALSE)


#get only gene, peaks and score columns
enhancers <- enhancers[, c("gene", "peak", "score")]
colnames(enhancers) <- c("gene", "cre", "score")
print("Saving peak-to-gene results")
# Write
write.csv(
  x = enhancers,
  file = path_out,
  row.names = FALSE,
  quote = FALSE
)
