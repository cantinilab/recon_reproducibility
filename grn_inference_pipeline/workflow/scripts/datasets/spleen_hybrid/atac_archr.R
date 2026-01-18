library(ArchR)

# Parse args
args <- commandArgs(trailingOnly = F)
project_path <- args[6]
matrix_path <- args[7]
var_path <- args[8]
obs_path <- args[9]
print(project_path)
print(matrix_path)
# Load the ArchRProject
project <- loadArchRProject(project_path)

# Extract peak matrix
peak_matrix <- getMatrixFromProject(project, "PeakMatrix")
sparse_peak_matrix <- as(assays(peak_matrix)$PeakMatrix, "TsparseMatrix")
# Save sparse peak matrix for python
writeMM(sparse_peak_matrix, file = matrix_path)

# Extract obs data (cells)
obs <- peak_matrix@colData
# save obs data for python 
write.table(obs, file = obs_path, sep = "\t", quote = F, row.names = T)

# Extract var data (peaks)
var = data.frame(seqnames=seqnames(peak_matrix@rowRanges),
  starts=start(peak_matrix@rowRanges)-1,
  ends=end(peak_matrix@rowRanges),
  names=c(rep(".", length(peak_matrix@rowRanges))),
  scores=c(rep(".", length(peak_matrix@rowRanges))),
  strands=strand(peak_matrix@rowRanges))

var_names <- paste0(var$seqnames, "-", var$starts, "-", var$ends)
write.table(var_names, file = var_path, sep = "\t", quote = F, row.names = F)
