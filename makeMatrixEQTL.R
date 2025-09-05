library(VariantAnnotation)
library(dplyr)
library(MatrixEQTL)

# Step 1: Read the mapping CSV file
num_mapped_reads <- read.csv("Num_mapped_reads.csv")

# Step 2: List all the VCF files in the directory (assuming they are all in the working directory)
vcf_files <- list.files(pattern = "*.vcf.gz")

# Step 3: Create an empty list to store genotype data
genotypes_list <- list()

# Step 4: Process each VCF file
for (vcf_file in vcf_files) {
  # Read VCF file
  vcf <- readVcf(vcf_file, genome = "hg19")
snp_info <- rowRanges(vcf)

snps_positions <- data.frame(
  snp = names(snp_info),
  chr = as.character(seqnames(snp_info)),
  pos = start(snp_info)
)

write.table(snps_positions, "snps_positions.txt", sep = "\t", row.names = FALSE, quote = FALSE)

  
  # Extract genotype information
  geno <- geno(vcf)$GT
  
  # Convert the genotypes to numeric
  convert_gt_to_numeric <- function(gt) {
    if (gt %in% c("0/0", "0|0")) return(0)
    if (gt %in% c("0/1", "1/0", "0|1", "1|0")) return(1)
    if (gt %in% c("1/1", "1|1")) return(2)
    return(NA)
  }
  
  geno_numeric <- apply(geno, c(1,2), convert_gt_to_numeric)
  geno_matrix <- as.data.frame(t(geno_numeric))
  
  # Step 5: Extract the sample ID from the VCF file name (assuming it matches the 'Run' column in the CSV)
sample_id <- sub("^.*_(SRR[0-9]+).*", "\\1", vcf_file)

tissue_type <- num_mapped_reads %>%
  filter(grepl(sample_id, Run)) %>%
  pull(cell_line)

# Safely assign or skip if no match found
if (length(tissue_type) == 0) {
  warning(paste("No tissue type found for", sample_id, "- skipping"))
  next
}
  # Step 7: Add sample's tissue type as a column to the genotype data
  geno_matrix$tissue_type <- tissue_type
  
  # Step 8: Store the genotype matrix in the list
  genotypes_list[[sample_id]] <- geno_matrix
}
# Step 1: Get all unique SNPs from all samples
library(data.table)

all_snps <- sort(unique(unlist(lapply(genotypes_list, colnames))))
library(Matrix)

# Pre-allocate matrix (sparse for memory efficiency)
sample_ids <- names(genotypes_list)
snp_ids <- all_snps
geno_mat <- Matrix(0, nrow = length(sample_ids), ncol = length(snp_ids),
                   dimnames = list(sample_ids, snp_ids), sparse = TRUE)

# Fill in matrix using row + column indexing
for (i in seq_along(sample_ids)) {
  sample <- genotypes_list[[i]]
  sample_snps <- intersect(setdiff(colnames(sample), "tissue_type"), snp_ids)
  geno_mat[i, sample_snps] <- as.numeric(sample[, sample_snps])
}
# Extract tissue types in the same order as sample_ids
tissue_types <- sapply(genotypes_list, function(df) df$tissue_type)
# Set tissue types as rownames of geno_mat
rownames(geno_mat) <- tissue_types

# Load the gene expression data
gene_expression <- read.delim("CCLE_RNAseq_genes_rpkm_20180929.gct", skip = 2)

# Step 1: Save gene names as row names
rownames(gene_expression) <- make.unique(as.character(gene_expression$Description))

# Step 2: Drop unnecessary columns
gene_expression$Name <- NULL
gene_expression$Description <- NULL

gene_expression <- log2(gene_expression + 1)
# Step 3: Ensure only numeric columns remain
gene_expression <- gene_expression[, sapply(gene_expression, is.numeric)]

# install.packages("biomaRt")  # if you haven't already
library(biomaRt)

# Read your gene list (from expression matrix)
genes <- row.names(gene_expression)  # or expr_data$gene, depending on your format

# Use biomaRt to get gene positions
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_positions <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

# Clean up and write to file
colnames(gene_positions) <- c("gene", "chromosome", "position")
write.table(gene_positions, "gene_positions.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Ensure the gene expression data is in the required format
gene_expression_matrix <- as.matrix(gene_expression)

expr_samples <- colnames(gene_expression_matrix)
# Get the intersecting sample IDs
shared_samples <- intersect(colnames(gene_expression_matrix), rownames(geno_mat))
# Subset and reorder both matrices
gene_expression_matrix_sub <- gene_expression_matrix[,shared_samples]
geno_mat_sub <- geno_mat[shared_samples, ]
geno_mat_sub <- t(geno_mat_sub)

# Apply MAF filter
calculate_maf <- function(genotypes) {
  valid_genotypes <- genotypes[!is.na(genotypes)]
  if (length(valid_genotypes) == 0) return(NA)
  total_alleles <- 2 * length(valid_genotypes)
  alt_alleles <- sum(valid_genotypes)
  alt_freq <- alt_alleles / total_alleles
  maf <- min(alt_freq, 1 - alt_freq)
  return(maf)
}

maf_values <- apply(geno_mat_sub, 1, calculate_maf)
maf_threshold <- 0.05
snp_keep <- names(maf_values)[which(maf_values > maf_threshold)]
geno_mat_sub <- geno_mat_sub[snp_keep, ]

# Now create SlicedData
snps <- SlicedData$new()
snps$CreateFromMatrix(as.matrix(geno_mat_sub))

# Create SlicedData objects for SNPs, gene expression, and covariates
snps <- SlicedData$new()
snps$CreateFromMatrix(as.matrix(geno_mat_sub))  # as.matrix() will coerce it

gene <- SlicedData$new()
gene$CreateFromMatrix(gene_expression_matrix_sub)
print(identical(rownames(gene_expression_matrix_sub), colnames(geno_mat_sub)))

# If there are covariates, load them (assuming you have a covariates file or data)
gene_variances <- apply(gene_expression_matrix_sub, 1, var)
gene_expression_matrix_filtered <- gene_expression_matrix_sub[gene_variances > 0, ]

# Run PCA on the expression matrix (samples as rows)
expr_pcs <- prcomp(t(gene_expression_matrix_filtered), center = TRUE, scale. = TRUE)

# Get the top N principal components (e.g., top 5)
top_pcs <- expr_pcs$x[, 1:32]
# Remove columns with zero variance
geno_var <- apply(geno_mat_sub, 1, var, na.rm = TRUE)
geno_mat_sub_filtered <- geno_mat_sub[geno_var > 0,]
geno_mat_sub_filtered[is.na(geno_mat_sub_filtered)] <- median(geno_mat_sub_filtered, na.rm = TRUE)

# Now you can safely run PCA
geno_pcs <- prcomp(t(geno_mat_sub_filtered), center = TRUE, scale. = TRUE)
top_geno_pcs<-geno_pcs$x[,1:20]
shared_samples <- intersect(rownames(top_pcs), rownames(top_geno_pcs))

top_pcs <- top_pcs[shared_samples,,drop=FALSE ]
top_geno_pcs <- top_geno_pcs[shared_samples,,drop=FALSE ]
all_covariates <- cbind(top_pcs, top_geno_pcs)



# Create SlicedData object for covariates using only expression PCs
cvrt <- SlicedData$new()
cvrt$CreateFromMatrix(t(all_covariates))

# Find common sample IDs across all

snpspos <- read.table("snps_positions.txt", header = TRUE)
genepos <- read.table("gene_positions.txt", header = TRUE)

# Run MatrixEQTL
output_file_name <- "output_matrixeqtl_results.txt"
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = 0.05,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = NULL,
  pvOutputThreshold.cis = 0,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = 1e6
)
cat("Number of eQTL results:", nrow(me$all$eqtls), "\n")

results <- read.table("output_matrixeqtl_results.txt", header = TRUE)
results$FDR <- p.adjust(results$p.value, method = "fdr")
write.table(results,"output_matrixeqtl_results.txt",sep = "\t", quote = FALSE, row.names = FALSE)
# The result is saved to the output file defined above, e.g. 'output_matrixeqtl_results.txt'
