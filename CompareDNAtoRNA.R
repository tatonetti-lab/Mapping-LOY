library(tidyr)
library(dplyr)
library(ggplot2)

# Read in the GCT file
gct_file <- "RNAseqdataClean.txt"
gct_data <- read.table(gct_file, header = TRUE, sep = "\t") # Skip the first 2 GCT metadata lines

# Extract the gene names (row names) and transpose the data for analysis
gene_positions <- gct_data$Description  # Assuming the "Name" column contains the gene IDs or positions
sample_data <- gct_data[, -1]  # Remove the first two columns ("Name" and "Description")

# Convert the data into long format
long_data <- sample_data %>%
  mutate(Gene = gene_positions) %>%
  pivot_longer(-Gene, names_to = "SampleID", values_to = "Expression")

# Filter for Y chromosome genes
y_gene_list <- read.table("genes_chrY_positions.txt", header = TRUE, stringsAsFactors = FALSE)
names(y_gene_list) <- c("Gene","Start","End")
filtered_data <- long_data %>% filter(Gene %in% y_gene_list$Gene)

# Step 2: Merge with filtered expression data
filtered_data <- filtered_data %>%
  left_join(y_gene_list, by = "Gene") %>%
  arrange(Start) %>% # Order by position
  mutate(Gene = factor(Gene, levels = unique(Gene))) # Ensure genes are ordered by position

sampleOrder<-read.table("orderedSamples.txt",sep="\t",header=TRUE)
run_to_cellline <- read.csv("Num_mapped_reads.csv", header = TRUE)  # Columns: Run, CellLine

sampleOrder <- sampleOrder %>%
  mutate(Run = gsub(".*_(SRR[0-9]+)_.*", "\\1", x)) %>%
  left_join(run_to_cellline, by = "Run") %>%
  distinct(x, .keep_all = TRUE) %>%
  rename(SampleID = x) %>%
  mutate(SampleIndex = dense_rank(SampleID))
filtered_data <- filtered_data %>%
  filter(!is.na(SampleID)) %>%
  mutate(SampleID = factor(SampleID, levels = sampleOrder$SampleID)) %>%
  mutate(SampleIndex = dense_rank(SampleID))  # Ensures sequential order without spaces

filtered_data <- filtered_data %>% filter(Start >= 2650741)

zero_expression_samples <- filtered_data %>%
  group_by(SampleID,SampleIndex) %>%
  summarize(TotalExpression = sum(Expression, na.rm = TRUE),.groups="drop") %>%
  arrange(TotalExpression)

# Ensure all samples from sampleOrder are represented
all_samples <- sampleOrder %>% select(SampleID, SampleIndex)

filtered_data_full <- filtered_data %>%
  right_join(all_samples, by = c("SampleID", "SampleIndex"))
 # Plot the data with black dots for Expression > 0
png("Gene_vs_Sample_BarsV2.png", width = 10, height = 12, units = "in", res = 300)

p <- ggplot(filtered_data, 
            aes(x = Start, xend = End, y = SampleIndex, yend = SampleIndex)) +
  geom_segment(color = "red", size = 1, na.rm = TRUE) +  # Use bars instead of dots
  xlab("Chromosome Position") +
  ylab("Samples") +
  ggtitle("Gene Expression Presence") +
  theme_minimal() +
  theme(text = element_text(size = 20)) +
  scale_x_continuous(
    limits = c(0, 3e7),  # Set x-axis limits
    breaks = seq(0, 3e7, length.out = 4)  # Add tick marks within the limits
  ) 

print(p)
dev.off()
