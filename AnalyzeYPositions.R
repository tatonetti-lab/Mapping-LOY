library(tidyr)
library(dplyr)
library(ggplot2)

# Read in the data
file1 <- read.table("Y_WGS_positions.txt", header = FALSE, sep = "\t")
names(file1)[1] <- "SampleID"

results<-file1

# Order by UniquePositions and assign SampleIndex
results_filtered<-results %>%
  filter(V2<3e7) %>%
  left_join(results_summary, by = "SampleID")

png("ChromosomeDensityHist.png",width=10,height=12,units="in",res=200)

r<-ggplot(results_filtered, aes(x = V2)) +
  geom_histogram(binwidth = 1e3, fill = "blue", alpha = 0.7) +
  scale_y_continuous(limits=c(0,80000),n.breaks=40)+
  xlab("Position on Chromosome") +
  ylab("Read Density") +
  ggtitle("Density of Reads Across Chromosome") +
  theme_minimal()+
  theme(text = element_text(size = 20))

print(r)
dev.off()

hist_data <- ggplot_build(r)$data[[1]]  # Extract histogram bin data
hist_data <- hist_data %>%
  select(x = xmin, xmax, count = y) %>%
  mutate(bin_mid = (x + xmax) / 2)  # Add bin midpoints for matching
hist_data <- hist_data %>%
  filter(is.finite(count)) # %>%
threshold <-2000
low_density_bins <- hist_data %>%
  filter(count < threshold)
annotation <- read.table("gene_chrY_positions.txt", header = TRUE, sep = "\t")  # Adjust file name and separator
sample_to_bins <- results_filtered %>%
  mutate(bin = cut(V2, breaks = c(hist_data$x, max(hist_data$xmax)))) %>%  # Assign bins
  group_by(bin) %>%
  summarise(SampleIDs = list(unique(SampleID)), .groups = "drop")  # Collect SampleIDs for each bin

low_density_bins <- low_density_bins %>%
  mutate(bin = cut(x, breaks = c(hist_data$x, max(hist_data$xmax)))) %>%
  left_join(sample_to_bins, by = "bin")

low_density_genes <- annotation %>%
  rowwise() %>%
  mutate(
    ReadDensity = sum(low_density_bins$count[
      (Start >= low_density_bins$x & Start <= low_density_bins$xmax) |
      (End >= low_density_bins$x & End <= low_density_bins$xmax)
    ], na.rm = TRUE),
    UniqueSamples = list(unique(unlist(low_density_bins$SampleIDs[
      (Start >= low_density_bins$x & Start <= low_density_bins$xmax) |
      (End >= low_density_bins$x & End <= low_density_bins$xmax)
    ])))
  ) %>%
  ungroup() %>%
  mutate(
    NumSamples = lengths(UniqueSamples)  # Count the number of unique samples
  ) %>%
  filter(ReadDensity>0) %>%
  arrange(ReadDensity)  # Order by read density

low_density_genes <- low_density_genes %>%
  mutate(UniqueSamples = sapply(UniqueSamples, function(x) paste(x, collapse = ",")))

# Save the results
write.table(
  low_density_genes,
  "low_density_2000genes.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


png("ChromosomeCDF.png",width=15,height=12,units="in",res=200)

t<-ggplot(hist_data, aes(x = count)) +
  stat_ecdf(geom = "point") +
  xlab("Read Density") +
  ylab("Cumulative Proportion") +
  ggtitle("Cumulative Distribution of Read Densities") #+

print(t)
dev.off()
bin_size <- 100000  # Adjust bin size as needed

results_filtered_s <- results_filtered %>%
  mutate(Bin = cut(V2, breaks = seq(0, 3e7, by = bin_size), labels = FALSE))

sparse_bins <- results_filtered_s %>%
  group_by(SampleID, Bin) %>%
  summarize(ReadCount = n(), .groups = "drop") %>%
  complete(SampleID, Bin, fill = list(ReadCount = 0))  # Fill missing bins with 0 reads

sparse_summary <- sparse_bins %>%
  group_by(SampleID) %>%
  summarize(SparseRegions = sum(ReadCount == 0), .groups = "drop")

results_ordered <- results_filtered_s %>%
  left_join(sparse_summary, by = "SampleID") %>%
  arrange(desc(SparseRegions)) %>%
  mutate(SampleIndex = as.integer(factor(SampleID, levels = unique(SampleID))))

write.table(results_ordered$SampleID,"orderedSamples.txt",sep="\t",row.names=FALSE)

png("ChromosomePosition.png", width = 10, height = 12, units = "in", res = 200)
p<-ggplot(results_ordered, aes(x = V2, y = SampleIndex)) +
  geom_point(na.rm = TRUE, size = 0.5) +
  xlab("Position on Chromosome") +
  ylab("Samples Ordered by Sparse Regions") +
  ggtitle("Samples Ordered by Sparse Regions") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),text = element_text(size = 10)) #+
print(p)
dev.off()

results_filtered <- results_filtered %>%
  mutate(Bin = cut(V2, breaks = seq(0, 3e7, by = 1e6), labels = FALSE))

bin_summary <- results_filtered %>%
  group_by(SampleID, Bin) %>%
  summarize(ReadsPerBin = n(), .groups = "drop")

png("ChromosomeBinSummary.png",width=10,height=12,units="in",res=200)

s<-ggplot(bin_summary, aes(x = Bin, y = ReadsPerBin, fill = SampleID)) +
  geom_col(position = "dodge") +
  xlab("Chromosome Bin") +
  ylab("Reads Per Bin") +
  ggtitle("Distribution of Reads Across Chromosome Bins") +
  theme_minimal() +
  theme(legend.position="none")

print(s)
dev.off()
