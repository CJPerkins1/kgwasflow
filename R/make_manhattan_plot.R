#Making pretty manhattan plots from kgwasflow outputs

#I'm going to have three outputs from kgwasflow for ecoli, locule #, and tube burst. 
#I want to make manhattan plots that show both nonsignificant and significant hits colored by chromosome

#Making plot using the aligned kmer output

library(dplyr)
library(Rsamtools)
library(stringr)
library(ggplot2)

setwd("C:/Users/cperk/Desktop/palanivelu_lab/kgwasflow")

# Formatting data ----------------

#read in the bam file

kgwasflow_output <- scanBam(file = file.path(getwd(), "data/", "locule_number_kmers_alignment.bam"))

#formatting

rname_factor <- kgwasflow_output[[1]]$rname
rname_vector <- levels(rname_factor)[rname_factor]

df <- data.frame(
  kmer_name = kgwasflow_output[[1]]$qname,
  chr = rname_vector,
  chr_coord = kgwasflow_output[[1]]$pos,
  p_value = kgwasflow_output[[1]]$qname
)

df <- df %>%
  mutate(
    kmer_name = str_extract(kmer_name, str_extract(kmer_name, "^[^_]*")),
    chr = str_replace(chr, "SL4.0", ""),
    p_value = as.numeric(str_extract(p_value, "(?<=_)[0-9\\.e\\-]+")),
    log10_p_value = -log10(p_value)
  )

# converting the coords to genome coords, getting a list of chromosome boundaries for the plot
chrom_lengths <- data.frame(
  chr = sprintf("ch%02d", 0:12),
  length = c(9643250, 90863682, 53473368, 65298490, 64459972, 65269487, 
             47258699, 67883646, 63995357, 68513564, 64792705, 54379777, 66688036)
  )

chromosome_info <- chrom_lengths %>%
  arrange(chr) %>%
  filter(chr !="ch00") %>%
  mutate(end = cumsum(length),
         start = lag(end, default = 0) + 1) %>%
  select(chr, length, start, end)

chromosome_info <- chromosome_info  %>%
  mutate(midpoint = (start + end) / 2)

df <- df %>%
  left_join(chromosome_info, by = c("chr")) %>%
  mutate(genomic_coord = start + chr_coord)

locule_df <- locule_df %>%
  left_join(chromosome_info, by = c("chr")) %>%
  mutate(genomic_coord = start + chr_coord)

# Making the plot ---------------

# Identify unique chromosomes in your data
unique_chromosomes <- unique(df$chr)

# Create a custom color palette with alternating cyan and magenta
custom_palette <- rep(c("cyan3", "darkmagenta"), length.out = length(unique_chromosomes))

# Map the custom color palette to the chromosomes
chromosome_color_mapping <- setNames(custom_palette, unique_chromosomes)

df %>%
  ggplot(aes(x = genomic_coord, y = log10_p_value, color = factor(chr))) +
  geom_point(size = 2) +
  geom_vline(data = chromosome_info, aes(xintercept = start), linewidth = 0.2, color = "gray", linetype = "dashed") +
  geom_vline(data = chromosome_info, aes(xintercept = end), linewidth = 0.2, color = "gray", linetype = "dashed") +
  scale_x_continuous(breaks = chromosome_info$midpoint,
                     labels = chromosome_info$chr,
                     limits = c(1, 772876783),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 18, 2),
                     labels = seq(0, 18, 2),
                     limits = c(0, 18)) +
  labs(title = "Significant genomic loci for locule number", x = "Chromosome", y = "-log10 p-value") +
  scale_color_manual(values = chromosome_color_mapping) +  # Use the custom color palette
  theme_bw() +
  theme(axis.title = element_text(size = 26, face = 'bold'),
        axis.text = element_text(size = 18, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 18, face = 'bold', color = 'black'),
        plot.title = element_text(size = 28, face = 'bold', margin = margin(0, 0, 10, 0)),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 1, color = 'black'),
        axis.ticks = element_line(linewidth = 1, color = 'black'), axis.ticks.length = unit(8, 'pt'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, 'cm'),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 18, face = 'bold', color = 'black'),
        legend.text = element_text(size = 14, face = 'bold', color = 'black'),
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside") +
  guides(color = FALSE)  # Remove the color legend



ggsave(filename = file.path(getwd(), "R_plots", "locule_number_gwas.png"),
       device = 'png',
       width = 14,
       height = 8,
       dpi = 400,
       units = 'in')
