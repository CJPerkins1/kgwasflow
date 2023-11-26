#  Samples - fastas

#  samples.tsv will include 5 columns, sample_name, library_name, fq1, fq2, sra
#  sample_name will include genome name
#  library_name will have the name of the library?
#  fq1 is the path to your first fastq file.gz
#  fq2 is the path to the second fastq file.gz
#  sra is optional option to include SRR number (I won't be doing this bc Cedar has already compiled all the fastq files)

# building a sheet that contains everything so that it's all aligned, and then I will pull the necessary columns from this sheet for the .tsv files I need for kgwasflow

library(dplyr)
library(stringr)
library(tidyr)
#---------------------

# Sample sheet for checking if kgwasflow works. Will use only 3 fastqs. 

# reading in accessions sheet containing the names for each accession, including CW
setwd('C:/Users/cperk/Desktop/palanivelu_lab/kgwasflow')
accessions <- read.csv('data/accessions.csv') %>%
  select(name_CW:name_original_razifard)

# reading in varitome metadata sheet containing all the seq info for the varitome div panel, including sras & library
varitome_metadata <- read.csv('data/solavar_metadata.txt')
var_sra <- varitome_metadata %>% # pulling out necessary columns, SRA#, sample name, SRX#
  select(Run, Sample.Name, Experiment) %>%
  mutate(Sample.Name = str_replace_all(Sample.Name, ' ', '')) #some of the samples have an unnecessary space

# alonge_metadata <- read.csv('data/alonge_metadata.txt')
# alonge_sra <- alonge_metadata %>%
  # select(Run, Sample.Name, Experiment)

#joining sras with cw

var_sheet <- var_sra %>%
  left_join(accessions, join_by(Sample.Name == name_original_razifard))

var_sheet <- var_sheet %>%
  arrange(desc(Run))

# There are multiple libraries for some of the accessions, figure out how to deal 
# with these. I'll choose a random one for the test.
qc <- var_sheet %>%
  group_by(name_CW) %>%
  summarize(n = n())

set.seed(13)
var_sheet <- var_sheet %>%
  group_by(name_CW) %>%
  slice_sample(n = 1)

#Keeping only necessary columns
var_sheet <- var_sheet %>%
  select(Run, Sample.Name, Experiment, name_CW) %>%
  drop_na()

#var_sheet has all necessary info for making paths and samples.tsv sheet. Now I need to add the pheno data so it all aligns nicely

#burst pheno:
diversity_panel <- read.csv('data/pollen_diversity_panel_phenotypes.csv')

#save columns of interest in new dataframe
dp_burst <- diversity_panel[, c('accession', 'adjusted_integral_increase')]
dp_burst <- dp_burst %>%
  rename(name_CW = accession)

#join with var_sheet:
var_sheet <- var_sheet %>%
  left_join(dp_burst, join_by(name_CW == name_CW))

#pistil-anther ratio
flower_phenotypes <- read.csv("data/flower_phenotypes.csv") %>%
  group_by(accession_id) %>%
  summarize(mean_anther_pistil_ratio = mean(anther_over_pistil)) %>%
  rename(name_CW = accession_id)

var_sheet <- var_sheet %>%
  left_join(flower_phenotypes, join_by(name_CW == name_CW))

#locule number
# Using the phenotyping data from the Varitome project, located here:
# https://solgenomics.net/ftp/varitome/GWAS/normalized_phen/
locule_num <- read.table(
  file = file.path(getwd(), "data", "varitome_locule_number.fam"),
  sep = '\t',
  header = FALSE
) %>%
  select(c(1, 6)) %>%
  rename(name_original_razifard = V1, locule = V6)

#add this data to the accessions sheet
locule_num <- locule_num %>%
  left_join(accessions, by = c("name_original_razifard")) %>%
  select(accession_id = name_CW, locule) %>%
  drop_na()

locule_num <- locule_num %>%
  rename(name_CW = accession_id)

var_sheet <- var_sheet %>%
  left_join(locule_num, join_by(name_CW))

#Drop all Nas from the var_sheet
var_sheet <- var_sheet %>%
  drop_na()

#var_sheet now has everything necessary to build paths for fastqs, build samples.tsv, build pheno.tsv for multiple different phenos

#---------------------------
# building paths for fastq inputs and making samples.tsv

#Path to cedar's fastq files
base_path <- "/xdisk/rpalaniv/cedar/kmers-gwas/genomes/varitome/fastq"

# create tsv dataframe
samples <- var_sheet%>%
  select(Run, Sample.Name, name_CW, Experiment)
#remove nas
samples <- na.omit(samples)

samples$fq1 <- file.path(base_path, paste0(samples$Experiment, '_', samples$Run, '_1.fastq.gz'))

samples$fq2 <- file.path(base_path, paste0(samples$Experiment, '_', samples$Run, '_2.fastq.gz'))

samples$sra <- samples$Run
samples$sample_name <- samples$name_CW
samples$library_name <- samples$Sample.Name

samples <- samples[, c('sample_name', 'library_name', 'fq1', 'fq2', 'sra')]

#Save as tsv in the config folder in kgwasflow dir

samples %>%
  write.table(
    file = file.path(getwd(), "config", "samples.tsv"),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
  )

#-----------------------------
#  Phenos
#  kgwasflow requires a tsv input with all your phenotypes in the format: pheno_name.pheno
#  The pheno files have one column with accession_id, and another column with phenotype_value

#burst integral increase pheno
burst_pheno <- var_sheet[, c("name_CW", 'adjusted_integral_increase')]

burst_pheno %>%
  write.table(
    file = file.path(getwd(), 'config', 'burst_integral.pheno'),
    sep = '\t',
    quote = F,
    row.names = F
  )

# pistil anther ratio pheno
pa_ratio <- var_sheet[, c("name_CW", 'mean_anther_pistil_ratio')]

pa_ratio %>%
  write.table(
    file = file.path(getwd(), 'config', 'pa_ratio.pheno'),
    sep = '\t',
    quote = F,
    row.names = F
  )

#locule num pheno
loc_num <- var_sheet[, c('name_CW', 'locule')]

loc_num %>%
  write.table(
    file = file.path(getwd(), 'config', 'locule_num.pheno'),
    sep = '\t',
    quote = F,
    row.names = F
  )

#----------------------------

# Make a phenos tsv file that contains all the phenotype names and paths to phenotypes, starting with just locule number, here is the info for the other phenos when necessary:
# 'burst_integral', '/xdisk/rpalaniv/cjperkins1/kgwasflow/phenos/burst_integral.pheno'
# 'pa_ratio', '/xdisk/rpalaniv/cjperkins1/kgwasflow/phenos/pa_ratio.pheno'

phenos <- data.frame(
  pheno_name = c('locule_num'),
  pheno_path = c('/xdisk/rpalaniv/cjperkins1/kgwasflow/phenos/locule_num.pheno')
)

phenos %>%
  write.table(
    file = file.path(getwd(), 'config', 'phenos.tsv'),
    sep = '\t',
    quote = F,
    row.names = F
  )

#pray this works....

#-----------------------------
#code for making samples.tsv for ecoli test dataset
library(dplyr)
library(stringr)

# building paths for fastq inputs

setwd('C:/Users/cperk/Desktop/palanivelu_lab/kgwasflow/R/ecoli_test')

ecoli_samples <- read.delim('samples.tsv')
samplesheet <- read.csv('samplesheet.csv')

ecoli_samples$fq1 <- samplesheet$fastq_1
ecoli_samples$fq2 <- samplesheet$fastq_2

write.table(ecoli_samples, "samples.tsv", sep = "\t", row.names = F)

#Path to cedar's fastq files
#base_path <- "/xdisk/rpalaniv/cedar/kmers-gwas/genomes/varitome/fastq"

# create tsv dataframe
#samples <- var_sras_aligned %>%
  #select(Run, Sample.Name, name_CW)
#remove nas
#samples <- na.omit(samples)

samples$fq1 <- file.path(base_path, paste0(srx_ids, '_', samples$Run, '_1.fastq.gz'))

samples$fq2 <- file.path(base_path, paste0(srx_ids, '_', samples$Run, '_2.fastq.gz'))

samples$sra <- samples$Run
samples$sample_name <- samples$Sample.Name
samples$library_name <- samples$name_CW

samples <- samples[, c('sample_name', 'library_name', 'fq1', 'fq2', 'sra')]