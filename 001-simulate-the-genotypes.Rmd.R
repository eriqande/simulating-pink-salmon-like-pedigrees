library(tidyverse)
library(gdropR)

# read in the allele freqs
lines <- read_lines("data/prior.txt")   
snp_freqs <- lines[str_detect(lines, "aFreq")] %>%
  str_split(" +") %>%
  "[["(1) %>%
  "["(-1) %>%
  as.numeric()

L <- length(snp_freqs)

m1 <- tibble(
  Chrom = "Unk",
  Locus = paste0("PinkSal", 1:L),
  Pos = 1:L,
  LocIdx = 1:L,
  a1 = snp_freqs,
  a2 = 1 - a1
) %>%
  pivot_longer(
    cols = c(a1, a2),
    names_to = "Allele",
    values_to = "Freq"
  ) %>%
  group_by(Locus) %>%
  mutate(AlleIdx = 1:n()) %>%
  ungroup() %>%
  select(Chrom, Locus, Pos, Allele, LocIdx, AlleIdx, Freq)

markers <- reindex_markers(m1)
markers$Pos <- NA_real_


# place markers randomly
genome <- chinook_chromosomes %>%
  mutate(
    idx = 1:n(),
    chrom = name1,
    num_bases = length,
    scaled_length = num_bases / max(num_bases)
  ) %>%
  select(idx, chrom, scaled_length, num_bases)


mapped_markers <- sprinkle_markers_into_genome(
  markers,
  genome
) %>%
  group_by(Chrom) %>%
  mutate(Locus = paste0(Chrom, "-", rep(1:(n()/2), each = 2))) %>%
  ungroup()


# Now get the pedigree from the inferred ped file
last_ped <- read_table2(
  "data/ped.txt", 
  col_names = c("gen", "Kid", "Pa", "Ma")
  ) %>%
  filter(gen == 24) %>%
  group_by(Kid) %>%
  slice(n()) %>%
  ungroup() %>%
  mutate(across(.fns = as.character)) %>%
  select(-gen)


pedigree <- last_ped %>%
  add_explicit_founder_parents() %>%
  add_pedigree_sex_column(fem_prob = 1.0)


# now simulate:
set.seed(123)
wide_genos <- simulate_linked_genotypes(
  mapped_markers,
  pedigree
) %>%
  arrange(as.integer(indiv)) %>%
  filter(as.integer(indiv) <= 6604)




long_genos <- wide_genos %>%
  pivot_longer(
    cols = c(-indiv, -sex),
    names_to = "Locus",
    values_to = "Genotype"
  )

err_genos <- long_genos %>%
  mutate(
    obs_geno = biallelic_geno_error(
      genos = Genotype, 
      het_miscall = 0.02, 
      hom_miscall = 0.005
    )
  )

# now, put the generation level of each individual
gen_levels <- pedigree %>%
  mutate(
    generation = 
      case_when(
        as.numeric(Kid) < 550 ~ 0L,
        as.numeric(Kid) >= 550 & as.numeric(Kid) < 3895 ~ 1L,
        as.numeric(Kid) >= 3895 & as.numeric(Kid) < 6605 ~ 2L,
        TRUE ~ NA_integer_
      )
  )

  