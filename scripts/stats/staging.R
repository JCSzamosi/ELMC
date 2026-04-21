


#### MAGs features & counts

##### Count within donors to get DB sizes

```{r}
## mags
all_mags = read.delim('data/clean/MAGids.txt', header = FALSE)
colnames(all_mags) = 'BinLab'
all_mags = (all_mags
            %>% mutate(Donor = str_replace(substr(BinLab, 1, 4),
                                           'Don', 'Donor'))
            %>% left_join(bin_tax, by = 'BinLab'))

### get counts per family & filter down to those present in both donors
all_fmag_cts = (all_mags
               %>% count(Donor, Family, name = 'size'))

mag_fam_both = (all_fmag_cts
                %>% select(Donor, Family)
                %>% unique()
                %>% count(Family)
                %>% filter(n > 1)
                %>% pull(Family))

### get counts per phylum & filter down to those present in both donors
all_pmag_cts = (all_mags
                %>% count(Donor, Phylum, name = 'size'))

mag_phy_both = (all_pmag_cts
                %>% select(Donor, Phylum)
                %>% unique()
                %>% count(Phylum)
                %>% filter(n > 1)
                %>% pull(Phylum))
```

```{r}
### get all mag coverage info
mag_cov = read.delim('data/clean/MAGCover.txt')

### join it with the taxon information
mag_cov = (mag_cov
           %>% mutate(Donor = str_replace(substr(feature, 1, 4),
                                          'Don', 'Donor'))
           %>% left_join(bin_fam_phy, by = c('feature' = 'BinLab')))

### Count how many mags in each family/phylum in each sample meet the cutoff,
### and fill in zeroes
fmag_count = (mag_cov
             %>% filter(coverage >= 75)
             %>% count(Sample, Donor, Family, name = 'count')
             %>% pivot_wider(names_from = 'Family', values_from = 'count',
                             values_fill = 0)
             %>% pivot_longer(-(Sample:Donor),
                              names_to = 'Family', values_to = 'count')
             %>% left_join(all_fmag_cts, by = c('Donor', 'Family')))

pmag_count = (mag_cov
              %>% filter(coverage >= 75)
              %>% count(Sample, Donor, Phylum, name = 'count')
              %>% pivot_wider(names_from = 'Phylum', values_from = 'count',
                              values_fill = 0)
              %>% pivot_longer(-(Sample:Donor),
                               names_to = 'Phylum', values_to = 'count')
              %>% left_join(all_pmag_cts, by = c('Donor', 'Phylum')))

### merge with in the sample data

fmag_count = (fmag_count
             %>% left_join(dat_map, by = 'Sample')
             %>% filter(!is.na(Group))
             %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                          (Group == 'AB' & Donor == 'DonorA') |
                                            (Group == 'BA' & Donor == 'DonorB') ~
                                            'first',
                                          TRUE ~ 'second'),
                        order = factor(order,
                                       levels = c('mixed', 'second', 'first')))
              %>% rename(feature = Family))

pmag_count = (pmag_count
              %>% left_join(dat_map, by = 'Sample')
              %>% filter(!is.na(Group))
              %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                           (Group == 'AB' & Donor == 'DonorA') |
                                             (Group == 'BA' & Donor == 'DonorB') ~
                                             'first',
                                           TRUE ~ 'second'),
                        order = factor(order,
                                       levels = c('mixed', 'second', 'first')))
              %>% rename(feature = Phylum))
```

#### Bins features & counts

```{r}
### get all bin coverage info
bin_cov = read.delim('data/clean/binCover.txt')

### join it with the taxon information
bin_cov = (bin_cov
           %>% mutate(Donor = str_replace(substr(feature, 1, 4),
                                          'Don', 'Donor'))
           %>% left_join(bin_tax, by = c('feature' = 'BinLab')))
bin_cov = (bin_cov %>% filter(!grepl('UnBin', feature)))

# Get the families and phyla that are present in both donors

bin_fam_both = (bin_cov
                %>% select(Donor, Family)
                %>% unique()
                %>% na.omit()
                %>% count(Family)
                %>% filter(n > 1)
                %>% pull(Family))
bin_phy_both = (bin_cov
                %>% select(Donor, Phylum)
                %>% unique()
                %>% na.omit()
                %>% count(Phylum)
                %>% filter(n > 1)
                %>% pull(Phylum))

# Get the number of bins in each donor at the family and phylum levels
all_fbin_cts = (bin_cov
                %>% filter(Family %in% bin_fam_both)
                %>% count(Donor, Family, name = 'size'))
all_pbin_cts = (bin_cov
                %>% filter(Phylum %in% bin_fam_both)
                %>% count(Donor, Phylum, name = 'size'))

### Count how many bins in each family/phylum in each sample meet the cutoff,
### and fill in zeroes
fbin_count = (bin_cov
             %>% filter(coverage >= 75)
             %>% count(Sample, Donor, Family, name = 'count')
             %>% pivot_wider(names_from = 'Family', values_from = 'count',
                             values_fill = 0)
             %>% pivot_longer(-(Sample:Donor),
                              names_to = 'Family', values_to = 'count')
             %>% left_join(all_fbin_cts, by = c('Donor', 'Family')))

pbin_count = (bin_cov
              %>% filter(coverage >= 75)
              %>% count(Sample, Donor, Phylum, name = 'count')
              %>% pivot_wider(names_from = 'Phylum', values_from = 'count',
                              values_fill = 0)
              %>% pivot_longer(-(Sample:Donor),
                               names_to = 'Phylum', values_to = 'count')
              %>% left_join(all_pbin_cts, by = c('Donor', 'Phylum')))

### merge with in the sample data

fbin_count = (fbin_count
             %>% left_join(dat_map, by = 'Sample')
             %>% filter(!is.na(Group))
             %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                          (Group == 'AB' & Donor == 'DonorA') |
                                            (Group == 'BA' & Donor == 'DonorB') ~
                                            'first',
                                          TRUE ~ 'second'),
                        order = factor(order,
                                       levels = c('mixed', 'second', 'first')))
              %>% rename(feature = Family))

pbin_count = (pbin_count
              %>% left_join(dat_map, by = 'Sample')
              %>% filter(!is.na(Group))
              %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                           (Group == 'AB' & Donor == 'DonorA') |
                                             (Group == 'BA' & Donor == 'DonorB') ~
                                             'first',
                                           TRUE ~ 'second'),
                        order = factor(order,
                                       levels = c('mixed', 'second', 'first')))
              %>% rename(feature = Phylum))
```

#### Species Data

```{r}
## species
spc_dat = read.delim('data/clean/speciesAbund.txt')
spc_dat = (spc_dat
           %>% mutate(feature = gsub('[a-z]__','',feature))
           %>% separate_wider_delim(feature, '|', 
                                    names = c('Kingdom', 'Phylum','Class', 
                                              'Order', 'Family', 'Genus', 
                                              'Species', 'Strain'), 
                                    too_few = 'align_start')
           %>% mutate(Sample = case_when(Sample == 'Wild116' ~ 'DonB',
                                         Sample == 'X608' ~ 'DonorA',
                                         TRUE ~ Sample)))

# Make sure families nested within phyla

chk_spc = (spc_dat
           %>% select(Family, Phylum)
           %>% unique()
           %>% count(Family)
           %>% filter(n > 1))
chk_spc

# There are two families that are not nested within phyla. For simplicity, we'll just remove them

spc_dat = (spc_dat
           %>% filter(!(Family %in% chk_spc$Family)))

# Separate donor samples from recipient samples
spc_don = (spc_dat
           %>% filter(startsWith(Sample, 'Don'),
                      coverage > 0.1))
spc_rec = (spc_dat
           %>% filter(!startsWith(Sample, 'Don'),
                      coverage > 0.1))

# Identify strains that are present in both donors

str_both = (spc_don
            %>% select(Sample, Strain)
            %>% unique()
            %>% count(Strain)
            %>% filter(n > 1)
            %>% pull(Strain))

# Remove those strains from the data

spc_rec = (spc_rec
           %>% filter(!(Strain %in% str_both)))
spc_don = (spc_don
           %>% filter(!(Strain %in% str_both)))

# Make it easy to associate a strain with a source donor. Remove all strains
# that don't have a source donor

spc_rec = (spc_don
           %>% mutate(within_don = str_replace(Sample, 'B', 'orB'))
           %>% select(within_don, Strain)
           %>% unique()
           %>% right_join(spc_rec, by = 'Strain')
           %>% select(-within_don, -Strain, -coverage, Strain, within_don)
           %>% filter(!is.na(within_don)))

# get db size for the two donors

spc_fam_size = (spc_don
            %>% mutate(within_don = str_replace(Sample, 'B', 'orB'))
            %>% count(within_don, Family, name = 'size'))

spc_phy_size = (spc_don
            %>% mutate(within_don = str_replace(Sample, 'B', 'orB'))
            %>% count(within_don, Phylum, name = 'size'))


# Figure out which families are present in both donors:
spc_fam_both = (spc_don
            %>% select(Sample, Family)
            %>% unique()
            %>% count(Family)
            %>% filter(n > 1)
            %>% pull(Family))
head(spc_fam_both)
length(spc_fam_both)


# Figure out which phyla are present in both donors:

spc_phy_both = (spc_don
                %>% select(Sample, Phylum)
                %>% unique()
                %>% count(Phylum)
                %>% filter(n > 1)
                %>% pull(Phylum))
head(spc_phy_both)
length(spc_phy_both)

# Get just the taxon info

spc_tax = (spc_dat
           %>% filter(Family %in% spc_fam_both)
           %>% select(Phylum, Family)
           %>% unique())
all(spc_phy_both %in% spc_tax$Phylum)
spc_tax

# Filter down to just the phyla with >1 family

spc_tax = (spc_tax
           %>% count(Phylum)
           %>% filter(n > 1)
           %>% left_join(spc_tax))
spc_tax

# Get counts for families

fspc_count = (spc_rec
              %>% filter(Family %in% spc_tax$Family)
              %>% select(within_don, Family, Sample)
              %>% count(Sample, Family, within_don, name = 'count')
              %>% left_join(dat_map, by = 'Sample')
              %>% pivot_wider(names_from = Family, values_from = count,
                              values_fill = 0)
              %>% pivot_longer(-(Sample:Cage), names_to = 'Family',
                               values_to = 'count')
              %>% left_join(spc_fam_size, by = c('within_don', 'Family')))
fspc_count = (fspc_count
             %>% filter(!is.na(Group))
             %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                          (Group == 'AB' & within_don == 'DonorA') |
                                            (Group == 'BA' & within_don == 'DonorB') ~
                                            'first',
                                          TRUE ~ 'second'),
                        order = factor(order,
                                       levels = c('mixed', 'second', 'first')))
              %>% rename(feature = Family))

# Get counts for phyla

pspc_count = (spc_rec
              %>% filter(Phylum %in% spc_tax$Phylum)
              %>% select(within_don, Phylum, Sample)
              %>% count(Sample, Phylum, within_don, name = 'count')
              %>% left_join(dat_map, by = 'Sample')
              %>% pivot_wider(names_from = Phylum, values_from = count,
                              values_fill = 0)
              %>% pivot_longer(-(Sample:Cage), names_to = 'Phylum',
                               values_to = 'count')
              %>% left_join(spc_phy_size, by = c('within_don', 'Phylum')))
pspc_count = (pspc_count
             %>% filter(!is.na(Group))
             %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                          (Group == 'AB' & within_don == 'DonorA') |
                                            (Group == 'BA' & within_don == 'DonorB') ~
                                            'first',
                                          TRUE ~ 'second'),
                        order = factor(order,
                                       levels = c('mixed', 'second', 'first')))
              %>% rename(feature = Phylum))
```

## Models

### MAGs

#### Family

```{r}
fmag_mod_out = loop_mods(fmag_count, feats = mag_fam_both)
fmag_mod_lst = fmag_mod_out[['mods']]
fmag_pr_df = fmag_mod_out[['prdf']]

fmag_contr_df = get_contr_pv(fmag_mod_lst)
```

##### Diagnositcs

```{r}
fmag_fr_plt = ggplot(fmag_pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~feature, scales = 'free')
fmag_fr_plt

fmag_abs_fr_plt = ggplot(fmag_pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~feature, scales = 'free')
fmag_abs_fr_plt

fmag_qq_plt = ggplot(fmag_pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~feature, scales = 'free')
fmag_qq_plt
```

These models are disasters. There just isn't enough data for this to be
reliable. I don't feel great about these estimates.

#### Phylum

```{r}
pmag_mod_out = loop_mods(pmag_count, feats = mag_phy_both)
pmag_mod_lst = pmag_mod_out[['mods']]
pmag_pr_df = pmag_mod_out[['prdf']]
pmag_contr_df = get_contr_pv(pmag_mod_lst)
```

##### Diagnostics

```{r}
pmag_fr_plt = ggplot(pmag_pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~feature, scales = 'free')
pmag_fr_plt

pmag_abs_fr_plt = ggplot(pmag_pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~feature, scales = 'free')
pmag_abs_fr_plt

pmag_qq_plt = ggplot(pmag_pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~feature, scales = 'free')
pmag_qq_plt
```

These are also disasters and I hate them.

#### Merge

```{r}
fmag_contr_df = (fmag_contr_df
                 %>% mutate(taxlev = 'Family')
                 %>% rename(Family = feature)
                 %>% left_join(fam_phy, by = 'Family')
                 %>% select(-Family, Family))
pmag_contr_df = (pmag_contr_df
                 %>% mutate(taxlev = 'Phylum',
                            Family = 'Full Phylum')
                 %>% rename(Phylum = feature)
                 %>% select(-Family, -Phylum, Phylum, Family))

all_contr_df = rbind(fmag_contr_df, pmag_contr_df)
all_contr_df = (all_contr_df
                %>% mutate(Family = relevel(factor(Family), 'Full Phylum')))
```

### Bins

#### Family

```{r}
fbin_mod_out = loop_mods(fbin_count, feats = bin_fam_both)
fmag_mod_lst = fmag_mod_out[['mods']]
fmag_pr_df = fmag_mod_out[['prdf']]

fmag_contr_df = get_contr_pv(fmag_mod_lst)
```

##### Diagnositcs

```{r}
fmag_fr_plt = ggplot(fmag_pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~feature, scales = 'free')
fmag_fr_plt

fmag_abs_fr_plt = ggplot(fmag_pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~feature, scales = 'free')
fmag_abs_fr_plt

fmag_qq_plt = ggplot(fmag_pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~feature, scales = 'free')
fmag_qq_plt
```

These models are disasters. There just isn't enough data for this to be
reliable. I don't feel great about these estimates.

#### Phylum

```{r}
pmag_mod_out = loop_mods(pmag_count, feats = mag_phy_both)
pmag_mod_lst = pmag_mod_out[['mods']]
pmag_pr_df = pmag_mod_out[['prdf']]
pmag_contr_df = get_contr_pv(pmag_mod_lst)
```

##### Diagnostics

```{r}
pmag_fr_plt = ggplot(pmag_pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~feature, scales = 'free')
pmag_fr_plt

pmag_abs_fr_plt = ggplot(pmag_pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~feature, scales = 'free')
pmag_abs_fr_plt

pmag_qq_plt = ggplot(pmag_pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~feature, scales = 'free')
pmag_qq_plt
```

These are also disasters and I hate them.

#### Merge

```{r}
fmag_contr_df = (fmag_contr_df
                 %>% mutate(taxlev = 'Family')
                 %>% rename(Family = feature)
                 %>% left_join(fam_phy, by = 'Family')
                 %>% select(-Family, Family))
pmag_contr_df = (pmag_contr_df
                 %>% mutate(taxlev = 'Phylum',
                            Family = 'Full Phylum')
                 %>% rename(Phylum = feature)
                 %>% select(-Family, -Phylum, Phylum, Family))

all_contr_df = rbind(fmag_contr_df, pmag_contr_df)
all_contr_df = (all_contr_df
                %>% mutate(Family = relevel(factor(Family), 'Full Phylum')))
```
### Species

## Results

### MAGs

```{r}
ggplot(filt_contr_df, aes(Family, estimate)) +
  geom_point(aes(colour = within_don), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = estimate - (SE * 1.96),
                     ymax = estimate + (SE * 1.96),
                     colour = within_don),
                 position = position_dodge(width = 0.5)) +
  facet_grid(numerator~Phylum, scales = 'free_x') +
  scale_color_manual(values = c("DonorA" = "#7b3294",
                                  "DonorB" = "#008837"),
                     name = 'Donor') +
  scale_x_discrete(limits = rev) +
  geom_hline(yintercept = 0, linetype = 2)+
  coord_flip()

# head(filt_contr_df)
fcd_pl_df = (filt_contr_df
             %>% select(taxlev, numerator, estimate, within_don))

# head(fcd_pl_df)

ggplot(fcd_pl_df, aes(estimate, colour = taxlev)) +
  geom_density() +
  facet_grid(numerator ~ within_don)

```

### Bins

### Species



# Read in the Species information

# Remove those from the remaining data
spc_present = (spc_dat
          %>% filter(!(startsWith(Sample, 'Don')),
                     !(Strain %in% spc_both),
                     coverage > 0.1)
          %>% left_join(spc_don_ref, by = 'Strain'))
dim(spc_present)
head(spc_present)

spc_counts = (spc_present
              %>% count(Sample, Donor, Kingdom, Phylum, Class, Order, Family,
                        name = 'nSpecies'))
head(spc_counts)

spc_don_cts = (spc_don
              %>% count(Donor, Kingdom, Phylum, Class, Order, Family,
                        name = 'Total'))
head(spc_don_cts)

spc_counts = (spc_counts
              %>% left_join(spc_don_cts, by = c('Donor', 'Kingdom', 'Phylum',
                                                'Class', 'Order', 'Family')))
head(spc_counts)

