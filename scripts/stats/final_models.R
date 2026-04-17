library(tidyverse)
library(glmmTMB)
library(emmeans)
library(AfterSl1p)
source('scripts/stats/functions.R')
theme_set(theme_bw())


dat = read.delim('data/clean/coloniz_model.txt')
head(dat)
summary(dat)

dat = (dat
       %>% mutate(Group = factor(Group, levels = c('AB', 'BA', 'ABAB')),
                  Donor = factor(Donor),
                  order = factor(order, levels = c('mixed', 'second', 'first')),
                  feature = factor(feature),
                  count = raw_count,
                  size = feature_size)
       )
summary(dat)

ggplot(dat, aes(count)) +
  geom_histogram() +
  facet_wrap(~feature, scales = 'free')

pr_df = data.frame(resid = NA,
                   fitted = NA, 
                   feature = NA)
mod_lst = list()
for (feat in unique(dat$feature)){
  mod = glmmTMB(raw_count/feature_size ~ order * Donor +
                  (1 | Cage) +
                  (1 | Sample),
                  # offset = feature_size,
                family = 'gaussian',
                data = filter(dat, feature == feat))
  pr_df = rbind(pr_df,
                data.frame(resid = resid(mod),
                           fitted = fitted(mod),
                           feature = feat))
  mod_lst[[feat]] = mod
}

pr_df = pr_df[-1,]
fr_plt = ggplot(pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~feature, scales = 'free')
fr_plt

abs_fr_plt = ggplot(pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~feature, scales = 'free')
abs_fr_plt

qq_plt = ggplot(pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~feature, scales = 'free')
qq_plt

# Extract p-values & coefficients

contr_df = data.frame(contrast = NA,
                      estimate = NA,
                      SE = NA, 
                      df = NA, 
                      t.ratio = NA,
                      p.value = NA,
                      within_don = NA,
                      feature = NA)
for (feat in names(mod_lst)){
  mod = mod_lst[[feat]]
  emA = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonorA'))
  contr_df = rbind(contr_df,
                   cbind(emA$contrasts, data.frame(within_don = 'DonorA',
                                                   feature = feat))
  )
  emB = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonorB'))
  contr_df = rbind(contr_df,
                   cbind(emB$contrasts, data.frame(within_don = 'DonorB',
                                                   feature = feat))
  )
}
contr_df = contr_df[-1,]
contr_df = separate_wider_delim(contr_df, contrast, ' - ', names = c('numerator',
                                                                     'denominator'))
contr_df


ggplot(contr_df, aes(numerator, estimate, colour = within_don)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = estimate - 1.96*SE,
                     ymax = estimate + 1.96*SE),
                 position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = 2) +
  coord_flip() +
  facet_wrap(~feature) +
  scale_color_manual(values = c("DonorA" = "#7b3294",
                                  "DonorB" = "#008837"),
                     name = 'Donor') +
  labs(y = 'Difference from mixed',
       x = 'Order')

# Get the p-values from this data frame:

contr_df

# For MAGs and Species, model within family

dat_map = (dat
           %>% select(Sample, Group, Category, Cage)
           %>% unique())

# Read in bin taxonomy information:

bin_tax = read.delim('data/clean/binQT.txt')

# Check if families are unique

chk_fams = (bin_tax
            %>% select(Phylum:Family)
            %>% unique())
dim(chk_fams)
n_distinct(chk_fams$Family)

# They are. I don't need to carry all the taxon info.

bin_tax = select(bin_tax, BinLab, Family, Phylum)
fam_phy = (bin_tax
           %>% select(Family, Phylum)
           %>% unique())

# Read in the full list of MAGs
all_mags = read.delim('data/clean/MAGids.txt', header = FALSE)
colnames(all_mags) = 'BinLab'

# Give each MAG its donor and taxon information
all_mags = (all_mags
            %>% mutate(Donor = substr(BinLab, 1, 4))
            %>% left_join(bin_tax, by = 'BinLab'))

# Get counts of the number of MAGs from each donor in each family
all_fmag_cts = (all_mags
               %>% count(Donor, Family, Phylum, name = 'Total'))

# Get counts of the number of MAGs from each donor in each phylum
all_pmag_cts = (all_mags
                %>% count(Donor, Phylum, name = 'Total'))

# Get families that are represented in both donors:
mag_fam_both = (all_fmag_cts
                %>% select(Donor, Family)
                %>% unique()
                %>% count(Family)
                %>% filter(n > 1)
                %>% pull(Family))

# Get phyla that are represented in both donors:
mag_phy_both = (all_pmag_cts
                %>% select(Donor, Phylum)
                %>% unique()
                %>% count(Phylum)
                %>% filter(n > 1)
                %>% pull(Phylum))

# Read in the MAG coverage information:
mag_cov = read.delim('data/clean/MAGCover.txt')

# Join it with the taxon information
mag_cov = (mag_cov
           %>% mutate(Donor = substr(feature, 1, 4))
           %>% left_join(bin_tax, by = c('feature' = 'BinLab')))

# Count how many mags in each family in each sample meet the cutoff, and fill in zeroes
fmag_count = (mag_cov
             %>% filter(coverage >= 75)
             %>% count(Sample, Donor, Family, name = 'nMAGs')
             %>% pivot_wider(names_from = 'Family', values_from = 'nMAGs',
                             values_fill = 0)
             %>% pivot_longer(-(Sample:Donor),
                              names_to = 'Family', values_to = 'nMAGs')
             %>% left_join(all_fmag_cts, by = c('Donor', 'Family')))

# Count how many mags in each phylum in each sample meet the cutoff, and fill in zeroes
pmag_count = (mag_cov
              %>% filter(coverage >= 75)
              %>% count(Sample, Donor, Phylum, name = 'nMAGs')
              %>% pivot_wider(names_from = 'Phylum', values_from = 'nMAGs',
                              values_fill = 0)
              %>% pivot_longer(-(Sample:Donor),
                               names_to = 'Phylum', values_to = 'nMAGs')
              %>% left_join(all_pmag_cts, by = c('Donor', 'Phylum')))

# Add in the sample data

fmag_count = (fmag_count
             %>%left_join(dat_map, by = 'Sample')
             %>% filter(!is.na(Group))
             %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                          (Group == 'AB' & Donor == 'DonA') |
                                            (Group == 'BA' & Donor == 'DonB') ~
                                            'first',
                                          TRUE ~ 'second')))

pmag_count = (pmag_count
              %>% left_join(dat_map, by = 'Sample')
              %>% filter(!is.na(Group))
              %>% mutate(order = case_when(Group == 'ABAB' ~ 'mixed',
                                           (Group == 'AB' & Donor == 'DonA') |
                                             (Group == 'BA' & Donor == 'DonB') ~
                                             'first',
                                           TRUE ~ 'second')))

# Model by family:

fmag_pr_df = data.frame(resid = NA,
                   fitted = NA, 
                   family = NA)
fmag_mod_lst = list()
for (fam in mag_fam_both){
  mod = glmmTMB(nMAGs/Total ~ order + Donor +
                  (1 | Cage) +
                  (1 | Sample),
                family = 'gaussian',
                data = filter(fmag_count, Family == fam))
  fmag_pr_df = rbind(fmag_pr_df,
                data.frame(resid = resid(mod),
                           fitted = fitted(mod),
                           family = fam))
  fmag_mod_lst[[fam]] = mod
}
fmag_pr_df = fmag_pr_df[-1,]
head(fmag_pr_df)

fmag_fr_plt = ggplot(fmag_pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~family, scales = 'free')
fmag_fr_plt

fmag_abs_fr_plt = ggplot(fmag_pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~family, scales = 'free')
fmag_abs_fr_plt

# These are disassters

fmag_qq_plt = ggplot(fmag_pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~family, scales = 'free')
fmag_qq_plt

# Extract p-values & coefficients

fmag_contr_df = data.frame(contrast = NA,
                      estimate = NA,
                      SE = NA, 
                      df = NA, 
                      t.ratio = NA,
                      p.value = NA,
                      within_don = NA,
                      taxlevel = NA,
                      Family = NA)
for (fam in names(fmag_mod_lst)){
  mod = fmag_mod_lst[[fam]]
  emA = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonA'))
  fmag_contr_df = rbind(fmag_contr_df,
                   cbind(emA$contrasts, data.frame(within_don = 'DonA',
                                                   taxlevel = 'Family',
                                                   Family = fam))
  )
  emB = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonB'))
  fmag_contr_df = rbind(fmag_contr_df,
                   cbind(emB$contrasts, data.frame(within_don = 'DonB',
                                                   taxlevel = 'Family',
                                                   Family = fam))
  )
}
fmag_contr_df = fmag_contr_df[-1,]
fmag_contr_df = (fmag_contr_df 
                %>% separate_wider_delim(contrast, ' - ', 
                                    names = c('numerator', 'denominator'))
                %>% left_join(fam_phy, by = 'Family'))
head(fmag_contr_df)

# Model by phylum:

pmag_pr_df = data.frame(resid = NA,
                   fitted = NA, 
                   phylum = NA)
pmag_mod_lst = list()
for (phy in mag_phy_both){
  mod = glmmTMB(nMAGs/Total ~ order + Donor +
                  (1 | Cage) +
                  (1 | Sample),
                family = 'gaussian',
                data = filter(pmag_count, Phylum == phy))
  pmag_pr_df = rbind(pmag_pr_df,
                data.frame(resid = resid(mod),
                           fitted = fitted(mod),
                           phylum = phy))
  pmag_mod_lst[[phy]] = mod
}
pmag_pr_df = pmag_pr_df[-1,]
head(pmag_pr_df)

pmag_fr_plt = ggplot(pmag_pr_df, aes(fitted, resid)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~phylum, scales = 'free')
pmag_fr_plt

pmag_abs_fr_plt = ggplot(pmag_pr_df, aes(fitted, abs(resid))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = 'lm') +
  facet_wrap(~phylum, scales = 'free')
pmag_abs_fr_plt

# These are disassters

pmag_qq_plt = ggplot(pmag_pr_df, aes(sample = resid)) +
  geom_qq() +
  geom_qq_line() +
  facet_wrap(~phylum, scales = 'free')
pmag_qq_plt

# Extract p-values & coefficients

pmag_contr_df = data.frame(contrast = NA,
                      estimate = NA,
                      SE = NA, 
                      df = NA, 
                      t.ratio = NA,
                      p.value = NA,
                      within_don = NA,
                      taxlevel = NA,
                      Phylum = NA)
for (phy in names(pmag_mod_lst)){
  mod = pmag_mod_lst[[phy]]
  emA = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonA'))
  pmag_contr_df = rbind(pmag_contr_df,
                   cbind(emA$contrasts, data.frame(within_don = 'DonA',
                                                   taxlevel = 'Phylum',
                                                   Phylum = phy))
  )
  emB = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonB'))
  pmag_contr_df = rbind(pmag_contr_df,
                   cbind(emB$contrasts, data.frame(within_don = 'DonB',
                                                   taxlevel = 'Phylum',
                                                   Phylum = phy))
  )
}
pmag_contr_df = pmag_contr_df[-1,]
pmag_contr_df = (pmag_contr_df 
                %>% separate_wider_delim(contrast, ' - ', 
                                    names = c('numerator', 'denominator')))
head(pmag_contr_df)
head(fmag_contr_df)

pmag_contr_df = (pmag_contr_df
                 %>% mutate(Family = 'Full Phylum')
                 %>% select(-Phylum, Phylum))

head(pmag_contr_df)

all_contr_df = rbind(fmag_contr_df, pmag_contr_df)
all_contr_df = (all_contr_df
                %>% mutate(Family = relevel(factor(Family), 'Full Phylum')))

ggplot(all_contr_df, aes(Family, estimate)) +
  geom_point(aes(colour = within_don), position = position_dodge(width = 0.5)) +
  facet_grid(numerator~Phylum, scales = 'free_x') +
  scale_color_manual(values = c("DonA" = "#7b3294",
                                  "DonB" = "#008837"),
                     name = 'Donor') +
  scale_x_discrete(limits = rev) +
  geom_hline(yintercept = 0, linetype = 2)+
  coord_flip()



# Read in the Species information
spc_dat = read.delim('data/clean/speciesAbund.txt')
spc_dat = (spc_dat
           %>% mutate(feature = gsub('[a-z]__','',feature))
           %>% separate_wider_delim(feature, '|', 
                                    names = c('Kingdom', 'Phylum','Class', 
                                              'Order', 'Family', 'Genus', 
                                              'Species', 'Strain'), 
                                    too_few = 'align_start')
           %>% mutate(Sample = case_when(Sample == 'Wild116' ~ 'DonB',
                                         Sample == 'X608' ~ 'DonA',
                                         TRUE ~ Sample)))

# Figure out which strains are present in both donors:
spc_don = (spc_dat
           %>% filter(Sample %in% c('DonA', 'DonB'),
                 coverage >= 0.1)
           %>% rename(Donor = Sample))
dim(spc_don)
head(spc_don)
spc_both = (spc_don
            %>% count(Strain)
            %>% filter(n > 1)
            %>% pull(Strain))
head(spc_both)
length(spc_both)

spc_don_ref = (spc_don
           %>% filter(!(Strain %in% spc_both))
           %>% select(Donor, Strain))

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

