get_contr_pv = function(mod_lst){
  # Get contrasts and p-values from models produced by loop_mods()
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
    emB = emmeans(mod, specs = trt.vs.ctrl1 ~ order, at = list(Donor = 'DonorB'))
    df = rbind(cbind(emA$contrasts, 
                     data.frame(within_don = 'DonorA', feature = feat)),
               cbind(emB$contrasts, 
                     data.frame(within_don = 'DonorB', feature = feat)))
    df = df[,colnames(contr_df)]
    contr_df = rbind(contr_df, df)
  }
  contr_df = contr_df[-1,]
  contr_df = separate_wider_delim(contr_df, contrast, ' - ', names = c('numerator',
                                                                       'denominator'))
  return(contr_df)
}

use_loop_mods = function(dat, ft, tl){
  df = filter(dat, FeatureType == ft, TaxLev == tl)
  outlst = loop_mods(df)
  outlst[['prdf']]$FeatureType = ft
  outlst[['prdf']]$TaxLev = tl
  return(outlst)
}

loop_mods = function(dat){
  # Loop over all features/taxa to run the models and collect outputs
  pr_df = data.frame(resid = NA,
                     fitted = NA, 
                     feature = NA)
  mod_lst = list()
  for (feat in unique(dat$feature)){
    mod = run_mod(filter(dat, feature == feat))
    pr_df = rbind(pr_df,
                  data.frame(resid = resid(mod),
                             fitted = fitted(mod),
                             feature = feat))
    mod_lst[[feat]] = mod
  }
  pr_df = pr_df[-1,]
  return(list(mods = mod_lst, prdf = pr_df))
}


run_mod = function(dat, re = TRUE){
  # Run the model that we use everywhere
  form = 'count/size ~ order * Donor'
  if (re){
    form = paste(form, '+ (1 | Cage) + (1 | Sample)')
  }
  mod = glmmTMB(formula(form),
                family = 'gaussian',
                data = dat)
}

sum_bins = function(df){
  breaks = n_distinct(df$bins)
  labs = paste('bin', 1:breaks, sep = '')
  df = (df
        %>% group_by(bins)
        %>% summarize(mean = mean(raw_count),
                      sd = sd(raw_count),
                      n = length(raw_count))
        %>% ungroup()
        %>% mutate(bins = factor(bins),
                   bins = order_levs(bins, mean),
                   bins = factor(bins, 
                                 labels = paste('bin', 1:breaks, sep = ''))))
  return(df)
}

mk_brk_df = function(df, breaks){
  bin_df = (df
            %>% select(raw_count)
            %>% mutate(bins = cut(raw_count, breaks = breaks)))
  bin_df = sum_bins(bin_df)
  return(bin_df)
}

plt_mean_sd = function(df, loglog = 'n'){
  plt = ggplot(df, aes(mean, sd)) +
    geom_point(alpha = 0.4) +
    xlim(0, 5e4)
  if (loglog == 'y'){
    plt = plt + 
      scale_y_log10() +
      scale_x_log10(limits = c(1, 5e4))
  }
  return(plt)
}

plt_bs = function(df){
  plt = ggplot(df, aes(bins, n)) +  
    geom_point() +  rotate_ticks() +
    ylab('Samples per bin') +
    xlab('bins (sorted by mean)') +
    ylim(0, 200)
  
  return(plt)
}

mk_qq = function(df, grp){
  cols = paste(c('fit', 'res'), grp, sep = '_')
  plt = ggplot(df, aes(sample = .data[[cols[2]]])) +
    geom_qq() +
    geom_qq_line()
  return(plt)
}

mk_fr = function(df, grp, abs = FALSE, guide = 'legend'){
  cols = paste(c('fit', 'res'), grp, sep = '_')
  if (abs){
    plt = ggplot(df, aes(x = .data[[cols[1]]],
                         y = abs(.data[[cols[2]]])))
  } else {
    plt = ggplot(df, aes(x = .data[[cols[1]]],
                         y = .data[[cols[2]]]))
    
  }
  
  plt = plt + geom_point(aes(colour = feature)) +
    scale_colour_brewer(palette = 'Dark2', guide = guide) +
    geom_smooth(method = 'lm')
  return(plt)
}

mk_plt_row = function(df, grp){
  rw = plot_grid(mk_qq(fit_df, grp), mk_fr(fit_df, grp, guide = 'none'),
                 mk_fr(fit_df, grp,TRUE, guide = 'none'),
                 ncol = 3, 
                 labels = grp,
                 label_size = 10)
  return(rw)
}

find_don_shared = function(df, col){
  shared = (df
            %>% rename(feature = {{ col }})
            %>% select(Donor, feature)
            %>% na.omit()
            %>% unique()
            %>% count(feature)
            %>% filter(n > 1)
            %>% pull(feature))
  return(shared)
}

get_counts = function(df, both_vect, taxlev, ftype, size_df){

  ft_count = (df
              %>% rename(feature = {{ taxlev }})
              %>% filter(feature %in% both_vect)
              %>% count(Sample, Donor, feature, name = 'count')
              %>% pivot_wider(names_from = 'feature', values_from = 'count',
                              values_fill = 0)
              %>% pivot_longer(-(Sample:Donor),
                               names_to = 'feature', values_to = 'count')
              %>% left_join((size_df
                             %>% rename(feature = {{ taxlev }})),
                             by = c('feature', 'Donor'))
              %>% mutate(FeatureType = ftype,
                         TaxLev = taxlev))
  return(ft_count)
}
