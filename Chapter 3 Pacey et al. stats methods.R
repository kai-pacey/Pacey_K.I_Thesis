##### Estimating available Acropora biomass using weight size relationships: r code summary 

library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
#install.packages("rstanarm")
library(rstanarm)   
library(brms)       
library(coda)       
library(bayesplot)  
library(rstan)      
library(emmeans)    
library(broom)      
library(DHARMa)     
library(tidybayes)  
library(ggeffects)  
library(tidyverse)  
library(broom.mixed)  


#### initial models all data 

alldata %>% 
  filter(dia.1.cm < 110)

prior1 <- prior(normal(1, 1), nlpar = "b1", lb = 0) + 
  prior(normal(1, 1), nlpar = "b2", lb = 0) 

intial.log  <- brm(bf(weight.g ~ b1*dia.1.cm^b2 , 
                      b1 ~ 1,
                      b2 ~ 1,
                      nl = TRUE),
                   data = alldata %>% 
                     filter(dia.1.cm < 110), 
                   prior = prior1, 
                   iter = 12000, 
                   warmup = 5000,
                   thin = 12, 
                   cores = 4,
                   chains = 4,
                   family = lognormal(),
                   #control = list(adapt_delta = 0.99,
                   #              max_treedepth = 12), 
                   file = "initial.log.p.2-",
                   sample_prior = "yes",
                   save_all_pars = T, 
                   seed = 1234
)

plot(conditional_effects(intial.log), points = T)
pp_check(intial.log)
summary(intial.log)
prior_summary(intial.log)
bayes_R2(intial.log)

### lognormal vs normal

prior1 <- prior(normal(3, 3), nlpar = "b1", lb = 0) + 
  prior(normal(2, 3), nlpar = "b2", lb = 0) 

initial.norm  <- brm(bf(weight.g ~ b1*dia.1.cm^b2 , 
                        b1 ~ 1,
                        b2 ~ 1,
                        #b3 ~ growth.form.6 * type, 
                        nl = TRUE),
                     data = alldata %>% 
                       filter(dia.1.cm < 110), 
                     prior = prior1, 
                     iter = 12000, 
                     warmup = 5000,
                     thin = 12, 
                     cores = 4,
                     chains = 4,
                     family = gaussian(),
                     #control = list(adapt_delta = 0.99,
                     #              max_treedepth = 12), 
                     file = "initial.norm.new.p.5",
                     sample_prior = "yes",
                     save_all_pars = T, 
                     seed = 1234
)

plot(conditional_effects(initial.norm), points = T)
pp_check(initial.norm)
summary(initial.norm)
prior_summary(initial.norm)
bayes_R2(initial.norm)


loo(initial.norm, intial.log)

### power vs linear vs exponential

prior1 <- prior(normal(5, 5), nlpar = "b1", lb = 0) + 
  prior(normal(0, 3), nlpar = "b2", lb = 0) 

all.data.exp  <- brm(bf(weight.g ~ b1*exp(b2*dia.1.cm) , 
                        b1 ~ 1,
                        b2 ~ 1 ,
                        #b3 ~ growth.form.6 * type, 
                        nl = TRUE),
                     data = alldata %>% 
                       filter(dia.1.cm < 110), 
                     prior = prior1, 
                     iter = 12000, 
                     warmup = 5000,
                     thin = 12, 
                     cores = 4,
                     chains = 4,
                     family = lognormal(),
                     #control = list(adapt_delta = 0.99,
                     #              max_treedepth = 12), 
                     file = "all.data.exp.p.10",
                     sample_prior = "yes",
                     save_all_pars = T, 
                     seed = 1234
)

plot(conditional_effects(all.data.exp), points = T)
pp_check(all.data.exp)
summary(all.data.exp)
prior_summary(all.data.exp)
bayes_R2(all.data.exp)

all.data.linear  <- brm(bf(weight.g ~ dia.1.cm ),
                        data = alldata %>% 
                          filter(dia.1.cm < 110), 
                        iter = 12000, 
                        warmup = 5000,
                        thin = 12, 
                        cores = 4,
                        chains = 4,
                        family = lognormal(),
                        #control = list(adapt_delta = 0.99,
                        #              max_treedepth = 12), 
                        file = "all.data.linear.p.6",
                        sample_prior = "yes",
                        save_all_pars = T,
                        seed = 1234
)

plot(conditional_effects(all.data.linear), points = T)
pp_check(all.data.linear)
summary(all.data.linear)
prior_summary(all.data.linear)
bayes_R2(all.data.linear)


loo(intial.log, all.data.linear, all.data.exp)

### random slope vs simple 

##### comparison model for loo

prior1 <- prior(normal(1, 1), nlpar = "b1", lb = 0) + 
  prior(normal(1, 1), nlpar = "b2", lb = 0) +
  prior(normal(1, 1), nlpar = "b3", lb = 0) 

alldata.rsl.brms.s$data  <- brm(bf(weight.g ~ b1*b2^b3 , 
                              b1 ~ 1,
                              b2 ~ dia.1.cm * (growth.form.6:dia.1.cm) * (type:dia.1.cm),
                              b3 ~ 1,
                              #b3 ~ growth.form.6 * type, 
                              nl = TRUE),
                           data = alldata %>% 
                             filter(dia.1.cm < 110), 
                           prior = prior1, 
                           iter = 12000, 
                           warmup = 5000,
                           thin = 12, 
                           cores = 4,
                           chains = 4,
                           family = lognormal(),
                           #control = list(adapt_delta = 0.99,
                           #              max_treedepth = 12), 
                           file = "alldata.rsl.brms.s.p.1",
                           sample_prior = "yes",
                           save_all_pars = T, 
                           seed = 1234
)

plot(conditional_effects(alldata.rsl.brms.s), points = T)
summary(alldata.rsl.brms)
prior_summary(alldata.rsl.brms.s)
bayes_R2(alldata.rsl.brms)
pp_check(alldata.rsl.brms, nsamples = 1000, type = "ecdf_overlay")

hypothesis(alldata.rsl.brms.s, "b2_dia.1.cm:growth.form.6Staghorn = b2_dia.1.cm:growth.form.6Bottlebrush")

hypothesis(alldata.rsl.brms.s, "b2_dia.1.cm:growth.form.6Caespitose = b2_dia.1.cm:growth.form.6Corymbose")


prior1 <- prior(normal(1, 1), nlpar = "b1", lb = 0) + 
  prior(normal(1, 1), nlpar = "b2", lb = 0) +
  prior(normal(1, 1), nlpar = "b3", lb = 0) 

alldata.rsl.brms.s.COMP5_30  <- brm(bf(weight.g ~ b1*b2^b3 , 
                              b1 ~ 1,
                              b2 ~ dia.1.cm * (growth.form.6:dia.1.cm) * (type:dia.1.cm),
                              b3 ~ 1,
                              #b3 ~ growth.form.6 * type, 
                              nl = TRUE),
                           data = alldata.rsl.brms.s$data %>% 
                             filter(dia.1.cm < 110) %>% 
                             filter(dia.1.cm <30 & dia.1.cm > 5), 
                           prior = prior1, 
                           iter = 12000, 
                           warmup = 5000,
                           thin = 12, 
                           cores = 4,
                           chains = 4,
                           family = lognormal(),
                           #control = list(adapt_delta = 0.99,
                           #              max_treedepth = 12), 
                           file = "alldata.rsl.brms.s.p.1.5_30",
                           sample_prior = "yes",
                           save_all_pars = T, 
                           seed = 1234
)


hypothesis(alldata.rsl.brms.s.COMP5_30, "b2_dia.1.cm:growth.form.6Caespitose = b2_dia.1.cm:growth.form.6Corymbose")

#alldata.rsl.brms.s$data %>% write_csv("all.data.model.csv")


intial.log # <- see above 

alldata.rsl.brms.comp.2.s  <- brm(bf(weight.g ~ b1*b2^b3 , 
                                     b1 ~ 1,
                                     b2 ~ dia.1.cm * (growth.form.6:dia.1.cm),
                                     b3 ~ 1,
                                     #b3 ~ growth.form.6 * type, 
                                     nl = TRUE),
                                  data = alldata %>% 
                                    filter(dia.1.cm < 110), 
                                  prior = prior1, 
                                  iter = 12000, 
                                  warmup = 5000,
                                  thin = 12, 
                                  cores = 4,
                                  chains = 4,
                                  family = lognormal(),
                                  #control = list(adapt_delta = 0.99,
                                  #              max_treedepth = 12), 
                                  file = "alldata.rsl.brms.comp.2.s.p.1",
                                  sample_prior = "yes",
                                  save_all_pars = T, 
                                  seed = 1234
)

plot(conditional_effects(alldata.rsl.brms.comp.2.s), points = T)


alldata.rsl.brms.comp.3.s  <- brm(bf(weight.g ~ b1*b2^b3 , 
                                     b1 ~ 1,
                                     b2 ~ dia.1.cm * (type:dia.1.cm),
                                     b3 ~ 1,
                                     #b3 ~ growth.form.6 * type, 
                                     nl = TRUE),
                                  data = alldata %>% 
                                    filter(dia.1.cm < 110), 
                                  prior = prior1, 
                                  iter = 12000, 
                                  warmup = 5000,
                                  thin = 12, 
                                  cores = 4,
                                  chains = 4,
                                  family = lognormal(),
                                  #control = list(adapt_delta = 0.99,
                                  #              max_treedepth = 12), 
                                  file = "alldata.rsl.brms.comp.3.s.p",
                                  sample_prior = "yes",
                                  save_all_pars = T, 
                                  seed = 1234
)


loo(intial.log, alldata.rsl.brms.comp.2.s, alldata.rsl.brms.comp.3.s, alldata.rsl.brms.s)


#### Weight size models 

### Corymbose 

prior1 <- prior(normal(0.1, 0.1), nlpar = "b1", lb = 0) + 
  prior(normal(0.5, 0.1), nlpar = "b2", lb = 0) 



corymbose.axb.p <- brm(bf(weight.g ~ b1*dia.1.cm^b2, 
                          b1 ~ 1,
                          b2 ~ 1,
                          nl = TRUE),
                       data = CaesCombBot.data.wfrag, 
                       prior = prior1, 
                       iter = 12000, 
                       warmup = 5000,
                       thin = 12, 
                       cores = 4,
                       family = lognormal(),
                       #control = list(adapt_delta = 0.99,
                       #               max_treedepth = 12), 
                       file = "corymbose.axb.p.23",
                       sample_prior = "yes",
                       save_all_pars = T, 
                       seed = 1234
)

corymbose.axb.p.l5_30 <- brm(bf(weight.g ~ b1*dia.1.cm^b2, 
                          b1 ~ 1,
                          b2 ~ 1,
                          nl = TRUE),
                       data = corymbose.axb.p$data %>% filter(dia.1.cm < 30 & dia.1.cm > 5), 
                       prior = prior1, 
                       iter = 12000, 
                       warmup = 5000,
                       thin = 12, 
                       cores = 4,
                       family = lognormal(),
                       #control = list(adapt_delta = 0.99,
                       #               max_treedepth = 12), 
                       file = "corymbose.axb.5_35",
                       sample_prior = "yes",
                       save_all_pars = T, 
                       seed = 1234
)


posterior::as_draws_rvars(corymbose.axb.p$fit)

prior1 <- prior(normal(0.1, 0.1), nlpar = "b1", lb = 0) + 
  prior(normal(3, 5), nlpar = "b2", lb = 0) +
  prior(normal(0.5, 0.1), nlpar = "b3", lb = 0) 

corymbose.axb.p.b3 <- brm(bf(weight.g ~ b1*b2^b3, 
                             b1 ~ 1,
                             b2 ~ dia.1.cm + (type:dia.1.cm),
                             b3 ~ 1,
                             nl = TRUE),
                          data = CaesCombBot.data.wfrag, 
                          prior = prior1, 
                          iter = 12000, 
                          warmup = 5000,
                          thin = 12, 
                          cores = 4,
                          family = lognormal(),
                          #control = list(adapt_delta = 0.99,
                          #               max_treedepth = 12), 
                          file = "corymbose.axb.p.b3",
                          sample_prior = "yes",
                          save_all_pars = T, 
                          seed = 1234
)



plot(conditional_effects(corymbose.axb.p.b3), points = T)
summary(corymbose.axb.p)
prior_summary(corymbose.axb.p)
bayes_R2(corymbose.axb.p)
pp_check(corymbose.axb.p, nsamples = 1000, type = "ecdf_overlay")

curve(exp(1.90*x^0.4), from = 0, to = 50)

### Digitate 

prior1 <- prior(normal(0.5, 0.5), nlpar = "b1", lb = 0) + 
  prior(normal(0.7, 0.5), nlpar = "b2", lb = 0) 

digitate.axb.p <- brm(bf(weight.g ~ b1*dia.1.cm^b2, 
                         b1 ~ 1,
                         b2 ~ 1,
                         #b3 ~ 1 + growth.form.6,
                         nl = TRUE),
                      data = Digitate.col.wfrag, 
                      prior = prior1, 
                      iter = 12000, 
                      warmup = 5000,
                      thin = 12, 
                      cores = 4,
                      family = "lognormal",
                      #control = list(adapt_delta = 0.99, 
                      #              max_treedepth = 12), 
                      file = "digitate.axb.p.14",
                      sample_prior = "yes",
                      save_all_pars = T, 
                      seed = 1234
)

plot(conditional_effects(digitate.axb.p), points = T)
summary(digitate.axb.p)
 
bayes_R2(digitate.axb.p)
pp_check(digitate.axb.p, nsamples = 1000, type = "ecdf_overlay")

## est
posterior::as_draws_rvars(digitate.axb.p$fit)


prior1 <- prior(normal(0.1, 0.1), nlpar = "b1", lb = 0) + 
  prior(normal(3, 5), nlpar = "b2", lb = 0) +
  prior(normal(0.5, 0.1), nlpar = "b3", lb = 0) 

digitate.axb.p.b3 <- brm(bf(weight.g ~ b1*b2^b3, 
                            b1 ~ 1,
                            b2 ~ dia.1.cm + (type:dia.1.cm),
                            b3 ~ 1,
                            nl = TRUE),
                         data = Digitate.col.wfrag, 
                         prior = prior1, 
                         iter = 12000, 
                         warmup = 5000,
                         thin = 12, 
                         cores = 4,
                         family = "lognormal",
                         #control = list(adapt_delta = 0.99, 
                         #              max_treedepth = 12), 
                         file = "digitate.axb.p.14.b3",
                         sample_prior = "yes",
                         save_all_pars = T, 
                         seed = 1234
)

plot(conditional_effects(digitate.axb.p.b3), points = T)


### Staghorn 

prior1 <- prior(normal(0.7, 0.5), nlpar = "b1", lb = 0) + 
  prior(normal(1, 0.7), nlpar = "b2", lb = 0) 

staghorn.axb.p <- brm(bf(weight.g ~ b1*dia.1.cm^b2, 
                         b1 ~ 1,
                         b2 ~ 1,
                         #b3 ~ 1 + growth.form.6,
                         nl = TRUE),
                      data = staghorn.col.only.dat.wfrag, 
                      prior = prior1, 
                      iter = 12000, 
                      warmup = 5000,
                      thin = 12, 
                      cores = 4,
                      family = lognormal(),
                      #control = list(adapt_delta = 0.99, 
                      #              max_treedepth = 12), 
                      file = "staghorn.axb.p.2",
                      sample_prior = "yes", 
                      seed = 1234
)

plot(conditional_effects(staghorn.axb.p), points = T)
summary(staghorn.axb.p)
prior_summary(staghorn.axb.p)
bayes_R2(staghorn.axb.p)
pp_check(staghorn.axb.p, nsamples = 1000, type = "ecdf_overlay")

## est
posterior::as_draws_rvars(staghorn.axb.p$fit)


prior1 <- prior(normal(0.1, 0.1), nlpar = "b1", lb = 0) + 
  prior(normal(3, 5), nlpar = "b2", lb = 0) +
  prior(normal(0.5, 0.1), nlpar = "b3", lb = 0)

staghorn.axb.p.b3 <- brm(bf(weight.g ~ b1*b2^b3, 
                            b1 ~ 1,
                            b2 ~ dia.1.cm + (type:dia.1.cm),
                            b3 ~ 1,
                            nl = TRUE),
                         data = staghorn.col.only.dat.wfrag, 
                         prior = prior1, 
                         iter = 12000, 
                         warmup = 5000,
                         thin = 12, 
                         cores = 4,
                         family = lognormal(),
                         #control = list(adapt_delta = 0.99, 
                         #              max_treedepth = 12), 
                         file = "staghorn.axb.p.2.b3",
                         sample_prior = "yes", 
                         seed = 1234
)

plot(conditional_effects(staghorn.axb.p.b3), points = T)

weight.size.mode.plot.stag <- staghorn.col.only.dat.wfrag %>%
  modelr::data_grid(dia.1.cm = seq_range(dia.1.cm, n = 508)) %>%
  add_epred_draws(staghorn.axb.p) %>%
  ggplot(aes(x = dia.1.cm, y = .epred/1000)) +
  stat_lineribbon(aes(y = .epred/1000), colour = "darkred", alpha = 0.9, .width = 0.95) +
  geom_point(data = staghorn.col.only.dat.wfrag, aes(dia.1.cm, weight.g/1000), alpha = 0.45, colour = "grey36", size = 3) +
  stat_lineribbon(aes(y = .epred/1000), colour = "darkred", alpha = 0.95, geom = "line", size = 1.5, .width = 0.95) +
  scale_fill_brewer(palette = "Reds") +
  scale_color_brewer(palette = "Set2") + 
  theme_classic() + 
  labs(x = "Maximum diameter (cm)", y = "Coral piece weight (kg)", fill = "Prob. level") + 
  theme(axis.text = element_text(size = 21, colour = "black"),
        axis.title = element_text(size = 22, colour = "black"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 17), 
        legend.position = "none", 
        legend.box.background = element_rect(colour = "black", size = 1.1), 
        panel.border = element_rect(colour = "black", fill=NA, size = 1.5)) + 
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(seq(from = 3, to = 30, by = 3))) + 
  scale_x_continuous(expand = c(0,0), 
                     #limits = c(0, max(CaesCombBot.data$dia.1.cm)),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) + 
  coord_cartesian(xlim = c(0, max(staghorn.col.only.dat.wfrag$dia.1.cm))) +
  annotate("text", x = 2.5, y = c(28, 28*((0.95-0.03)), 28*((0.90-0.05)), 28*(0.85-0.08)), label = staghorn.label, size = 8, hjust = 0, parse = T) #+ 
#annotate("text", x = 2.5, y = 28000, label = paste("y %~% italic(lognormal)~(`1.53` + `2.57`~italic(x)^0.43)"), size = 5, parse = T, hjust = 0) +
#annotate("text", x = 2.5, y = 27000, label = "italic(n)==508", size = 5, hjust = 0, parse = T) 
weight.size.mode.plot.stag + ggsave("Staghorn weight_size plot 095.pdf", units = "cm", width = 20, height = 22)

weight.size.mode.plot.stag +
  stat_function(fun = function(x) exp(1.58*x^0.42)/1000, size = 2, colour = "black", linetype = "longdash", alpha = 0.65)

posterior::as_draws_rvars(staghorn.axb.p$fit)

### Table 
prior1 <- prior(normal(0.5, 0.5), nlpar = "b1", lb = 0) + 
  prior(normal(0.7, 0.5), nlpar = "b2", lb = 0) 

table.axb.p <- brm(bf(weight.g ~ b1*dia.1.cm^b2, 
                      b1 ~ 1,
                      b2 ~ 1,
                      #b3 ~ 1 + growth.form.6,
                      nl = TRUE),
                   data = table.wfrag, 
                   prior = prior1, 
                   iter = 12000, 
                   warmup = 5000,
                   thin = 12, 
                   cores = 4,
                   family = lognormal(),
                   #control = list(adapt_delta = 0.99,
                   #              max_treedepth = 12), 
                   file = "table.axb.p.3",
                   sample_prior = "yes", 
                   seed = 1234
)

plot(conditional_effects(table.axb.p), points = TRUE)
summary(table.axb.p) 
prior_summary(table.axb.p)
bayes_R2(table.axb.p)
pp_check(staghorn.axb.p, nsamples = 1000, type = "ecdf_overlay")

## est
posterior::as_draws_rvars(table.axb.p$fit)

?as_draws_rvars()

prior1 <- prior(normal(0.1, 0.1), nlpar = "b1", lb = 0) + 
  prior(normal(3, 5), nlpar = "b2", lb = 0) +
  prior(normal(0.5, 0.1), nlpar = "b3", lb = 0)

table.axb.p.b3 <- brm(bf(weight.g ~ b1*b2^b3, 
                         b1 ~ 1,
                         b2 ~ dia.1.cm + (type:dia.1.cm),
                         b3 ~ 1,
                         nl = TRUE),
                      data = table.wfrag, 
                      prior = prior1, 
                      iter = 12000, 
                      warmup = 5000,
                      thin = 12, 
                      cores = 4,
                      family = lognormal(),
                      #control = list(adapt_delta = 0.99,
                      #              max_treedepth = 12), 
                      file = "table.axb.p.3.b3",
                      sample_prior = "yes", 
                      seed = 1234
)

plot(conditional_effects(table.axb.p.b3), points = TRUE)



weight.size.mode.plot.tabl <- table.wfrag %>%
  modelr::data_grid(dia.1.cm = seq_range(dia.1.cm, n = 599)) %>%
  add_epred_draws(table.axb.p) %>%
  ggplot(aes(x = dia.1.cm, y = .epred/1000)) +
  stat_lineribbon(aes(y = .epred/1000), colour = "purple4", alpha = 0.9, .width = 0.95) +
  geom_point(data = table.wfrag, aes(dia.1.cm, weight.g/1000), alpha = 0.45, colour = "grey36", size = 3) +
  stat_lineribbon(aes(y = .epred/1000), colour = "purple4", alpha = 0.95, geom = "line", size = 1.5, .width = 0.95) +
  scale_fill_brewer(palette = "Purples") +
  scale_color_brewer(palette = "Set2") + 
  theme_classic() + 
  labs(x = "Maximum diameter (cm)", y = "", fill = "Prob. level") + 
  theme(axis.text = element_text(size = 21, colour = "black"),
        axis.title = element_text(size = 22, colour = "black"), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 17), 
        legend.position = "none", 
        legend.box.background = element_rect(colour = "black", size = 1.1), 
        panel.border = element_rect(colour = "black", fill=NA, size = 1.5)) + 
  scale_y_continuous(expand = c(0,0), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)) + 
  scale_x_continuous(expand = c(0,0), 
                     #limits = c(0, max(CaesCombBot.data$dia.1.cm)),
                     breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  coord_cartesian(xlim = c(0, max(table.wfrag$dia.1.cm))) +
  annotate("text", x = 2.5, y = c(53.2, 53.2*((0.95-0.03)), 53.2*((0.90-0.05)), 53.2*(0.85-0.08)), label = table.label, size = 8, hjust = 0, parse = T) #+ 
# annotate("text", x = 2.5, y = 53000, label = paste("y %~% italic(lognormal)~(`4.33` + `7.71`~italic(x)^0.36)"), size = 6, parse = T, hjust = 0) +
#annotate("text", x = 2.5, y = 51000, label = "italic(n)==599", size = 6, hjust = 0, parse = T) 
weight.size.mode.plot.tabl + ggsave("table weight_size plot.2 095.pdf", units = "cm", width = 20, height = 22)


weight.size.mode.plot.tabl + 
  stat_function(fun = function(x) exp(2.13*x^0.35)/1000, size = 2, colour = "black", linetype = "longdash", alpha = 0.65)

#### Estimating biomass 

##### Corymbose/caespitose #####

transect.data <- read_csv("transects data.csv", col_types = "ccccccncccccccnncc") 

unique(transect.data$growth.form.6)

corymbose.transect <- transect.data %>% 
  filter(growth.form.6 == "Corymbose" | growth.form.6 ==  "Caespitose") %>% 
  mutate(dia.1.cm = pmax(diameter.1, diameter.2)) %>% 
  drop_na(dia.1.cm)

### modelr 

library(modelr)

corymbose.weight.prediction.new <- add_predictions(corymbose.transect, corymbose.axb.p, var = "weight.g.est.data" )

corymbose.weight.prediction.1.new <- corymbose.weight.prediction.new %>% 
  select(region, reef, site, transect, depth, type, growth.form.6, dia.1.cm, weight.g.est.data)

### quick graph 

#ggplot(corymbose.weight.prediction.1.new, aes(weight.g.est.data[,"Estimate"], group = reef)) + geom_histogram()

### Redo estimate using model equation 

corymbose.weight.prediction.1.new.eq <- corymbose.weight.prediction.1.new %>% 
  mutate(eq.weight.g = exp(1.90*dia.1.cm^0.40))



##### Digitate ##### 

digitate.transect <- transect.data %>% 
  filter(growth.form.6 == "Digitate") %>% 
  mutate(dia.1.cm = pmax(diameter.1, diameter.2)) %>% 
  drop_na(dia.1.cm)

### modelr 

library(modelr)

digitate.weight.prediction.new <- add_predictions(digitate.transect, digitate.axb.p, var = "weight.g.est.data" )

digitate.weight.prediction.1.new <- digitate.weight.prediction.new %>% 
  select(region, reef, site, transect, depth, type, growth.form.6, dia.1.cm, weight.g.est.data)

#View(digitate.weight.prediction.1)

digitate.weight.prediction.1.new.eq <- digitate.weight.prediction.1.new %>% 
  mutate(eq.weight.g = exp(1.9*dia.1.cm^0.42))

#View(digitate.weight.prediction.1.new.eq)

##### Staghorn #####

staghorn.transect <- transect.data %>% 
  filter(growth.form.6 == "Staghorn" | growth.form.6 == "Bottlebrush" ) %>% 
  mutate(dia.1.cm = pmax(diameter.1, diameter.2)) %>% 
  drop_na(dia.1.cm)

### modelr 


staghorn.weight.prediction.new <- add_predictions(staghorn.transect, staghorn.axb.p, var = "weight.g.est.data" )

staghorn.weight.prediction.1.new <- staghorn.weight.prediction.new %>% 
  select(region, reef, site, transect, depth, type, growth.form.6, dia.1.cm, weight.g.est.data)

#View(staghorn.weight.prediction.1)

staghorn.weight.prediction.1.new.eq <- staghorn.weight.prediction.1.new %>% 
  mutate(eq.weight.g = exp(1.6*dia.1.cm^0.42))

#View(staghorn.weight.prediction.1.new.eq)

##### Table ##### 

table.transect <- transect.data %>% 
  filter(growth.form.6 == "Table") %>% 
  mutate(dia.1.cm = pmax(diameter.1, diameter.2)) %>% 
  drop_na(dia.1.cm)

### modelr 

table.weight.prediction.new <- add_predictions(table.transect, table.axb.p, var = "weight.g.est.data" )

table.weight.prediction.1.new <- table.weight.prediction.new %>% 
  select(region, reef, site, transect, depth, type, growth.form.6, dia.1.cm, weight.g.est.data)

#View(table.weight.prediction.1)

table.weight.prediction.1.new.eq <- table.weight.prediction.1.new %>% 
  mutate(eq.weight.g = exp(2.10*dia.1.cm^0.35))

#View(table.weight.prediction.1.new.eq)

##### Combining sets #####

transect.data.est.wfrag.eq.prep <- bind_rows(tibble(corymbose.weight.prediction.1.new.eq), 
                                             as_tibble(digitate.weight.prediction.1.new.eq), 
                                             as_tibble(staghorn.weight.prediction.1.new.eq), 
                                             as_tibble(table.weight.prediction.1.new.eq))



weight.g.est.data <- transect.data.est.wfrag.eq.prep$weight.g.est.data %>% 
  as_tibble() %>% 
  rename(model.est = V1,
         model.se = V2, 
         lower.est = V3,
         upper.est = V4)

transect.data.est.wfrag.eq <- transect.data.est.wfrag.eq.prep %>% 
  select(-weight.g.est.data) %>% 
  tibble(weight.g.est.data) %>% 
  select(-model.se, -lower.est, -upper.est)

mean(transect.data.est.wfrag.eq$eq.weight.g)
mean(transect.data.est.wfrag.eq$model.est)
mean(transect.data.est.wfrag.eq$model.est.old)

### Table 2



library(plotrix)

error.eq <- transect.data.est.wfrag.eq %>% 
  mutate("growth.form.new" = recode(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn"), 
         "transect.new" = recode(transect, 
                                 `1(p1)` = "1", 
                                 `1(p2)` = "1")) %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  ungroup() %>%  
  group_by(region, growth.form.new, reef, site, transect.new) %>% 
  summarise("mean.weight.equation" = mean(eq.weight.g),
            "mean.weight.model" = mean(model.est),
            "weight.se.equation" = std.error(eq.weight.g), 
            "weight.se.model" = std.error(model.est)) %>% 
  mutate("av.weight.pm.pt.equation" = mean.weight.equation/50,
         "av.weight.pm.pt.model" = mean.weight.model/50)

## total biomass

total.bio <- transect.data.est.wfrag.eq %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  filter(!dia.1.cm > 110) %>% 
  #filter(ID<5336) %>% 
  mutate(growth.form.new = recode_factor(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn")) %>% 
  ungroup() %>% 
  summarise(sum = sum(eq.weight.g),
    min = min(eq.weight.g), 
            mean = mean(eq.weight.g), 
            std.er = std.error(eq.weight.g),
            max = max(eq.weight.g), 
            n = n(), 
            group = "biomass (g)")  %>% 
  arrange(-mean) %>% 
  relocate(min, .before = max)

View(total.bio)


### summary table

library(plotrix)



table.summary.trans.g <- transect.data.est.wfrag.eq %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  filter(!dia.1.cm > 110) %>% 
  #filter(ID<5336) %>% 
  mutate(growth.form.new = recode_factor(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn")) %>% 
  ungroup() %>% 
  #group_by(growth.form.new) %>% 
  summarise(min = min(eq.weight.g), 
            mean = mean(eq.weight.g), 
            std.er = std.error(eq.weight.g),
            max = max(eq.weight.g), 
            n = n(), 
            group = "biomass (g)")  %>% 
  arrange(-mean) %>% 
  relocate(min, .before = max)

View(table.summary.trans.g)

table.summary.trans.cm <- errorbar.data %>% 
  filter(estimate.g > 0 & dia.1.cm > 0) %>% 
  filter(!dia.1.cm > 110) %>% 
  #filter(ID<5336) %>%   
  mutate(growth.form.new = recode_factor(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn")) %>% 
  ungroup() %>% 
  group_by(growth.form.new) %>% 
  summarise(min = min(dia.1.cm),
            mean = mean(dia.1.cm), 
            std.er = std.error(dia.1.cm),
            max = max(dia.1.cm), 
            n = n(), 
            group = "diameter (cm)")

#### Biomass estimate for Stucco 

transect.data.est.wfrag.eq %>% 
  mutate("growth.form.new" = recode(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn"), 
         "transect.new" = recode(transect, 
                                 `1(p1)` = "1", 
                                 `1(p2)` = "1")) %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  ungroup() %>%  
  group_by(region, growth.form.new, reef) %>% 
  summarise("mean.weight.equation" = mean(eq.weight.g),
            "mean.weight.model" = mean(model.est),
            "weight.se.equation" = std.error(eq.weight.g), 
            "weight.se.model" = std.error(model.est)) %>% 
  mutate("av.weight.pm.pt.equation" = mean.weight.equation/50,
         "av.weight.pm.pt.model" = mean.weight.equation/50) %>% 
  filter(reef == "Stucco Site")

transect.data.est.wfrag.eq %>% 
  mutate("growth.form.new" = recode(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn"), 
         "transect.new" = recode(transect, 
                                 `1(p1)` = "1", 
                                 `1(p2)` = "1")) %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  ungroup() %>%  
  group_by(reef) %>% 
  summarise("mean.weight.equation" = mean(eq.weight.g),
            "mean.weight.model" = mean(model.est),
            "weight.se.equation" = std.error(eq.weight.g), 
            "weight.se.model" = std.error(model.est)) %>% 
  mutate("av.weight.pm.pt.equation" = mean.weight.equation/50,
         "av.weight.pm.pt.model" = mean.weight.model/50) %>% 
  filter(reef == "Stucco Site")


### square metre of Stucco Reef 
4740000

### all corals combined 
5.85 + 4.88 + 6.48 + 5.25

22.46 * 4740000
106460400/1000
106460.4/140000 

0.7604314 * 100

76.04314

### corymbose corals 
5.85

5.85 * 4740000
27729000/1000
27729/140000

0.1980643 * 100

19.80643

#### ANOVA 

### Model validation 

## replaced errorbar data
errorbar.model.data <- transect.data.est.wfrag.eq %>% 
  mutate("growth.form.new" = recode(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn"), 
         "transect.new" = recode(transect, 
                                 `1(p1)` = "1", 
                                 `1(p2)` = "1")) %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  ungroup() %>%  
  group_by(region, growth.form.new, reef, site, transect.new) %>% 
  summarise("mean.weight" = mean(eq.weight.g), 
            "weight.se" = std.error(eq.weight.g)) %>% 
  mutate("av.weight.pm.pt" = mean.weight/50) %>% 
  write_csv("error.bar.data.eq.csv")

unique(errorbar.model.data$reef)

stucco <- errorbar.model.data %>% 
  filter(reef == "Stucco Site") %>% 
  group_by(growth.form.new) %>% 
  summarise(mean = mean(av.weight.pm.pt))

model.1.norm <- brm(av.weight.pm.pt ~ growth.form.new * region + (1|reef/site), 
                    data = errorbar.model.data, 
                    iter = 12000, 
                    warmup = 5000, 
                    thin = 12, 
                    cores = 4,
                    family = gaussian(), 
                    save_pars = save_pars(all = T, group = T),
                    control = list(adapt_delta = 0.95), 
                    file = "model.1.norm.eq.2new", 
                    seed = 1234)

summary(model.1.norm.d.) 
pp_check(model.2.norm.d.1)
as.mcmc(model.1.norm.a.)
plot(as.mcmc(model.1.norm.a.))

View(errorbar.model.data)

model.2.lognorm <- brm(av.weight.pm.pt ~ growth.form.new * region + (1|reef/site), 
                       data = errorbar.model.data, 
                       iter = 12000, 
                       warmup = 5000, 
                       thin = 12, 
                       cores = 4,
                       family = lognormal(), 
                       save_pars = save_pars(all = T, group = T),
                       control = list(adapt_delta = 0.95#,
                                      #               max_treedepth = 12
                       ), 
                       file = "model.2.lognorm.eq.2new", 
                       seed = 1234)

plot(as.mcmc(model.2.lognorm))

loo(model.1.norm, model.2.lognorm) # lognormal 

model.3.region.s <- brm(av.weight.pm.pt ~ region + (1|reef/site), 
                        data = errorbar.model.data, 
                        iter =12000, 
                        warmup = 5000, 
                        thin = 12, 
                        cores = 4,
                        family = lognormal(), 
                        save_pars = save_pars(all = T, group = T),
                        control = list(adapt_delta = 0.95#,
                                       # max_treedepth = 12
                        ), 
                        file = "model.3.region.eq.2new", 
                        seed = 1234)

model.4.growthform.s$data <- brm(av.weight.pm.pt ~ growth.form.new + (1|reef/site), 
                            data = errorbar.model.data, 
                            iter = 12000, 
                            warmup = 5000, 
                            thin = 12, 
                            cores = 4,
                            family = lognormal(), 
                            save_pars = save_pars(all = T, group = T),
                            control = list(adapt_delta = 0.95#,
                                           #     max_treedepth = 12
                            ), 
                            file = "model.4.growthform.eq.2new", 
                            seed = 1234)

model.5.novary.s <- brm(av.weight.pm.pt ~ growth.form.new, 
                        data = errorbar.model.data, 
                        iter = 12000, 
                        warmup = 5000, 
                        thin = 12, 
                        cores = 4,
                        family = lognormal(), 
                        save_pars = save_pars(all = T, group = T),
                        control = list(adapt_delta = 0.95#,
                                       #      max_treedepth = 12
                        ), 
                        file = "model.5.novary.eq.2new", 
                        seed = 1234)

#install.packages("BayesFactor")
#library(BayesFactor)

model.null <- brm(av.weight.pm.pt ~ 1, 
                       data = errorbar.model.data, 
                       iter = 12000, 
                       warmup = 5000, 
                       thin = 12, 
                       cores = 4,
                       family = lognormal(), 
                       save_pars = save_pars(all = T, group = T),
                       control = list(adapt_delta = 0.95#,
                                      #               max_treedepth = 12
                       ), 
                       file = "model.anova.null", 
                       seed = 1234)

loo(model.1.norm, model.2.lognorm)

loo(model.2.lognorm, model.null)

LOO(model.2.lognorm, model.3.region.s, model.4.growthform.s, model.5.novary.s, model.null, moment_match = T)

#bayes_factor()

#model.5.novary$fit@stan_args

#loo(model.1.norm, model.2.lognorm, model.3.region.s, model.4.growthform.s, model.5.novary.s)

#bayes_factor( model.1.norm, model.2.lognorm, model.3.region.s, model.4.growthform.s, model.5.novary.s)

preds <- posterior_predict(model.2.lognorm, ndraws=1000, summary=F)
errorbar.model.data.resids <- createDHARMa(simulatedResponse = t(preds), 
                                           observedResponse = errorbar.model.data$av.weight.pm.pt,
                                           fittedPredictedResponse = apply(preds, 2, median), 
                                           integerResponse = T)

plot(errorbar.model.data.resids)
mcmc_trace(model.2.lognorm)


anova.region.em <- emmeans(model.2.lognorm, pairwise~ region)$contrast %>% 
  gather_emmeans_draws() %>% 
  mutate(Fit = .value)

anova.em <- emmeans(model.2.lognorm, pairwise~region + growth.form.new)$contrast %>% 
  gather_emmeans_draws() %>% 
  mutate(Fit = .value)

anova.growth.form.em <- emmeans(model.2.lognorm, pairwise~ growth.form.new)$contrast %>% 
  gather_emmeans_draws() %>% 
  mutate(Fit = .value)

View(anova.region.em %>%  
       group_by(contrast) %>% 
       summarise(P=sum(Fit>0)/n())) 

## swap model parameters 

#swap.data <- errorbar.model.data
#swap.data$region <- fct_relevel(errorbar.model.data$region, "Cairns", after = 1 )

#model.swap.1 <- brm(av.weight.pm.pt ~ growth.form.new * region + (1|reef/site), 
 #                   data = swap.data, 
  #                  iter = 50000, 
   #                 warmup = 22000, 
    #                thin = 50, 
     #               cores = 4,
      #              family = lognormal(), 
       #             save_pars = save_pars(all = T, group = T),
        #            control = list(adapt_delta = 0.95,
         #                          max_treedepth = 12
          #          ), 
           #         file = "model.2.swap.1", 
            #        seed = 1234)

#summary(model.swap.1)

#hypothesis(model.swap.1 , 'regionCairns < regionMackay')
#hypothesis(model.swap.1 , 'regionCairns > regionMackay')

#swap.data.2 <- errorbar.model.data
#swap.data.2$region <- fct_relevel(errorbar.model.data$region, "Mackay", after =  0)

#model.swap.2<- brm(av.weight.pm.pt ~ growth.form.new * region + (1|reef/site), 
 #                  data = swap.data.2, 
  #                 iter = 50000, 
   #                warmup = 22000, 
    #               thin = 50, 
     #              cores = 4,
      #             family = lognormal(), 
       #            save_pars = save_pars(all = T, group = T),
        #           control = list(adapt_delta = 0.95,
         #                         max_treedepth = 12
          #         ), 
           #        file = "model.2.swap.2", 
            #       seed = 1234)

#summary(model.swap.2)

#hypothesis(model.swap.2, 'regionLizardIsland < regionCairns')
#hypothesis(model.swap.2, 'regionLizardIsland > regionCairns')

#### acropora kg per km 2 calculations

errorbar.model.data.2 <- transect.data.est.wfrag.eq %>% 
  mutate("growth.form.new" = recode(growth.form.6, Caespitose = "Corymbose", Bottlebrush = "Staghorn"), 
         "transect.new" = recode(transect, 
                                 `1(p1)` = "1", 
                                 `1(p2)` = "1")) %>% 
  filter(eq.weight.g > 0 & dia.1.cm > 0) %>% 
  ungroup() %>%  
  group_by(region, reef, site, transect.new) %>% 
  summarise("mean.weight" = mean(eq.weight.g), 
            "weight.se" = std.error(eq.weight.g), 
            "av.weight.pm.pt" = mean.weight/50) %>%
  ungroup() %>% 
  summarise("av.weight.pm.pt" = mean.weight/50) 


errorbar.model.data.2 %>% 
  ungroup() %>% 
  summarise(sum = sum(av.weight.pm.pt))

### proportion fragments

alldata.rsl.brms.s$data %>% 
  group_by(type) %>% 
  summarise(n = n())

type                 n
<fct>            <int>
1 Fragment          3177
2 Colony             131
3 Colony with base   664

3177+131+664

#### HELP

errorbar.model.data %>% 
  
#### Weight 

  
  
  
  
