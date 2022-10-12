#################################
## IMPORT FUNCTIONS & PACKAGES ##
#################################
library(readr)
library(viridis)
library(R.filesets)
source('hawkes_functions.R')
source('hawkes_utilities.R')
#################
## IMPORT DATA ##
#################

# set magnitude of completeness and time domain
M0 = 2.99
T1 = 0; T2 = 357

# import data
list.input <- input.file.to.list('user_input.txt')
dd.ama <- list.input$catalog


# histogram plot
pl.hist <- ggplot(dd.ama, aes(x = time_date)) + 
  geom_histogram(binwidth = 60*60*24*7, color = 'black', fill = 'white') + 
  geom_star(data = dd.ama[dd.ama$Mw >= 5,], mapping = aes(x = time_date, y = 0), size = 2, fill = 'red') + 
  xlab('time') + 
  ylab('counts') + 
  annotate('text', x = dd.ama$time_date[nrow(dd.ama) - 5], y = 150, label = '(a)') + 
  theme_classic() 

# time vs magnitude scatter plot
pl.mgs <- ggplot(dd.ama[dd.ama$Mw < 5,], aes(x = time_date, y = Mw)) + 
  geom_point(size = 0.2) + 
  geom_star(data = dd.ama[dd.ama$Mw > 5,], mapping = aes(x = time_date, y = Mw), size = 2, fill = 'red') + 
  xlab('time') + 
  ylab('Mw') + 
  annotate('text', x = dd.ama$time_date[nrow(dd.ama) - 5], y = 5, label = '(b)') + 
  theme_classic() 


# calculate cumulative frequencies
CumSum <- sapply(T1:T2, \(x) sum(dd.ama$time.diff < x))
CumSum5 <- sapply(T1:T2, \(x) sum(dd.ama[dd.ama$Mw > 5, ]$time.diff < x))


# plot cumulative
pl.cum <- ggplot() + 
  geom_vline(data = dd.ama[dd.ama$Mw >= 5,], mapping = aes(xintercept = time.diff), 
             linetype = 2, color = 'darkgrey', size = 0.2) + 
  geom_line(aes(x = T1:T2, y = CumSum5*50), color = 'red', linetype = 2) + 
  geom_line(aes(x = T1:T2, y = CumSum)) + 
  geom_star(data = dd.ama[dd.ama$Mw >= 5,], mapping = aes(x = time.diff, y = 0), size = 2, fill = 'red') + 
  xlab('days') + 
  annotate('text', x = dd.ama$time.diff[nrow(dd.ama) - 5], y = 250, label = '(c)') +
  annotate('text', x = dd.ama$time.diff[nrow(dd.ama) - 5], y = 900, label = 'Mw > 3') +
  annotate('text', x = dd.ama$time.diff[nrow(dd.ama) - 5], y = 550, label = 'Mw > 5', color = 'red') +
  theme_classic() + 
  scale_y_continuous("cumulative count", 
                     sec.axis = sec_axis(~ . /50, name = 'cumulative count (Mw > 5)'))


# FIGURE 1
pdf(file = 'figure1.pdf')
multiplot(pl.hist, pl.mgs, pl.cum, layout = matrix(c(1,2,3,3), byrow = TRUE, ncol = 2))
dev.off()

########################
## MCMC MODEL FITTING ##
########################

time.st <- Sys.time()
# generate MCMC posterior samples
# change the sims argument to get different numbers of posterior samples
MCMC.ama <- sampleETASposterior(ts = dd.ama$time.diff, 
                                 magnitudes = dd.ama$Mw, 
                                 M0 = 2.99, 
                                 T=T2, sims = 10000, 
                                burnin = 5000, approx = TRUE)
Sys.time() - time.st
# save 
saveRDS(MCMC.ama, file = 'mcmc/MCMC.ama.samples.Rds')
MCMC.ama <- loadRDS('mcmc/MCMC.ama.samples.Rds')

##################################
## MCMC POSTERIOR DISTRIBUTIONS ##
##################################

# parameters name
par.names <- c('mu', 'K', 'alpha', 'c', 'p')

# extrat information for plotting purposes
post.MCMC <- data.frame(value = c(MCMC.ama[,1], MCMC.ama[,2], MCMC.ama[,3], MCMC.ama[,4], MCMC.ama[,5]),
                        param = rep(par.names, each = nrow(MCMC.ama)),
                        type = 'MCMC - posterior')

# data.frame containing priors
prior.MCMC <- rbind(data.frame(x = seq(min(MCMC.ama[,1]), max(MCMC.ama[,1]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dgamma(x, 0.1, 0.1),
                             param = 'mu'),
                    data.frame(x = seq(min(MCMC.ama[,2]), max(MCMC.ama[,2]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'K'),
                    data.frame(x = seq(min(MCMC.ama[,3]), max(MCMC.ama[,3]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'alpha'),
                    data.frame(x = seq(min(MCMC.ama[,4]), max(MCMC.ama[,4]),
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'c'),
                    data.frame(x = seq(min(MCMC.ama[,5]), max(MCMC.ama[,5]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 1, 10),
                             param = 'p')
) %>%
  mutate(type = 'MCMC - prior')



# FIGURE 7
pdf(file = 'figure7.pdf')
ggplot(post.MCMC, aes(x = value)) + geom_density(aes(color = type,
                                                     linetype = type)) + 
  geom_line(data = prior.MCMC, mapping = aes(x, pdf, color = type,
                                             linetype = type)) + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
dev.off()


###########################
## INLABRU MODEL FITTING ##
###########################

# data for Inlabru
data.bru <- data.frame(ts = dd.ama$time.diff, 
                       magnitudes = dd.ama$Mw) %>%
  mutate(idx.p = 1:nrow(dd.ama))

# gamma copula transformation
gamma.t <- function(x, a, b){
  bru_forward_transformation(qgamma, x, a, b)
}
# uniform copula transformation
unif.t <- function(x, a, b){
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t <- function(x, m, s){
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

##################################
## BINNING STRATEGY SENSITIVITY ##
##################################
# This is done considering the priors of the Inlabru replicate case

# code to find log-gaussian parameters to replicate the bayesianETAS prior for K

## find distribution of Kb as assumed by bayesianETAS priors
k.be = runif(1000000, min = 0, max = 10)
c.be = runif(1000000, min = 0, max = 10)
p.be = runif(1000000, min = 1, max = 10)

# respective
k.bru <- k.be/(c.be*(p.be-1))

# values
m.ln <- -1
sd.ln <- 2.03
# check if quantiles are close to each other
round(rbind(quantile(k.bru, c(0.01, 0.25, 0.5, 0.75, 0.99)),
            qlnorm(c(0.01, 0.25, 0.5, 0.75, 0.99), m.ln, sd.ln)), 4)
# mean and variance of the log-normal distribution
mean.ln(m.ln, sd.ln)
sd.ln(m.ln, sd.ln)

# set link functions for Inlabru replicate 
link.f.be <- list(mu = \(x) gamma.t(x, 0.1, 0.1), 
                  K = \(x) loggaus.t(x, -1, 2.03), 
                  alpha = \(x) unif.t(x, 0, 10), 
                  c_ = \(x) unif.t(x, 0, 10), 
                  p = \(x) unif.t(x, 1, 10))

# set data.frame with binning's parameters values
param.bin_exp <- expand.grid(coef_exp = c(1, 2, 3, 5, 7, 10),
                             N.exp = c(3, 10),
                             delta = c(0.1, 0.2, 0.5))


# set initial values for the parameters in the internal scale to get reasonable initial values (default is 0)
th.init <- list(th.mu = 0.5,
                th.K = 0.5,
                th.alpha = -2,
                th.c = -2,
                th.p = -2) 

# options for inlabru 
bru.opt.list <- list(bru_verbose = 4, # type of visual output 
                     bru_max_iter = 100, # maximum number of iterations
                     #bru_method = list(max_step = 0.5),
                     inla.mode = 'experimental', # type of inla algorithm
                     bru_initial = th.init) # parameters initial values

# fit the model for each binning's parameters combination
for(i in 1:nrow(param.bin_exp)){
  cat('coef.t. = ', param.bin_exp$coef_exp[i],
      'delta.t. = ', param.bin_exp$delta[i], 
      'N.max. = ', param.bin_exp$N.exp[i], '\n')
  # fit the model 
  fit_e <- Hawkes.bru(sample.s = data.bru, # data 
                      M0 = M0, # magnitude of completeness
                      T1 = 0, T2 = T2, # time domain
                      link.functions = link.f.be, # link functions
                      coef.t. = param.bin_exp$coef_exp[i], # binning parameter (delta)
                      delta.t. = param.bin_exp$delta[i], # binning parameter (Delta)
                      N.max. = param.bin_exp$N.exp[i], # binning parameter (n.max)
                      bru.opt = bru.opt.list) # bru options
  # save
  saveRDS(fit_e, file = paste0('fits/fit_binning/',
                               'fit_c',param.bin_exp$coef_exp[i],
                               '_d', param.bin_exp$delta[i],
                               '_N', param.bin_exp$N.exp[i], '.Rds') )
  cat('Completed: ', i/nrow(param.bin_exp), '\n')
}

# load the models
fit_bin_list <- lapply(list.files(path = 'fits/fit_binning/'), \(fore.path)
                       loadRDS(paste0('fits/fit_binning/', fore.path)))

# check convergence - Max deviation from previous: xx% of SD has to be smaller than 1
sapply(fit_bin_list, \(fit) tail(fit$bru_iinla$log, 3) )
# extract information about models (maximum number of iterations, time, if converged)
param.bin_exp$max.iter <- sapply(fit_bin_list, \(fit) max(fit$bru_iinla$track$iteration) )
param.bin_exp$time <- round(sapply(fit_bin_list, \(fit) 
                                   as.numeric(sum(fit$bru_iinla$timings$Time), units = 'mins')), 2)
param.bin_exp$converged <- param.bin_exp$max.iter < 100
# Table 5
param.bin_exp[order(param.bin_exp$time),]


# extract posterior distributions only for delta = 2
post.list <- lapply(which(param.bin_exp$coef_exp == 2), \(idx) 
                    extract.post.df(fit_bin_list[[idx]], link.f.be) %>%
                      mutate(c.exp = param.bin_exp$coef_exp[idx],
                             N.exp = param.bin_exp$N.exp[idx],
                             delta = param.bin_exp$delta[idx]))

# FIGURE 10
pdf('figure10.pdf')
ggplot(bind_rows(post.list) %>% 
         filter(param %in% c('K', 'c', 'p')), 
       aes(x,y,color = factor(delta), linetype = factor(N.exp))) + 
  geom_line() + 
  facet_wrap(facets = vars(param), scales = 'free') + 
  scale_y_log10() + 
  labs(color = ~ Delta, linetype = ~ n[max]) +
  xlab('value') + 
  ylab('log10(pdf)') +
  theme_classic()
dev.off()


# this is the Inlabru replicate model
fit_be <- fit_bin_list[[which(param.bin_exp$delta == 0.2 &
                                param.bin_exp$N.exp == 3 &
                                param.bin_exp$coef_exp == 2
)]]


# extract posterior distribution and set priors
bru.be.post <- extract.post.df(fit_be, link.f.be) %>%
  mutate(prior.type = 'be',
         binning = 'def',
         prior = case_when(param == 'mu' ~ dgamma(x,0.1,0.1),
                           param == 'K' ~ dlnorm(x, meanlog = -1, sdlog = 2.03),
                           param == 'alpha' ~ dunif(x, 0, 10),
                           param == 'c' ~ dunif(x, 0, 10),
                           param == 'p' ~ dunif(x, 1, 10)))


# FIGURE 8
pdf('figure8.pdf')
ggplot(bru.be.post, aes(x,y)) + 
  geom_line(aes(color = 'Inlabru rep - posterior',
                linetype = 'Inlabru rep - posterior')) +
  geom_line(aes(y = prior, color = 'Inlabru rep - prior',
                linetype = 'Inlabru rep - prior'))
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) + 
  xlab('value') + 
  ylab('log10(pdf)') +
  theme_classic()
dev.off()

#######################
## PRIOR SENSITIVITY ##
#######################

# logarithm standard deviation values
v.par <- c(1, 1.5, 2, 2.5)

# table 6
rbind(data.frame(sd.log = v.par[1], 
                 mean = mean(rlnorm(10000, meanlog = 0, sdlog = v.par[1])),
                 sd = sd(rlnorm(10000, meanlog = 0, sdlog = v.par[1])),
                 q0.025 = qlnorm(0.025, meanlog = 0, sdlog = v.par[1]),
                 q0.5 = qlnorm(0.5, meanlog = 0, sdlog = v.par[1]),
                 q0.975 = qlnorm(0.975, meanlog = 0, sdlog = v.par[1])),
      data.frame(sd.log = v.par[2], 
                 mean = mean(rlnorm(10000, meanlog = 0, sdlog = v.par[2])),
                 sd = sd(rlnorm(10000, meanlog = 0, sdlog = v.par[2])),
                 q0.025 = qlnorm(0.025, meanlog = 0, sdlog = v.par[2]),
                 q0.5 = qlnorm(0.5, meanlog = 0, sdlog = v.par[2]),
                 q0.975 = qlnorm(0.975, meanlog = 0, sdlog = v.par[2])),
      data.frame(sd.log = v.par[3], 
                 mean = mean(rlnorm(10000, meanlog = 0, sdlog = v.par[3])),
                 sd = sd(rlnorm(10000, meanlog = 0, sdlog = v.par[3])),
                 q0.025 = qlnorm(0.025, meanlog = 0, sdlog = v.par[3]),
                 q0.5 = qlnorm(0.5, meanlog = 0, sdlog = v.par[3]),
                 q0.975 = qlnorm(0.975, meanlog = 0, sdlog = v.par[3])),
      data.frame(sd.log = v.par[4], 
                 mean = mean(rlnorm(10000, meanlog = 0, sdlog = v.par[4])),
                 sd = sd(rlnorm(10000, meanlog = 0, sdlog = v.par[4])),
                 q0.025 = qlnorm(0.025, meanlog = 0, sdlog = v.par[4]),
                 q0.5 = qlnorm(0.5, meanlog = 0, sdlog = v.par[4]),
                 q0.975 = qlnorm(0.975, meanlog = 0, sdlog = v.par[4]))
      
) 

# create list of link functions
link.f.lg.list <- lapply(1:length(v.par), \(i)  
                         list(mu = \(x) loggaus.t(x, 0, force(v.par[i])), 
                              K = \(x) loggaus.t(x, 0, force(v.par[i])), 
                              alpha = \(x) loggaus.t(x, 0, force(v.par[i])), 
                              c_ = \(x) loggaus.t(x, 0, force(v.par[i])), 
                              p = \(x) 1 + loggaus.t(x, 0, force(v.par[i])))
)

# inverse transformation such that inv.lgaus.t(loggaus.t(x)) = x
inv.lgaus.t <- function(x, a, b){
  qnorm(plnorm(x, a, b))
}

# posterior mean
post.par.mean <- c(link.f.be$mu(fit_be$summary.fixed$mean[1]),
                   link.f.be$K(fit_be$summary.fixed$mean[2]),
                   link.f.be$alpha(fit_be$summary.fixed$mean[3]),
                   link.f.be$c_(fit_be$summary.fixed$mean[4]),
                   link.f.be$p(fit_be$summary.fixed$mean[5]))

# fit the model for different priors
for(i in 4:length(v.par)){
  # select the link functions
  link_f <- link.f.gamma.list[[i]]
  # select parameter
  logsd.par <- v.par[i]
  # set initial values (equals to the posterior mean)
  th.init_ <- list(th.mu = inv.lgaus.t(post.par.mean[1], 0, logsd.par),
                   th.K = inv.lgaus.t(post.par.mean[2], 0, logsd.par),
                   th.alpha = inv.lgaus.t(post.par.mean[3],0, logsd.par),
                   th.c = inv.lgaus.t(post.par.mean[4], 0, logsd.par),
                   th.p = inv.lgaus.t(post.par.mean[5] - 1, 0,  logsd.par))
  # fit the model
  fit_ <- Hawkes.bru(sample.s = data.bru, M0 = M0,
                      T1 = 0, T2 = T2, 
                      coef.t. = 2, delta.t. = 0.1, N.max. = 10, 
                      link.functions = link_f,
                      bru.opt = list(bru_verbose = 3,
                                     bru_max_iter = 110,
                                     inla.mode = 'experimental',
                                     bru_initial = th.init_))
  # store the model
  saveRDS(fit_, file = paste0('fits/fit_priors/lnorm/', 
                              'fit_par', logsd.par, '_N.max10.Rds') )
}

# this two steps are needed only if the user does not fit the models
# retrieve file names
l.files <- list.files(path = 'fits/fit_priors/lnorm/')

# load fitted models
fit.prior.list <- foreach(i = 1:length(l.files)) %do% {
  fit_ <- loadRDS(paste0('fits/fit_priors/lnorm/', l.files[i]))
  g.p <- parse_number(l.files[i])
  list(fit = fit_,
       logsd.par = g.p)
}

# extract posterior distributions
post.pr.list <- lapply(1:length(fit.prior.list), \(idx) 
                       extract.post.df(fit.prior.list[[idx]]$fit, 
                                       link.f.gamma.list[[idx]]) %>%
                         mutate(logsd.p = fit.prior.list[[idx]]$logsd.par))

# merge them for plotting
post.pr.bind <- bind_rows(post.pr.list) 

# add prior information
post.pr.bind <- post.pr.bind %>% 
  mutate(prior = case_when(param == 'mu' ~ dlnorm(x, 0, logsd.p),
                           param == 'K' ~ dlnorm(x, 0, logsd.p),
                           param == 'alpha' ~ dlnorm(x,0, logsd.p),
                           param == 'c' ~ dlnorm(x,0, logsd.p),
                           param == 'p' ~ dlnorm(x - 1, 0, logsd.p)))

# FIGURE 11
pdf('figure11.pdf')
ggplot(post.pr.bind, 
       aes(x,y, color = factor(gamma.p))) + 
  geom_line(aes(linetype = 'Posterior')) + 
  geom_line(aes(y = prior, linetype = 'Prior')) + 
  labs(color = ~ sigma[log], linetype = '') + 
  scale_y_log10() +
  ylab('log10(pdf)') + 
  xlab('value') + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  scale_color_viridis(discrete = TRUE)
dev.off()



########################
## INLABRU GAMMA CASE ##
########################

# set link functions
link.gamma <-   list(mu = \(x) gamma.t(x, 0.1, 1), 
                     K = \(x) gamma.t(x, 1, 0.5), 
                     alpha = \(x) gamma.t(x, 1, 0.5), 
                     c_ = \(x) gamma.t(x, 0.1, 1), 
                     p = \(x) 1 + gamma.t(x, 0.1, 0.5)) 
# fit the model
fit_gamma <- Hawkes.bru(sample.s = data.bru, M0 = M0,
                        T1 = 0, T2 = T2,
                        link.functions = link.gamma,
                        coef.t. = 2, delta.t. = 0.1, N.max. = 10,
                        bru.opt = bru.opt.list.gamma)
saveRDS(fit_gamma, file = 'fits/fit_gamma.Rds')
fit_gamma <- loadRDS('fits/fit_gamma.Rds')

##  Prior and posterior data.frame.
bru.gamma.post <- extract.post.df(fit_gamma, link.gamma) %>%
  mutate(prior.type = 'Inlabru - gamma',
         prior = case_when(param == 'mu' ~ dgamma(x, 0.1, 1),
                           param == 'K' ~ dgamma(x, 1, 0.5),
                           param == 'alpha' ~ dgamma(x, 1, 0.5),
                           param == 'c' ~ dgamma(x, 0.1, 1),
                           param == 'p' ~ dgamma(x - 1, 0.1, 0.5)))

## posterior comparison
bru.be.post <- bru.be.post %>%
  mutate(model = 'Inlabru rep',
         prior.n = 'Prior',
         post.n = 'Posterior')

bru.gamma.post <- bru.gamma.post %>%
  mutate(model = 'Inlabru gamma',
         prior.n = 'Prior',
         post.n = 'Posterior')


bru.post.bind <- bind_rows(bru.be.post, bru.gamma.post) 

# FIGURE 9 just remove or add scale_y_log10() to obtain the two.
pdf('figure9.pdf')
ggplot(bru.post.bind, aes(x)) + 
  geom_line(aes(y = y, color = model, linetype = post.n)) + 
  geom_line(aes(y = prior, color = model, linetype = prior.n)) +
  scale_y_log10() +
  facet_wrap(facets = vars(param), scales = 'free',
             labeller =  label_parsed) +
  xlab('value') + 
  ylab('log10(pdf)') + 
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
dev.off()

######################
## GOOODNESS-OF-FIT ##
######################

# extract 10000 values from the posterior of the parameters
# this is done in 10 batches of 1000 
# Inlabru replicate case
bru.sample.be <- foreach(i = 1:10, .combine = rbind) %do% {
  sample.p <- generate(fit_be, data.frame(1),  ~ data.frame(mu = link.f.be$mu(th.mu),
                                                            K = link.f.be$K(th.K),
                                                            alpha = link.f.be$alpha(th.alpha),
                                                            c = link.f.be$c_(th.c),
                                                            p = link.f.be$p(th.p)), 
                       n.samples = 1000)
  
  bind_rows(sample.p)
}
saveRDS(bru.be.sample, file = 'posterior_samples/bru.sample.be.Rds')
bru.be.sample <- loadRDS('posterior_samples/bru.sample.be.Rds')
# Inlabru gamma case
bru.sample.gamma <- foreach(i = 1:10, .combine = rbind) %do% {
  sample.p <- generate(fit_gamma, data.frame(1),  ~ data.frame(mu = link.gamma$mu(th.mu),
                                                               K = link.gamma$K(th.K),
                                                               alpha = link.gamma$alpha(th.alpha),
                                                               c = link.gamma$c_(th.c),
                                                               p = link.gamma$p(th.p)), 
                       n.samples = 1000)
  
  bind_rows(sample.p)
}
saveRDS(bru.sample.gamma, file = 'posterior_samples/bru.sample.gamma.Rds')
bru.sample.gamma <- loadRDS('posterior_samples/bru.sample.gamma.Rds')

# take last 10000 posterior samples
mcmc.sample.p <- tail(MCMC.ama, 10000)

# initialize empty matrices - elements would be Lambda(th)  
# different rows correspond to different posterior samples
# different columns correspond to different observations th
expect.mcmc <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))
expect.bru <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))
expect.bru.gamma <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))

# for each posterior sample
for(i in 1:nrow(bru.sample.rep)){
  print(i/nrow(bru.sample.rep))
  # find value of Lambda(th) for each time
  expect.mcmc[i,] <- sapply(dd.ama$time.diff, \(x) mcmc.sample.p[i,1]*(x) + sum(exp(
    log.Lambda_h(th = mcmc.sample.p[i,], ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                 mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )
  
  bru.p <- as.numeric(bru.sample.rep[i,])
  expect.bru[i,] <- sapply(dd.ama$time.diff, \(x) bru.p[1]*(x) + sum(exp(
    log.Lambda_h2(th = bru.p, ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                  mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )  
  bru.p.gamma <- as.numeric(bru.sample.gamma[i,])
  expect.bru.gamma[i,] <- sapply(dd.ama$time.diff, \(x) bru.p.gamma[1]*(x) + sum(exp(
    log.Lambda_h2(th = bru.p.gamma, ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                  mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )
}
saveRDS(expect.mcmc, file = 'posterior_samples/expect.mcmc.Rds')
saveRDS(expect.bru, file = 'posterior_samples/expect.bru.Rds')
saveRDS(expect.bru.gamma, file = 'posterior_samples/expect.bru.gamma.Rds')

expect.mcmc <- loadRDS('posterior_samples/expect.mcmc.Rds')
expect.bru <- loadRDS('posterior_samples/expect.bru.Rds')
expect.bru.gamma <- loadRDS('posterior_samples/expect.bru.gamma.Rds')

# extract median and quantiles value of Lambda(th) 
df.res <- data.frame(days = rep(dd.ama$time.diff,3),
                     Lambda.med = c(apply(expect.mcmc, 2, median), 
                                    apply(expect.bru, 2, median),
                                    apply(expect.bru.gamma, 2, median)),
                     Lambda.low = c( apply(expect.mcmc, 2, \(x) quantile(x, 0.025)), 
                                     apply(expect.bru, 2, \(x) quantile(x, 0.025)),
                                     apply(expect.bru.gamma, 2, \(x) quantile(x, 0.025))),
                     Lambda.up = c(apply(expect.mcmc, 2, \(x) quantile(x, 0.975)), 
                                   apply(expect.bru, 2, \(x) quantile(x, 0.975)),
                                   apply(expect.bru.gamma, 2, \(x) quantile(x, 0.975))),
                     model = rep(c('BayesianETAS', 'Inlabru - rep', 'Inlabru - gamma'), 
                                 each = nrow(dd.ama)),
                     cumfreq = rep(sapply(dd.ama$time.diff, \(x) sum(dd.ama$time.diff < x)), 3))


# initialize matrices for Lambda(th)
xx <- seq(0,1000,length.out = 500)
cs.mcmc <- matrix(NA, ncol = length(xx), nrow = nrow(expect.mcmc))
cs.bru <-  matrix(NA, ncol = length(xx), nrow = nrow(expect.bru))
cs.bru.gamma <-  matrix(NA, ncol = length(xx), nrow = nrow(expect.bru.gamma))

# for each of them calculate cumulative frequencies 
for(i in 1:nrow(cs.mcmc)){
  cs.mcmc[i,] <- sapply(xx, \(x) sum(expect.mcmc[i,] <= x))
  cs.bru[i,] <- sapply(xx, \(x) sum(expect.bru[i,] <= x))
  cs.bru.gamma[i,] <- sapply(xx, \(x) sum(expect.bru.gamma[i,] <= x))
}

# store in a data.frame for plotting
df.cs <- data.frame(Lambda = xx, 
                    cs.med = c( apply(cs.mcmc, 2, median), 
                                apply(cs.bru, 2, median),
                                apply(cs.bru.gamma, 2, median)),
                    cs.low = c( apply(cs.mcmc, 2, \(x) quantile(x, 0.025)), 
                                apply(cs.bru, 2, \(x) quantile(x, 0.025)),
                                apply(cs.bru.gamma, 2, \(x) quantile(x, 0.025))),
                    cs.up = c( apply(cs.mcmc, 2, \(x) quantile(x, 0.975)), 
                               apply(cs.bru, 2, \(x) quantile(x, 0.975)),
                               apply(cs.bru.gamma, 2, \(x) quantile(x, 0.975))),
                    model = rep(c('BayesianETAS', 'Inlabru - rep', 'Inlabru - gamma'), 
                                each = ncol(cs.mcmc))
)

# put together
df.total <- data.frame(x = c( df.res$days, df.cs$Lambda ),
                       cumfreq = c( df.res$cumfreq, rep(NA, nrow(df.cs)) ),
                       xx = c( rep(NA, nrow(df.res)), df.cs$Lambda ),
                       med = c( df.res$Lambda.med, df.cs$cs.med ),
                       low = c( df.res$Lambda.low, df.cs$cs.low ),
                       up = c( df.res$Lambda.up, df.cs$cs.up ),
                       model = c( df.res$model, df.cs$model ),
                       plot = c( rep('days', nrow(df.res)), rep('Lambda', nrow(df.cs)) )
)

# FIGURE 2 - TOP ROW
plot.gof.bru <- ggplot(df.total[df.total$model != 'BayesianETAS',], aes(x = x)) + 
  geom_point(aes(y = cumfreq), size = 1) +
  geom_line(aes(y = xx), linetype = 2) + 
  geom_line(aes(y = med, linetype = model, color = model)) +
  geom_ribbon(aes(ymin = low, ymax = up, 
                  linetype = model, color = model, fill = model), alpha = 0.2) +
  ylab('Cumulative frequencies') + 
  facet_wrap(facets = vars(plot), scales = 'free',
             strip.position = "bottom", 
             labeller = as_labeller(c(days = "days", Lambda = "Lambda") )) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  xlab(NULL) +
  geom_text(data = data.frame(x = c(250,800), y = 250, plot = c('days', 'Lambda'),
                              lab = c('(a)', '(b)')),
            mapping = aes(x, y, label = lab)) + 
  scale_color_manual(values = gg_color_hue(3)[c(3, 1)]) +
  scale_fill_manual(values = gg_color_hue(3)[c(3, 1)]) + 
  #scale_linetype_manual(values = c(1,2)) + 
  theme(legend.title = element_blank())

# FIGURE 2 - BOTTOM ROW
plot.gof.mcmc <- ggplot(df.total[df.total$model != 'Inlabru - gamma',], aes(x = x)) + 
  geom_point(aes(y = cumfreq), size = 1) +
  geom_line(aes(y = xx), linetype = 2) + 
  geom_line(aes(y = med, linetype = model, color = model)) +
  geom_ribbon(aes(ymin = low, ymax = up, 
                  linetype = model, color = model, fill = model), alpha = 0.2) +
  ylab('Cumulative frequencies') + 
  facet_wrap(facets = vars(plot), scales = 'free',
             strip.position = "bottom", 
             labeller = as_labeller(c(days = "days", Lambda = "Lambda") )) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  xlab(NULL) +
  geom_text(data = data.frame(x = c(250,800), y = 250, plot = c('days', 'Lambda'),
                              lab = c('(c)', '(d)')),
            mapping = aes(x, y, label = lab)) + 
  scale_color_manual(values = gg_color_hue(3)[c(2, 1)]) +
  scale_fill_manual(values = gg_color_hue(3)[c(2 ,1)]) + 
  scale_linetype_manual(values = c(1,3)) + 
  theme(legend.title = element_blank())

# FIGURE 2
pdf(file = 'figure2.pdf')
multiplot(plot.gof.bru, plot.gof.mcmc)
dev.off()

#############################
# EXPECTED NUMBER OF EVENTS #
#############################

# ML estimator of magnitude distribution parameter
beta.ml <- 1/mean(data.bru$magnitudes - M0)

# find posterior sample of expected number of events for MCMC
MCMC.N <- vapply(1:nrow(MCMC.ama), \(x) 
                 MCMC.ama[x,1]*(T2) + sum(exp(log.Lambda_h(MCMC.ama[x,], 
                                                           ti = data.bru$ts, 
                                                           mi = data.bru$magnitudes, M0 = M0,
                                                           T1 = T1, T2 = T2) ) ), 0)

# simulated parameters from inlabru model in ETAS scale
par.etas.sim <- data.frame(mu = link.f.be$mu(par.samps[1,]),
                           K = link.f.be$K(par.samps[2,]),
                           alpha = link.f.be$alpha(par.samps[3,]),
                           c = link.f.be$c_(par.samps[4,]), 
                           p = link.f.be$p(par.samps[5,]))
# posterior sample of expected number of events for inlabru
bru.N <- vapply(1:nrow(par.etas.sim), \(x) 
                par.etas.sim[x,1]*(T2 - T1) + sum(exp(log.Lambda_h2(par.etas.sim[x,], ti = data.bru$ts, 
                                                                    mi = data.bru$magnitudes, M0 = M0,
                                                                    T1 = T1, T2 = T2) ) ), 0)

###################
# BRANCHING RATIO #
###################
# simulate 10000 values from magnitude distribution
mm <- rexp(10000, beta.ml) + M0

# branching ratio posterior samples for MCMC
MCMC.branch <- vapply(1:nrow(MCMC.ama), 
                      \(x) mean(MCMC.ama[x,2]*exp(MCMC.ama[x,3]*(mm - M0))), 0)

# branching ratio posterior samples for inlabru
bru.branch <- vapply(1:nrow(par.etas.sim), \(x) 
                     mean(exp(log.Lambda_h2(th = par.etas.sim[x,],
                                            ti = 1,
                                            mi = mm,
                                            M0 = M0, T1 = 0,
                                            T2 = Inf))), 0)



# create data.frame storing the expected number of events and branching ratio posterior samples for plotting
df.comp <- rbind(data.frame(x = bru.N,
                            implementation = 'inlabru - rep',
                            quantity = 'Lambda(T[1], T[2])',
                            true = nrow(data.bru)), 
                 data.frame(x = MCMC.N,
                            implementation = 'bayesianETAS',
                            quantity = 'Lambda(T[1], T[2])',
                            true = nrow(data.bru)),
                 data.frame(x = bru.branch,
                            implementation = 'inlabru - rep',
                            quantity = 'BR',
                            true = NA), 
                 data.frame(x = MCMC.branch,
                            implementation = 'bayesianETAS',
                            quantity = 'BR',
                            true = NA))

# FIGURE 3
pdf('figure3.pdf')
ggplot(df.comp, aes(x = x, color = implementation, linetype = implementation)) + 
  geom_density(size = 1.2) + 
  geom_vline(aes(xintercept = true), color = 'black', size = 1.2, linetype = 3) + 
  facet_wrap(facets = vars(quantity), scales = 'free', 
             labeller = label_parsed) + 
  xlab('value') + 
  theme_classic() + 
  theme(legend.position = 'bottom')
dev.off()



########################################
# RETROSPECTIVE FORECASTING EXPERIMENT #
########################################

# generate posterior sample of ETAS parameters

par.samps <- generate(fit_be, data.frame(), ~ c(th.mu, th.K, th.alpha, th.c, th.p),
                      n.samples = 10000)

# transform them in ETAS scale
par.s.etas <- cbind(mu = link.f.be$mu(par.samps[1,]),
                    K = link.f.be$K(par.samps[2,]),
                    alpha = link.f.be$alpha(par.samps[3,]),
                    c = link.f.be$c_(par.samps[4,]),
                    p = link.f.be$p(par.samps[5,]))
# set number of periods to be forecasts
n.day = 120
# set starting date
f.start <- data.bru$ts[1] + 1e-6
# initialize matrix of forecasting periods
t.lims <- matrix(NA, nrow = n.day, ncol = 3)
# set magnitude above which a forecasting period is splitted
mag.split <- 5.5

# for each period
for(day in 1:n.day){
  # set starting and end time of the forecasting period 
  if(day == 1){
    f.T1 <- f.start
    f.T2 <- f.start + 1
  }
  else{
    f.T1 = f.T2
    f.T2 = f.T1 + 1
  }
  # select data in the forecasting period
  data.T1.T2 <- data.bru[data.bru$ts > f.T1 & data.bru$ts < f.T2, ]
  # check if any event above mag.split
  if(any(data.T1.T2$magnitudes > mag.split)){
    f.T2 = data.T1.T2$ts[which.max(data.T1.T2$magnitudes > mag.split)] + 1e-6
  }
  # store time forecasting interval
  t.lims[day, ] <- c(f.T1,  ( (f.T2 + f.T1)/2 ), f.T2)
  # produce forecast for the period
  single.fore <- cat_forecast(theta.samp = par.s.etas,
                              fore.T1 = f.T1,
                              fore.delta = (f.T2 - f.T1),
                              M0 = M0, beta.p = beta.ml,
                              data.input = data.bru,
                              folder.path = 'fore_cat/',
                              fore.name = paste0('fore.day.', day))
}


# retrieve daily forecasts and get quantiles of the number of events per day
N.quant.sim <- foreach(day = 1:n.day, .combine = rbind) %do% {
  print(day)
  f.cat <- read.table(file = paste0('fore_cat/fore.day.', day ,'.txt'),
                      header = TRUE)
  n.sim <- vapply(1:nrow(par.s.etas), \(x) sum(f.cat$cat.idx == x), 0)
  n.true <- sum(data.bru$ts > t.lims[day, 1] & 
                  data.bru$ts <= t.lims[day,3] )
  data.frame(lower = quantile(n.sim, 0.025),
             median = quantile(n.sim, 0.5),
             upper = quantile(n.sim, 0.975),
             true = n.true)
}

# set rownames
rownames(N.quant.sim) <- NULL
# create data.frame of mids points of each time period for plotting
t.lims.df <- data.frame(t.mid = t.lims[,2])

periods <- 1:n.day
# FIGURE 4 TOP ROW
plot.foreN.natural <- ggplot(fore.N, aes(x = t.lims.df$t.mid, 
                                         y = median)) + 
  geom_line(color = 'red') + 
  geom_ribbon(aes(xmin = t.lims.df$t.mid, 
                  xmax = t.lims.df$t.mid, 
                  ymin = q0.025, 
                  ymax = q0.975),
              alpha = 0.2, color = 'orange', fill = 'orange') + 
  geom_point(aes(x = t.lims.df$t.mid, 
                 y = true)) + #, size = 0.5)) + 
  #scale_y_log10() + 
  xlab('periods') + 
  ylab('N') + 
  annotate('text', x = 85, y = 300, label = '(a)')+
  theme_classic()

# FIGURE 4 BOTTOM ROW
# Remove periods with 0 events
idx.good <- fore.N$true > 0 
plot.foreN.log <- 
  ggplot(fore.N[idx.good,], aes(x = t.lims.df$t.mid[idx.good], 
                                y = log(median))) + 
  geom_line(color = 'red') + 
  geom_ribbon(aes(xmin = t.lims.df$t.mid[idx.good], 
                  xmax = t.lims.df$t.mid[idx.good], 
                  ymin = log(q0.025), 
                  ymax = log(q0.975)),
              alpha = 0.2, color = 'orange', fill = 'orange') + 
  geom_point(aes(x = t.lims.df$t.mid[idx.good], 
                 y = log(true))) + #, size = 0.5)) + 
  #scale_y_log10() + 
  xlab('periods') + 
  ylab('log(N)') + 
  annotate('text', x = 85, y = log(300), label = '(b)')+
  theme_classic()

# FIGURE 4
pdf('figure4.pdf')
multiplot(plot.foreN.natural, plot.foreN.log)
dev.off()

#########################
# SIMULATION EXPERIMENT #
#########################

# retrieve posterior median used for generating the synthetic catalogs
par.median <- fit_be$summary.fixed$`0.5quant`

# transform posterio median in ETAS scale
par.s <- c(link.f.be$mu(par.median[1]),
           link.f.be$K(par.median[2]),
           link.f.be$alpha(par.median[3]),
           link.f.be$c_(par.median[4]),
           link.f.be$p(par.median[5]))

# select the fixed data points in each simulations
data.known <- data.bru[order(data.bru$magnitudes, decreasing = TRUE),][1:3,]
data.known <- data.known[order(data.known$ts),]

# set number of synthetic catalogs
n.samp <- 10000
for(i in 1:n.samp){
  # generate a synthetic catalogue
  samp.etas.list <- sample.ETAS(theta = par.s, beta.p = beta.ml, M0 = 2.99, 
                                T1 = data.bru$ts[1], T2 = T2,
                                Ht = data.known)
  samp.etas <- bind_rows(samp.etas.list)
  # initialize list of samples
  if(i == 1){
    sim.list <- list(samp.etas)
    save(sim.list, file = 'sim.list.Rds')
  }
  # add the synthetic catalog to the list of samples
  else{
    sim.list <- append(sim.list, list(samp.etas))
    # save the list every 1000 synthetic catalogs
    if(i %% 1000 == 0){
      cat('saving \n')
      save(sim.list, file = 'sim.list.Rds')  
    }
  }
  
}
# save the entire list (if the total number of synthetic catalogs is not multiple of 1000)
save(sim.list, file = 'sim.list.Rds')
load('sim.list.Rds')
# retrieve number of events per catalog
N.sim <- vapply(sim.list, \(x) nrow(x), 0)
# have a look at the distrbution of the number of events per synthetic catalog
ggplot() + geom_density(aes(x = N.sim)) + 
  xlim(0,5000) + 
  geom_vline(xintercept = nrow(data.bru))
#summary(N.sim)

# extract catalogs corresponding to the 0.025, 0.5, and 0.975 quantile of the distribution of the number of events
# per catalog
sim.cat.0.025 <- sim.list[[which.min(abs(N.sim - quantile(N.sim, 0.025)))]]
sim.cat.0.5 <- sim.list[[which.min(abs(N.sim - quantile(N.sim, 0.5)))]]
sim.cat.0.975 <- sim.list[[which.min(abs(N.sim - quantile(N.sim, 0.975)))]]
sim.quant.list <- list(sim.cat.0.025, 
                       sim.cat.0.5, 
                       sim.cat.0.975) 
# chosen quantiles
quant.v <- c(0.025, 0.5, 0.975)
# plot names
plot.names <- c('q = 0.025, N = ', 
                'q = 0.5, N = ',
                'q = 0.975, N = ')
# plot of the selected sequences
plot.list <- lapply(1:length(sim.quant.list), \(x) 
                    ggplot(sim.quant.list[[x]], aes(x = ts)) + 
                      geom_histogram(bins = 50) + 
                      labs(title = paste0(plot.names[x], nrow(sim.quant.list[[x]]))) +
                      theme_classic() + 
                      xlab('days') 
)

# have a look
multiplot(plotlist = plot.list) 

# initialize list of posterior samples of ETAS parameters obtained with MCMC
MCMC.post.list <- list()
# initialize vector of computational times
mcmc.time <- c()
for(i in seq_len(length(sim.quant.list))){
  # take the data
  data.sample <- sim.quant.list[[i]]
  # order the observations per time
  data.sample <- data.sample[order(data.sample$ts),]
  # initialize time
  start.t <- Sys.time()
  # obtain posterior samples
  MCMC.post.list[[i]] <- sampleETASposterior(ts = data.sample$ts,
                                             magnitudes = data.sample$magnitudes,
                                             M0 = M0, T = T2, 
                                             sims = 5000, burnin = 5000,
                                             approx = TRUE)
  # store computational time
  mcmc.time[i] <- Sys.time() - start.t
}
# save the list
save(MCMC.post.list, file = 'mcmc/post.list.Rds')
load('mcmc/post.list.Rds')
# set parameters names
par.names <- c('mu', 'K', 'alpha', 'c', 'p')

# store posterior samples as data.frame 
MCMC.post.list <- lapply(1:length(MCMC.post.list), \(x) 
                         data.frame(MCMC.post.list[[x]]) %>% 
                           mutate(Nobs = nrow(sim.quant.list[[x]]),
                                  time = mcmc.time[x]))
# merge them
MCMC.post.df <- bind_rows(MCMC.post.list)
# set column names of the data.frame
colnames(MCMC.post.df) <- c('mu', 'K', 'alpha', 'c', 'p', 'Nobs', 'time')
# store everything in a unique data.frame for plotting
MCMC.post.df <- rbind(data.frame(value = MCMC.post.df$mu,
                                 Nobs = MCMC.post.df$Nobs,
                                 param = 'mu'),
                      data.frame(value = MCMC.post.df$K,
                                 Nobs = MCMC.post.df$Nobs,
                                 param = 'K'),
                      data.frame(value = MCMC.post.df$alpha,
                                 Nobs = MCMC.post.df$Nobs,
                                 param = 'alpha'),
                      data.frame(value = MCMC.post.df$c,
                                 Nobs = MCMC.post.df$Nobs,
                                 param = 'c'),
                      data.frame(value = MCMC.post.df$p,
                                 Nobs = MCMC.post.df$Nobs,
                                 param = 'p'))

# FIGURE 5
pdf('figure5.pdf')
ggplot(MCMC.post.df, aes(x = value, color = factor(Nobs), linetype = factor(Nobs))) + 
  geom_density() +
  scale_color_viridis(discrete = TRUE) + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) + 
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  ylab('value') + 
  xlab('pdf') + 
  labs(color = 'N', linetype = 'N')
dev.off()

# initialize list of inlabru posteriors
bru.post.list <- list()
# initialize vector of computational times
bru.times <- c()
for(i in seq_len(length(sim.quant.list))){
  # take the data
  data.sample <- sim.quant.list[[i]]
  # order the observations per time
  data.sample <- data.sample[order(data.sample$ts),]
  # create event identifier column
  data.sample$idx.p <- seq_len(nrow(data.sample))
  # initialize time
  start.t <- Sys.time()
  # fit model
  fit_samp <- Hawkes.bru2(sample.s = data.sample, # data 
                          M0 = M0, # magnitude of completeness
                          T1 = 0, T2 = T2, # time domain
                          link.functions = link.f.be, # link functions
                          coef.t. = 2, # binning parameter (delta)
                          delta.t. = 0.1, # binning parameter (Delta)
                          N.max. = 3, # binning parameter (n.max)
                          bru.opt = bru.opt.list) # bru options
  # store computational time
  bru.times[i] <- Sys.time() - start.t
  # store posterior of the parameters
  bru.post.list[[i]] <- rbind(data.frame(inla.tmarginal(link.f.be$mu, fit_samp$marginals.fixed$th.mu),
                                     param = 'mu', data = quant.v[i]),
                          data.frame(inla.tmarginal(link.f.be$K, fit_samp$marginals.fixed$th.K),
                                     param = 'K', data = quant.v[i]),
                          data.frame(inla.tmarginal(link.f.be$alpha, fit_samp$marginals.fixed$th.alpha),
                                     param = 'alpha', data = quant.v[i]),
                          data.frame(inla.tmarginal(link.f.be$c_, fit_samp$marginals.fixed$th.c),
                                     param = 'c', data = quant.v[i]),
                          data.frame(inla.tmarginal(link.f.be$p, fit_samp$marginals.fixed$th.p),
                                     param = 'p', data = quant.v[i]))
}
save(bru.post.list, file = 'post.list.synth.Rds')
# merge all posteriors
bru.post.df <- bind_rows(bru.post.list)
bru.post.df$N.obs = as.integer(quantile(N.sim, bru.post.df$data))

# FIGURE 6
pdf('figure6.pdf')
ggplot(bru.post.df, aes(color = factor(N.obs), linetype = factor(N.obs))) + 
  geom_line(mapping = aes(x,y)) + 
  facet_wrap(facets = vars(param),
             scales = 'free', labeller = label_parsed) +
  scale_color_viridis(discrete = TRUE) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  xlab('value') + 
  ylab('pdf') + 
  labs(color = 'N', linetype = 'N')
dev.off()




