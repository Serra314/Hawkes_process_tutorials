# function to obtain mean and log of a lognormal distribution
mean.ln <- function(mean.log, sd.log){
  exp(mean.log + (1/2)*(sd.log^2) )
}
sd.ln <- function(mean.log, sd.log){
  sqrt( exp(2*mean.log + sd.log^2) * (exp(sd.log^2) - 1) )
}

par.names <- c('mu', 'K', 'alpha', 'c', 'p')
# function to extract posterior of the parameters
extract.post.df <- function(fit_, link.f, p.name = par.names){
  foreach(i = 1:length(link.f), .combine = rbind) %do% {
    distr = inla.tmarginal(link.f[[i]], fit_$marginals.fixed[[i]])
    
    data.frame(x = distr[,1], y = distr[,2], param = p.name[i])
  }
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
