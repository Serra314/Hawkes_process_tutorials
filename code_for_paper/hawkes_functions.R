library(ggplot2)
library(inlabru)
library(INLA)
library(foreach)
#library(bayesianETAS)
library(dplyr)
#library(ggstar)
library(foreach)





# time triggering function - used by bayesianETAS
gt <- function(th, t, ti, mi, M0){
  # set 0 
  output <- rep(0,length(ti))
  before <- ti < t
  # update only ti < t
  if(sum(before) > 0){
    log.out <- log(th[2]) + th[3]*(mi[before] - M0)  + log(th[5] - 1) + (th[5] - 1)*log(th[4]) - th[5]*log(t - ti[before] + th[4])
    output[before] <- exp(log.out)  
  }
  else{
    output
  }
  output
}

# time triggering function - used by Inlabru
gt.2 <- function(th, t, ti, mi, M0){
  output <- rep(0,length(ti))
  t.diff <- t - ti
  neg <- t.diff <= 0
  if(sum(!neg) > 0){
    log.out <- log(th[2]) + th[3]*(mi[!neg] - M0)  - th[5]*log(1 + t.diff[!neg]/th[4])
    output[!neg] <- exp(log.out)  
  }
  else{
    output
  }
  output
}

# conditional intensity - used by bayesianETAS
lambda_ <- function(th, t, ti.v, mi.v, M0){
  if(is.null(ti.v) | all(ti.v > t)){
    th[1]
  }
  th[1] + sum(gt(th, t, ti.v, mi.v, M0))
}

# conditional intensity - used by Inlabru
lambda_2 <- function(th, t, ti.v, mi.v, M0){
  if(is.null(ti.v) | all(ti.v > t)){
    th[1]
  }
  th[1] + sum(gt.2(th, t, ti.v, mi.v, M0))
}

# integrated triggering function - used by bayesianETAS
log.Lambda_h <- function(th, ti, mi, M0, T1, T2){
  th <- as.numeric(th)
  T.low <- pmax(ti, T1)#sapply(ti, \(x) max(T1, x))
  log(th[2]) + th[3]*(mi - M0) + log((th[4]/(T.low - ti + th[4]))^(th[5] - 1) - (th[4]/(T2 - ti + th[4]))^(th[5] - 1))
}


# substitute sapply(ti, \(x) max(T1, x)) with pmax(ti, T1)

# integrated triggering function - used by Inlabru
log.Lambda_h2 <- function(th, ti, mi, M0, T1, T2){
  th <- as.numeric(th)
  T.low <- pmax(T1, ti)#sapply(ti, \(x) max(T1, x))
  
  gamma.l <- (T.low - ti)/th[4]
  gamma.u <- (T2 - ti)/th[4]
  w.l <- (1 + gamma.l)^(1-th[5])  
  w.u <- (1 + gamma.u)^(1-th[5])
  # output
  log(th[2]) + th[3]*(mi - M0) + log(th[4]) - log(th[5] - 1) + log1p(w.l - 1) + log1p(-w.u/w.l) 
  
}

# find breaks point for grid
breaks_exp <- function(tt_, T2_, coef_ = 2, delta_, N_exp_ = 10){
  
  tt_breaks <- tt_ + delta_*((1 + coef_)^(0:N_exp_))
  tt_breaks <- tt_breaks[tt_breaks < T2_]
  if(T2_ - tt_ < delta_){
    return(c(tt_, T2_))
  }
  if(T2_ - tt_breaks[length(tt_breaks)] < delta_){
    tt_breaks[length(tt_breaks)] = T2_
  }
  if(tt_breaks[length(tt_breaks)] < T2_){
    tt_breaks <- c(tt_breaks, T2_)
  }
  return(c(tt_,tt_breaks))
} 


time.grid <- function(data.point, coef.t, delta.t,
                            T2., displaygrid = FALSE, N.exp.){
  
  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # spatial bins
  if(displaygrid){
    Plot_grid(xx = xx., yy = yy., delta_ = delta., n.layer = n.layer., 
              bdy_ =  bdy., min.edge = min.edge.)
  }

  # time bins
  # find bins break points
  t_b <- breaks_exp(tt., T2., coef_ = coef.t, delta_ = delta.t, N_exp_ = N.exp.)
  time.bins <- data.frame(t.start = t_b[-length(t_b)], 
                          t.end = t_b[-1]) %>%
    mutate(t.bin.name = paste0(round(t.start,3),'-',round(t.end,3)))
  if(nrow(time.bins) - 1 == 0){
    time.bins$t.ref_layer = paste0('last-',idx.p)  
  }
  else{
    time.bins$t.ref_layer = c(1:(nrow(time.bins) - 1), paste0('last-',idx.p))
  } 
  cbind(time.bins, data.point, row.names = NULL)
}

It_df <- function(param_, time.df){
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_c <- param_[4]
  param_p <- param_[5]
  T.l <- pmax(tth, T1b) #sapply(1:length(tth), \(x) max(tth[x], T1b[x]))
  fun.l <- (1 + (T.l - tth)/param_c)^(1-param_p)
  fun.u <- (1 + (T2b - tth)/param_c)^(1-param_p)
  ( param_c/ (param_p - 1) )* ( fun.l - fun.u )
}


compute.grid <- function(param., list.input_){
  
  It.vec <- It_df(param_ = param., time.df = list.input_$time.sel)
  
  It.vec[list.input_$Imapping]
}



# function to fit Hawkes process model
Hawkes.bru <- function(sample.s, M0, T1, T2, link.functions = NULL, 
                       coef.t., delta.t., N.max., bru.opt){
  # Expected number of background events
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  # this is the expression of log(Lambda0)
  form.0 <- counts ~ log(link.functions$mu(th.mu)) + log(T2 - T1) 
  # first likelihood
  lik.0 <- inlabru::like(formula = form.0,
                         data = df.0,
                         family = 'poisson',
                         options = list(E = df.0$exposures))
  
  cat('Start creating grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    time.grid(data.point = sample.s[idx,], 
              coef.t = coef.t., 
              delta.t = delta.t., 
              T2. = T2, N.exp. = N.max.
              )
  }
  df.j$counts <- 0
  df.j$exposures <- 1
  df.j$part = 'triggered'
  
  t.names <- unique(df.j$t.ref_layer)
  time.sel <- df.j[vapply(t.names, \(bname) match(TRUE, df.j$t.ref_layer == bname), 0L), , drop = FALSE]
  Imapping <- match(df.j$t.ref_layer, t.names)
  list.input <- list(df_grid = df.j, M0 = M0, Imapping = Imapping, time.sel = time.sel)
  
  cat('Finished creating grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, list.input_, ncore_ = ncore){
    theta_ <- c(0, 
                link.functions$K(th.K[1]), 
                link.functions$alpha(th.alpha[1]), 
                link.functions$c_(th.c[1]), 
                link.functions$p(th.p[1]))
    
    #cat('theta - LogL', theta_, '\n')
    comp. <- compute.grid(param. = theta_, list.input_ = list.input_)
    #print(sum(is.na(comp.list$It)))
    #print(sum(is.infinite(comp.list$It)))
    out <- theta_[3]*(list.input_$df_grid$magnitudes - list.input_$M0) + log(theta_[2] + 1e-100) + log(comp. + 1e-100) 
    out
  }
  
  # creating formula for past events contributions to integrated lambda
  form.j.part <- counts ~ logLambda.h.inla(th.K = th.K, th.alpha = th.alpha,
                                           th.c = th.c, th.p = th.p,
                                           list.input_ = list.input)
  # second for triggered part of the integral
  lik.j.part <- inlabru::like(formula = form.j.part,
                              data = df.j,
                              family = 'poisson',
                              options = list(E = df.j$exposures)) 
  
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, M0){
    
    if(is.null(link.functions)){
      th.p <- c(th.mu[1], th.K[1], th.alpha[1], th.c[1], th.p[1])  
    }
    else{
      th.p <- c(link.functions$mu(th.mu[1]),
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))
    }
    
    
    mean(unlist(mclapply(tt, \(x) {
      th_x <- th < x 
      log(lambda_2(th = th.p, t = x, ti.v = th[th_x], 
                   mi.v = mh[th_x], M0 = M0))
    },
    mc.cores = 5)))
  }
  
  # creating formula for summation part
  form.s.part <- counts ~ loglambda.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                         th.c = th.c, th.p = th.p, tt = sample.s$ts, 
                                         th = sample.s$ts, 
                                         mh = sample.s$magnitudes, M0 = M0)
  
  lik.s.part <- inlabru::like(formula = form.s.part,
                              data = df.s,
                              family = 'poisson',
                              options = list(E = df.s$exposures)) 
  
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0 , prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
  bru(lik.0, lik.j.part, lik.s.part, components = cmp.part, 
      options = bru.opt)
  
  
}






Hawkes.bru2 <- function(sample.s, M0, T1, T2, link.functions = NULL, 
                        coef.t., delta.t., N.max., bru.opt){

  # Expected number of background events
  df.0 <- data.frame(counts = 0, exposures = 1, part = 'background')
  # this is the expression of log(Lambda0)
  
  
  cat('Start creating grid...', '\n')
  time.g.st <- Sys.time()
  df.j <- foreach(idx = 1:nrow(sample.s), .combine = rbind) %do% {
    time.grid(data.point = sample.s[idx,], 
              coef.t = coef.t., 
              delta.t = delta.t., 
              T2. = T2, N.exp. = N.max.
    )
  }
  df.j$counts <- 0
  df.j$exposures <- 1
  df.j$part = 'triggered'
  
  t.names <- unique(df.j$t.ref_layer)
  time.sel <- df.j[vapply(t.names, \(bname) match(TRUE, df.j$t.ref_layer == bname), 0L), , drop = FALSE]
  Imapping <- match(df.j$t.ref_layer, t.names)
  
  cat('Finished creating grid, time ', Sys.time() - time.g.st, '\n')   
  
  logLambda.h.inla <- function(th.K, th.alpha, th.c, th.p, list.input_, ncore_ = ncore){
    theta_ <- c(0, 
                link.functions$K(th.K[1]), 
                link.functions$alpha(th.alpha[1]), 
                link.functions$c_(th.c[1]), 
                link.functions$p(th.p[1]))
    
    #cat('theta - LogL', theta_, '\n')
    comp. <- compute.grid(param. = theta_, list.input_ = list.input_)
    #print(sum(is.na(comp.list$It)))
    #print(sum(is.infinite(comp.list$It)))
    out <- theta_[3]*(list.input_$df_grid$magnitudes - list.input_$M0) + log(theta_[2] + 1e-100) + log(comp. + 1e-100) 
    out
  }
  
  # creating formula for past events contributions to integrated lambda
  # third is for the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = 'SL')
  
  loglambda.inla <- function(th.mu, th.K, th.alpha, th.c, th.p, tt, th, mh, M0){
    
    if(is.null(link.functions)){
      th.p <- c(th.mu[1], th.K[1], th.alpha[1], th.c[1], th.p[1])  
    }
    else{
      th.p <- c(link.functions$mu(th.mu[1]),
                link.functions$K(th.K[1]),
                link.functions$alpha(th.alpha[1]),
                link.functions$c_(th.c[1]),
                link.functions$p(th.p[1]))
    }
    
    
    mean(unlist(mclapply(tt, \(x) {
      th_x <- th < x 
      log(lambda_2(th = th.p, t = x, ti.v = th[th_x], 
                   mi.v = mh[th_x], M0 = M0))
    },
    mc.cores = 5)))
  }
  
  list.input <- list(df_grid = df.j, M0 = M0, Imapping = Imapping, time.sel = time.sel,
                     sample.s = sample.s)
  data.input = bind_rows(df.0, df.s, df.j)
  list.input <- append(list.input, 
                       list(idx.bkg = data.input$part == 'background',
                            idx.trig = data.input$part == 'triggered',
                            idx.sl = data.input$part == 'SL'))
  predictor.fun <- function(th.mu, th.K, th.alpha, th.c, th.p, 
                            list.input, T1, T2, M0){
    
    out <- rep(0, nrow(data.input))
    out[list.input$idx.bkg] <- log(link.functions$mu(th.mu[1])) + log(T2 - T1)
    out[list.input$idx.trig] <- logLambda.h.inla(th.K = th.K, th.alpha = th.alpha,
                                                            th.c = th.c, th.p = th.p,
                                                            list.input_ = list.input)
    out[list.input$idx.sl] <- loglambda.inla(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                                   th.c = th.c, th.p = th.p, 
                                                   tt = list.input$sample.s$ts, 
                                                   th = list.input$sample.s$ts, 
                                                   mh = list.input$sample.s$magnitudes, 
                                                   M0 = M0)
    out
  }
  
  merged.form <- counts ~ predictor.fun(th.mu = th.mu, th.K = th.K, th.alpha = th.alpha,
                                        th.c = th.c, th.p = th.p, list.input = list.input,
                                        T1= T1, T2 = T2, M0 = M0)
  cmp.part <- counts ~ -1 + 
    th.mu(1, model = 'linear', mean.linear = 0 , prec.linear = 1) + 
    th.K(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.alpha(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.c(1, model = 'linear', mean.linear = 0, prec.linear = 1) +
    th.p(1, model = 'linear', mean.linear = 0, prec.linear = 1) 
  
  bru(formula = merged.form, components = cmp.part, data = data.input, family = 'Poisson', 
      options = append(bru.opt, list(E = data.input$exposure)))
  
}


## functions to get the density from a sample and build a copula transformation

q.K.mcmc <- function(pr, lower.tail = TRUE, log.p = FALSE){
  if(lower.tail){
    if(log.p){
      return(quantile(k.bru, exp(pr)))
    }
    else{
      return(quantile(k.bru, pr))
    }
  } 
  if(!lower.tail){
    if(log.p){
      return(quantile(k.bru, 1 - exp(pr)))
    }
    else{
      return(quantile(k.bru, 1 - pr))
    }
  } 
}


mcmc.k.t <- function(x){
  bru_forward_transformation(q.K.mcmc, x)
}

## sampling

## integrated time-triggering function

It <- function(theta, th, T2){
  gamma.u <- (T2 - th)/theta[4]
  ( theta[4]/(theta[5] - 1) )*(1 - (gamma.u + 1)^(1-theta[5]) )
}


## inverse of integrated time-triggering function
Inv.It <- function(theta, omega, th){
  th + theta[4]*( ( 1 - omega * (1/theta[4])*(theta[5] - 1) )^( -1/(theta[5] - 1) ) - 1)
}


## sampling times
sample.time <- function(theta, n.ev, th, T2){
  if(n.ev == 0){
    df <- data.frame(ts = 1, x = 1, y = 1, magnitudes = 1, gen = 0)
    return(df[-1,])
  }
  bound.l <- 0 #It(th.p, th, T)
  bound.u <- It(theta, th, T2)
  unif.s <- runif(n.ev, min = bound.l, max = bound.u)
  t.sample <- Inv.It(theta, unif.s, th)
  t.sample
}


# sampling triggered events
sample.triggered <- function(theta, beta.p, th, n.ev, M0, T1, T2){
  # if the number of events to be placed is zero returns an empty data.frame
  if(n.ev == 0){
    samp.points <- data.frame(x = 1, y = 1, ts = 1, magnitudes = 1)
    samp.points <- samp.points[-1,]
    return(samp.points)
  }
  else{
    
    # sample times
    samp.ts <- sample.time(theta, n.ev, th, T2)
    # sample magnitudes
    samp.mags <- rexp(n.ev, rate = beta.p) + M0

    # build output dataset
    samp.points <- data.frame(ts = samp.ts, magnitudes = samp.mags)
    # return only the ones with time different from NA (the one with NA are outside the interval T1, T2)
    # even though it should not happen given how we built sample.omori
    return(samp.points[!is.na(samp.points$ts),])
  }
}


sample.generation <- function(theta, beta.p, Ht, M0, T1, T2, ncore = 1){
  
  # number of parents
  n.parent <- nrow(Ht)
  # calculate the aftershock rate for each parent in history
  trig.rates <- exp(log.Lambda_h2(th = theta, 
                                  ti = Ht$ts, mi = Ht$magnitudes, 
                                  M0 = M0, T1 = T1, T2 = T2))
  # extract number of aftershock for each parent
  n.ev.v <- sapply(1:n.parent, function(x) rpois(1, trig.rates[x]))
  
  #if(sum(n.ev.v) > 1000){
  #  print('Warning : More than 1000 triggered events')
  #  app <- data.frame(x = 1, y = 1, ts = 1, mags = 1)
  #  app <- app[-1,]
  #  return(app)
  #}
  # if no aftershock has to be generated returns empty data.frame
  if(sum(n.ev.v) == 0){
    app <- data.frame(x = 1, y = 1, ts = 1, magnitudes = 1)
    app <- app[-1,]
    return(app)
  }
  
  # identify parent with number of aftershocks > 0 
  idx.p <- which(n.ev.v > 0)
  
  #print(sample.triggered(theta.v, beta.p, Sigma, Chol.M, n.ev.v[idx.p[1]], Ht[idx.p[1],], T1, T2, M0, bdy, crsobj))
  # sample (in parallel) the aftershocks for each parent 
  sample.list <- mclapply(idx.p, function(idx) 
    sample.triggered(theta = theta, beta.p = beta.p, th = Ht$ts[idx], 
                     n.ev = n.ev.v[idx], M0, T1, T2), mc.cores = ncore)
  
  # bind the data.frame in the list and return
  sample.pts <- bind_rows(sample.list) 
  sample.pts
}



sample.ETAS <- function(theta, beta.p, M0, T1, T2,  
                        Ht = NULL, ncore = 1){
  
  # if the upper extreme greater than lower
  if(T2 < T1){
    stop('Error - right-end of time interval greater than left-end')
  }
  
  n.bkg <- rpois(1, theta[1]*(T2 - T1))

  #cat('Background : ', n.bkg, '\n')
  # if no background events are generated initialize an empty data.frame
  if(n.bkg == 0){
    bkg.df <- data.frame(ts = 1, magnitudes = 1, gen = 0)
    bkg.df <- bkg.df[-1,]
  }
  else{
    # sample bkg events
    # if no bk.field.list element is passed it assumes uniform background rate
    
    # otherwise it samples using the information provided
      
    bkg.df <- data.frame(ts = runif(n.bkg, T1, T2), 
                         magnitudes = rexp(n.bkg, beta.p) + M0, 
                         gen = 1)
  }
  
  
  # if known events are provided
  if(!is.null(Ht)){
    #### TODO : the function has to generate all the points.
    # sample a generation from the known events
    gen.from.past <- sample.generation(theta, beta.p, Ht, M0, T1, T2, ncore)
    # if at least an aftershock is produced
    if(nrow(gen.from.past) > 0){
      # set generation
      gen.from.past$gen = 0
      # Merge first generation and background events
      Gen.list <- list(rbind(gen.from.past, bkg.df))  
    }
    else{
      Gen.list <- list(bkg.df)
    }
    
  } 
  else{
    Gen.list <- list(bkg.df)
  }
  # stop if we have no background events and no events generated from known observations
  if(nrow(Gen.list[[1]]) == 0){
    #print(exp(theta.v[1])*(T2 - T1)*(area(bdy)/1000000))
    #stop('No events generated - increase theta1')
    return(Gen.list)
  }
  
  # initialize flag and gen counter
  flag = TRUE
  gen = 1
  # this goes until the condition inside the loop is met
  while(flag){
    # set parents
    parents <- Gen.list[[gen]]
    #cat('Gen : ', nrow(parents), '\n')
    #print(c(T1,T2))
    #print(range(parents$ts))
    # generate aftershocks
    triggered <- sample.generation(theta, beta.p, parents, 
                                   M0, T1, T2, ncore)
    #print(nrow(triggered))
    # stop the loop if there are no more aftershocks
    if(nrow(triggered) == 0){
      flag = FALSE}
    else{
      # set generations
      triggered$gen = gen + 1
      # store new generation
      Gen.list[[gen + 1]] = triggered
      # update generation counter
      gen = gen + 1
    }
  }
  #append(list(Ht), Gen.list)
  Gen.list
}


## forecasts

cat_forecast <- function(theta.samp, fore.T1, fore.delta, M0, beta.p, data.input, 
                         folder.path, fore.name){
  data.ht <- data.input[data.input$ts < fore.T1,]
  if(nrow(data.ht) == 0){
    data.ht = NULL
  }
  start.i <- 1
  for(i in seq_len(nrow(theta.samp))){
    samp.etas.list <- sample.ETAS(theta = as.numeric(theta.samp[i,]), beta.p = beta.p, M0 = M0, 
                                  T1 = fore.T1, T2 = fore.T1 + fore.delta,
                                  Ht = data.ht)
    samp.etas <- bind_rows(samp.etas.list)
    if(nrow(samp.etas) == 0){
      if(i == start.i){
        start.i <- start.i + 1
      }
      next
    }
    samp.etas$cat.idx <- i
    if(i == start.i){
      write.table(samp.etas, file = paste0(folder.path, fore.name, '.txt'), 
                  col.names = TRUE, row.names = FALSE)
    } else { # else append the new rows to the file
      write.table(samp.etas, file = paste0(folder.path, fore.name, '.txt'), append = TRUE, 
                  col.names = FALSE, row.names = FALSE)
    }
  }
}

## functions to run authomatically

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


input.file.to.list <- function(input_path){
  con <- file(input_path)
  on.exit(close(con))
  par.raw <- readLines(con)
  for(i in 1:length(par.raw)){
    row.i <- par.raw[[i]]
    if(grepl('start.date', row.i)){
      eval(parse(text = row.i))
    } 
    else if(grepl('end.date', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.completeness', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.longitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max.longitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('min.latitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max.latitude', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.path', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.header', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.sep', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.skip', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('catalog.colnames', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_mu', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_K', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_alpha', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_c', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('a_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('b_p', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.mu.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.K.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.alpha.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.c.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('th.p.init', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_iter', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('max_step', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('coef.t', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('DELTA', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('Nmax', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('n.periods', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('period.length', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('start.date.fore', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('magnitude.update', row.i)){
      eval(parse(text = row.i))
    }
    else if(grepl('output.name', row.i)){
      eval(parse(text = row.i))
    }
  }
  # loading catalog
  catalog <- read.table(catalog.path, 
                        header = catalog.header, 
                        sep = catalog.sep, 
                        skip = catalog.skip)
  if(!catalog.header){
    colnames(catalog) <- catalog.colnames
  }
  if(!('time_string' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed time equal to "time_string"')
  }
  if(!('Lon' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed longitudes equal to "time_string"')
  }
  if(!('Lat' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed latitude equal to "time_string"')
  }
  if(!('magnitudes' %in% colnames(catalog))){
    stop('Error in the catalog column names: please set the column name of the observed magnitudes equal to "time_string"')
  }
  #
  # catalog preparation
  start.date <- as.POSIXct(start.date, format = '%Y-%m-%d %H:%M:%OS')
  end.date <- as.POSIXct(end.date)
  catalog <- catalog %>% 
    mutate(time_date = as.POSIXct(gsub(pattern = 'T', 
                                       replacement = ' ', 
                                       x = time_string),
                                       format = '%Y-%m-%d %H:%M:%OS')) %>%
    filter(time_date >= start.date,
           time_date <= end.date,
           Lon >= min.longitude,
           Lon <= max.longitude,
           Lat >= min.latitude,
           Lat <= max.latitude,
           magnitudes >= magnitude.completeness) %>%
    mutate(time_diff = as.numeric(difftime(time_date, start.date, units = 'days')))
  cat('Finish loading & preparing catalog', '\n')
  
  # create data.frame for inlabru
  data.bru <- data.frame(ts = catalog$time_diff,
                         magnitudes = catalog$magnitudes,
                         idx.p = seq_len(nrow(catalog)))
  # set up time interval and magnitude of completeness 
  T1 <- 0
  T2 <- ceiling(as.numeric(difftime(end.date, start.date, units = 'days')))
  M0 <- magnitude.completeness
  
  # priors
  link.f <- list(mu = \(x) gamma.t(x, a_mu, b_mu), 
                 K = \(x) loggaus.t(x, a_K, b_K), 
                 alpha = \(x) unif.t(x, a_alpha, b_alpha), 
                 c_ = \(x) unif.t(x, a_c, b_c), 
                 p = \(x) unif.t(x, a_p, b_p))
  
  # initial value
  th.init <- list(th.mu = th.mu.init,
                  th.K = th.K.init,
                  th.alpha = th.alpha.init,
                  th.c = th.c.init,
                  th.p = th.p.init)
  
  # options for inlabru 
  if(is.null(max_step)){
    bru.opt.list <- list(bru_verbose = 3, # type of visual output 
                         bru_max_iter = max_iter, # maximum number of iterations
                         num.threads = 5,
                         #bru_method = list(max_step = 0.5),
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
  } else {
    bru.opt.list <- list(bru_verbose = 3, # type of visual output 
                         bru_max_iter = max_iter, # maximum number of iterations
                         bru_method = list(max_step = max_step),
                         num.threads = 5,
                         inla.mode = 'experimental', # type of inla algorithm
                         bru_initial = th.init) # parameters initial values
    
  }
    # output list
  list(catalog = catalog,
       catalog.bru = data.bru,
       time.int = c(start.date, end.date),
       T12 = c(T1, T2),
       lat.int = c(min.latitude, max.latitude),
       lon.int = c(min.longitude, max.longitude),
       M0 = M0,
       link.functions = link.f,
       bru.opt.list = bru.opt.list,
       coef.t = coef.t,
       delta.t = DELTA,
       Nmax = Nmax,
       #model.fit = fit_etas,
       n.periods = n.periods,
       period.length = period.length,
       start.date.fore = start.date.fore,
       magnitude.update = magnitude.update,
       output.name = output.name
  )
}


Temporal.ETAS.fit <- function(input.list){
  cat('Start model fitting', '\n')
  fit_etas <- Hawkes.bru2(sample.s = input.list$catalog.bru, # data 
                          M0 = input.list$M0, # magnitude of completeness
                          T1 = input.list$T12[1], 
                          T2 = input.list$T12[2], # time domain
                          link.functions = input.list$link.functions, # link functions
                          coef.t. = input.list$coef.t, # binning parameter (delta)
                          delta.t. = input.list$delta.t, # binning parameter (Delta)
                          N.max. = input.list$Nmax, # binning parameter (n.max)
                          bru.opt = input.list$bru.opt.list) # bru options
  cat('Finish model fitting', '\n')
  fit_etas
}





## get parameters posterior distribution
get_posterior_param <- function(input.list){
  post.mu <- data.frame(inla.tmarginal(input.list$link.functions$mu,
                                       input.list$model.fit$marginals.fixed$th.mu),
                        param = 'mu')
  post.K <- data.frame(inla.tmarginal(input.list$link.functions$K,
                                      input.list$model.fit$marginals.fixed$th.K),
                       param = 'K')
  post.alpha <- data.frame(inla.tmarginal(input.list$link.functions$alpha,
                                          input.list$model.fit$marginals.fixed$th.alpha),
                           param = 'alpha')
  post.c <- data.frame(inla.tmarginal(input.list$link.functions$c_,
                                      input.list$model.fit$marginals.fixed$th.c),
                       param = 'c')
  post.p <- data.frame(inla.tmarginal(input.list$link.functions$p,
                                      input.list$model.fit$marginals.fixed$th.p),
                       param = 'p')
  post.df <- rbind(post.mu, post.K, post.alpha, post.c, post.p)
  post.plot <- ggplot(post.df, aes(x,y)) + 
    geom_line() + 
    facet_wrap(facets = vars(param), scales = 'free', labeller = label_parsed)+
    xlab('param') + 
    ylab('pdf')
  list(post.df = post.df,
       post.plot = post.plot)
}


post_sampling <- function(input.list, n.samp){
  post.samp <- generate(input.list$model.fit, data.frame(),
                        ~ c(input.list$link.functions$mu(th.mu),
                            input.list$link.functions$K(th.K),
                            input.list$link.functions$alpha(th.alpha),
                            input.list$link.functions$c_(th.c),
                            input.list$link.functions$p(th.p)), n.samples = n.samp)
  data.frame(mu = post.samp[1,],
             K = post.samp[2,],
             alpha = post.samp[3,],
             c = post.samp[4,],
             p = post.samp[5,])
}

# number of events
lambda.N <- function(th.mu, th.K, th.alpha, th.c, th.p, T1, T2, M0, Ht,
                     link.functions){
  theta_etas <- c(link.functions$mu(th.mu[1]),
                  link.functions$K(th.K[1]),
                  link.functions$alpha(th.alpha[1]),
                  link.functions$c_(th.c[1]),
                  link.functions$p(th.p[1]))
  
  theta_etas[1]*(T2 - T1) + sum(exp(log.Lambda_h2(th = theta_etas, 
                                                  ti = Ht$ts,
                                                  mi = Ht$magnitudes, 
                                                  M0 = M0,
                                                  T1 = T1, T2 = T2)))
}

get_posterior_N <- function(input.list){
  lambda.N.post <- predict(input.list$model.fit, # model fit 
                           data.frame(), # data (empty because the history of the process is passed to the function below directly)
                           ~ lambda.N(th.mu, th.K, th.alpha, th.c, th.p, 
                                      input.list$T12[1], input.list$T12[2], input.list$M0, 
                                      input.list$catalog.bru,
                                      input.list$link.functions)) # target function
  N.seq <- floor(lambda.N.post$q0.025 - lambda.N.post$q0.025*0.10):ceiling(lambda.N.post$q0.975 + lambda.N.post$q0.975*0.10)
  N.post.df <- predict(input.list$model.fit, data.frame(), 
                       ~ data.frame(N = N.seq, 
                                    pdf = dpois(N.seq, 
                                                lambda.N(th.mu, th.K, 
                                                         th.alpha, th.c, th.p, 
                                                         input.list$T12[1], input.list$T12[2], input.list$M0, 
                                                         input.list$catalog.bru,
                                                         input.list$link.functions))))
  N.post.plot <- ggplot(N.post.df, aes(x = N, y = mean)) + 
    geom_line() + 
    geom_vline(xintercept = nrow(input.list$catalog.bru), linetype = 2) + 
    ylab('pdf')
  
  N.post.plot.shaded <- N.post.plot + 
    geom_ribbon(aes(xmax = N, xmin = N, ymin = q0.025, ymax = q0.975), alpha = 0.2,
                fill = 'blue') 
  list(post.df = N.post.df,
       post.plot = N.post.plot,
       post.plot.shaded = N.post.plot.shaded)
}

produce.forecast <- function(input.list, par.sample, beta.par, folder.path, fore.name){
  n.periods = input.list$n.periods
  f.start <- as.numeric(difftime(input.list$start.date.fore, input.list$time.int[1], units = 'days'))
  t.lims <- matrix(NA, nrow = n.periods, ncol = 2)
  mag.split <- input.list$magnitude.update
  for(period in seq_len(n.periods)){
    if(period/n.periods %/% 0.25 == 0){
      cat('Percentage of simulated periods:', (period/n.periods)*100, '% \n')
    }
    if(period == 1){
      f.T1 <- f.start
      f.T2 <- f.start + input.list$period.length
    }
    else{
      f.T1 = f.T2
      f.T2 = f.T1 + input.list$period.length
    }
    data.T1.T2 <- input.list$catalog.bru[input.list$catalog.bru$ts > f.T1 & input.list$catalog.bru$ts < f.T2, ]
    if(any(data.T1.T2$magnitudes > mag.split)){
      f.T2 = data.T1.T2$ts[which.max(data.T1.T2$magnitudes > mag.split)] + 1e-6
    }
    t.lims[period, ] <- c(f.T1, f.T2)
    single.fore <- cat_forecast(theta.samp = par.sample,
                                fore.T1 = f.T1,
                                fore.delta = (f.T2 - f.T1),
                                M0 = input.list$M0, 
                                beta.p = beta.par,
                                data.input = input.list$catalog.bru,
                                folder.path = folder.path,
                                fore.name = paste0(fore.name, period))
  }
  cat('Forecast finished')
  list(t.lims = t.lims)
}



find.fore.tlim <- function(input.list){
  n.periods = input.list$n.periods
  f.start <- as.numeric(difftime(input.list$start.date.fore, input.list$time.int[1], units = 'days'))
  t.lims <- matrix(NA, nrow = n.periods, ncol = 2)
  mag.split <- input.list$magnitude.update
  for(period in seq_len(n.periods)){
    if(period == 1){
      f.T1 <- f.start
      f.T2 <- f.start + input.list$period.length
    }
    else{
      f.T1 = f.T2
      f.T2 = f.T1 + input.list$period.length
    }
    data.T1.T2 <- input.list$catalog.bru[input.list$catalog.bru$ts > f.T1 & input.list$catalog.bru$ts < f.T2, ]
    if(any(data.T1.T2$magnitudes > mag.split)){
      f.T2 = data.T1.T2$ts[which.max(data.T1.T2$magnitudes > mag.split)] + 1e-6
    }
    t.lims[period, ] <- c(f.T1, f.T2)
  }
  list(t.lims = t.lims)
}





summary.forecast.N <- function(input.list, n.rep, t.lims, folder.path, fore.name){
  N.quant.sim <- foreach(period = 1:input.list$n.periods, .combine = rbind) %do% {
    #print(period)
    f.cat <- read.table(file = paste0(folder.path, fore.name, period ,'.txt'),
                        header = TRUE)
    n.sim <- vapply(1:n.rep, \(x) sum(f.cat$cat.idx == x), 0)
    n.true <- sum(input.list$catalog.bru$ts > t.lims[period, 1] & 
                    input.list$catalog.bru$ts <= t.lims[period,2] )
    data.frame(q0.025 = quantile(n.sim, 0.025),
               median = quantile(n.sim, 0.5),
               q0.975 = quantile(n.sim, 0.975),
               true = n.true)
  }
  rownames(N.quant.sim) <- NULL
  N.quant.sim  
}








