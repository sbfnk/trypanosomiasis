#===============================================================================
# Trypanosomiasis multi-species model
#===============================================================================

host_data <- read.csv('~/Code/tryps/data/bipindi.csv', header=T, sep=',')
vector_data <- read.csv('~/Code/tryps/data/bipindi_vector.csv', header=T,
                        sep=',')


n_hosts <- length(host_data$mu)

#groups <- list(1, seq(n_hosts-1)+1)
groups <- list(seq(n_hosts))

n_groups <- length(groups)
n_vectors <- length(vector_data$mu)

fg <- NULL
for (i in seq(n_groups)) {
  fg[i] <- sum(host_data$f[groups[[i]]])
}

if (n_vectors == 1) {
  vector_data$n <- 1
}

b_hat <-
  c(0.105686, 7.72722, 3.89085, 0.026204, 1.39793, 1.66512, 1.17079, 1.75079,
    8.65013, 1.02053, 1.77992, 0.533177) 
bv <-
  c(0.308736)

xi <- 0

parms <- c(b=b_hat, n=host_data$n, f=host_data$f, mu=host_data$mu,
           gamma=host_data$gamma, bv=bv, nv=vector_data$n, muv=vector_data$mu,
           tau=vector_data$tau, alpha=vector_data$alpha, xi=xi, fg=fg,
           N=round(host_data$n*3540/7), Nv=10000)

# construct initial state vector

x0 <- round(c(I=host_data$M/host_data$N*host_data$n*3540/7,
              Cv=vector_data$M/vector_data$N/4*1000,
              Iv=vector_data$M/vector_data$N/4*1000,
              Gv=vector_data$M/vector_data$N/4*1000))

# construct state-change matrix

nu <- matrix(0, nrow=n_hosts + 3 * n_vectors * n_groups,
             ncol = 2 * n_hosts + (6 + 3 * (n_groups - 1)) * n_vectors * n_groups)


column <- 1
  
for (i in seq(n_hosts)) {
  # susceptible host infected by vector
  nu[i, column] <- +1
  column <- column + 1
  # infected host loses infectiousness
  nu[i, column] <- -1
  column <- column + 1
}

  
for (i in seq(n_groups)) {
  for (j in seq(n_vectors)) {
    # teneral vector infected by host
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 1, column] <- +1
    column <- column + 1
    # infection maturing in incubating vector
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 1, column] <- -1
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 2, column] <- +1
    column <- column + 1
    # incubating vector death
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 1, column] <- -1
    column <- column + 1
    # infected vector death
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 2, column] <- -1
    column <- column + 1
    # non-infectious bite on host
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 3, column] <- +1
    column <- column + 1
    # non-infected vector death
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 3, column] <- -1
    column <- column + 1
    # host switches
    for (k in seq(n_groups)) {
      if (k != i) {
        nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 1, column] <- -1
        nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (k - 1) + 1, column] <- +1 
        column <- column + 1
        nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 2, column] <- -1
        nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (k - 1) + 2, column] <- +1 
        column <- column + 1
        nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 3, column] <- -1
        nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (k - 1) + 3, column] <- +1 
        column <- column + 1
      }
    }    
  }
}

# construct propensity vector

a <- NULL

column <- 1
for (i in seq(n_groups)) {
  if (n_groups > 1) {
    istring <- i
  } else {
    istring <- ""
  }
  for (k in groups[[i]]) {
    if (n_hosts > 1) {
      kstring <- k
    } else {
      kstring <- ""
    }
    infstring <- paste("b", kstring, " / n", kstring, " * (", sep="")
    for (j in seq(n_vectors)) {
      if (n_vectors > 1) {
        jstring <- j
        vstring <- ((j-1) * n_groups + k)
      } else {
        jstring <- ""
        vstring <- ""
      }
      infstring <- paste(infstring, " + tau", jstring, " * f", kstring,
                         " / fg", istring, " * nv", jstring, 
                         " * Iv", vstring, " / Nv", vstring, sep="")
    }
    infstring <- paste(infstring, ") * (N", kstring, " - I", kstring, ")", sep="")
    a[column] <- infstring
    column <- column + 1
    a[column] <- paste("(mu", kstring, " + gamma", kstring, ") * I", kstring, sep="")
    column <- column + 1
  }
}

for (j in seq(n_vectors)) {
  if (n_vectors > 1) {
    jstring <- j
  } else {
    jstring <- ""
  }
  for (i in seq(n_groups)) {
    if (n_groups > 1) {
      istring <- i
    } else {
      istring <- ""
    }
    lambdastring <- paste("bv", jstring, " * tau", jstring, " / fg", istring,
                          " * (", sep="") 
    for (k in groups[[i]]) {
      if (n_hosts > 1) {
        kstring <- k
      } else {
        kstring <- ""
      }
      if (n_vectors > 1) {
        vstring <- ((j-1) * n_groups + i)
      } else {
        vstring <- ""
      }
      lambdastring <- paste(lambdastring, " + f", kstring, " * I", kstring,
                            " / N", kstring, sep="")
    }
    lambdastring <- paste(lambdastring, ")", sep="")
    tstring <- paste("(Nv", vstring, " - Cv", vstring, " - Iv", vstring,
                     " - Gv", vstring,")", sep="") 
    
    a[column] <- paste(lambdastring, " * ", tstring, sep="")
    column <- column + 1
    a[column] <- paste("alpha", jstring, " * Cv", vstring, sep="")
    column <- column + 1
    a[column] <- paste("muv", jstring, " * Cv", vstring, sep="")
    column <- column + 1
    a[column] <- paste("muv", jstring, " * Iv", vstring, sep="")
    column <- column + 1
    a[column] <- paste("(tau", jstring, " - ", lambdastring, ") *", tstring, sep="")
    column <- column + 1
    a[column] <- paste("muv", jstring, " * Gv", vstring, sep="")
    column <- column + 1
    for (l in seq(n_groups)) {
      if (l != i) {
        a[column] <- paste("xi", jstring, " * fg ", l, " * Cv", vstring, sep="")
        column <- column + 1
        a[column] <- paste("xi", jstring, " * fg ", l, " * Iv", vstring, sep="")
        column <- column + 1
        a[column] <- paste("xi", jstring, " * fg ", l, " * Gv", vstring, sep="")
        column <- column + 1
      }
    }
  }
}
