#===============================================================================
# Trypanosomiasis multi-species model
#===============================================================================

host_data <- read.csv('~/Code/tryps/data/bipindi.csv', header=T, sep=',')
vector_data <- read.csv('~/Code/tryps/data/bipindi_vector.csv', header=T,
                        sep=',')


n_hosts <- length(host_data$mu)

groups <- list(1, seq(n_hosts-1)+1)

n_groups <- length(groups)
n_vectors <- length(vector_data$mu)

fg <- NULL
for (i in seq(n_groups)) {
  fg[i] <- sum(host_data$f[groups[[i]]])
}

if (n_vectors == 1) {
  vector_data$n <- 1
}

b_hat <- ...
bv <- ...

xi <- 1

lambda <- ...
lambda_v <- ...

parms <- ...

x0 <- ...

# construct initial state vector

x0

for (i in seq(n_groups)) {
  a_host <- c(paste("lambda*(1-I", i, "")...
  for (j in seq(n_vectors)) {
    for (k in seq(n_groups)) {
      
    }
  }
}

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
    nu[n_hosts + 3 * (j - 1) * n_groups + 3 * (i - 1) + 1, column] <- -1
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

for (i in seq(n_groups)) {
  for (k in groups[[i]]) {
    infstring <- paste("b", k, "/n", k, " * (")
    for (j in seq(n_vectors)) {
      infstring <- paste(infstring, "+tau", j, "*f", k, "/fg", i, "*nv", j,
                         "*Iv", ((j-1) * n_groups + k))
    }
    infstring <- paste(infstring, ")")
    a[k] <- infstring
  }
}

for (j in seq(n_vectors)) {
  for (i in seq(n_groups)) {
    infstring <- paste("bv", j, "*tau", j, 
  }
}
