get_ngm <- function(data, vector, params, ignorezeroes = TRUE, gambiense = TRUE,
                    nongambiense = FALSE, convert_area= FALSE,
                    vector_prevalence = FALSE)
{  
  library('BB');
  library('mnormt');
  
  source('joinfactors.R');
  source('betaffoiv.R');

  # store variables
  rN <- data$N
  vN <- vector$N

  rM <- rep(0, length(rN))
  vM <- rep(0, length(vN))
  if (!is.null(gambiense)) {
    rM <- rM + data$pos_tbg
    vM <- vM + vector$pos_tbg
  }

  if (!is.null(nongambiense)) {
    rM <- rM + data$pos_tbng
    vM <- vM + vector$pos_tbng
  }

  rprev <- rM/rN

  if (!is.null(vector_prevalence)) {
    vprev <- vM/vN
  } else {
    vprev <- NULL
  }

  rmu <- data$mortality
  vmu <- vector$mortality
  rgamma <- data$rec_rate
  vgamma <- vector$rec_rate
  rlambda <- rprev/(1-rprev)*(rmu + rgamma)
  vlambda <- vprev/(1-vprev)*(vmu + vgamma)
  theta <- data$theta
  rabundance <- data$abundance
  vdensity <- vector$density

  biting_rate <- vector$biting_rate

  if (!is.null(convert_area)) {
    area_convert <- params$area_convert
  } else {
    area_convert <- NULL
  }

  if (!is.null(ignorezeroes)) {
    rN <- rN[rM>0]
    rmu <- rmu[rM>0]
    rgamma <- rgamma[rM>0]
    vN <- vN[vM>0]
    vmu <- vmu[vM>0]
    vgamma <- vgamma[vM>0]
  }

  if (!is.null(ignorezeroes)) {
    rM <- rM[rM>0]
    vM <- vM[vM>0]
  }

  factor <- joinfactors(theta, biting_rate, area_convert)
  res <- betaffoiv(rlambda, vdensity, rabundance, factor, vprev, rprev, vmu)
  beta <- res$par

  NGM <- matrix(0,length(rgamma)+length(vgamma), length(rgamma)+length(vgamma))

  NGM[1:length(vgamma),(length(vgamma)+1):(length(vgamma)+length(rgamma))] <-
    1/vmu %*% t(beta) * factor * vdensity/rabundance

  for (i in 1:length(vgamma)) {
    NGM[(length(vgamma)+1):(length(vgamma)+length(rgamma)), i] <-
      beta/(rgamma+rmu)
  }

  calcprev <- sum(beta * factor * rprev) / sum(vmu + beta * factor * rprev)

  R0 <- sqrt(sum(NGM[2:(length(rgamma)+1), 1] * NGM[1, 2:(length(rgamma)+1)]))

  species <- NGM[i+1, 1]*NGM[1, i+1]

  return (list(NGM = NGM, R0 = R0, species = species, calcprev = calcprev,
              beta = beta))
}
