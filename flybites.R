flybites <- function(NGM, humans=c(1), pigs=c(4), sheep=c(2), antelope=c(24))
  {
    proj_labels <- c("humans", "pigs", "sheep", "antelope", "golden cat",
                     "other")

    all <- 1:(nrow(NGM)) 
    others <- setdiff(all, c(humans,pigs,sheep,antelope))
 
    sum_ngm <- sum(NGM)
    
    human_ngm <- sum(NGM[humans,])
    pigs_ngm <- sum(NGM[pigs,])
    sheep_ngm <- sum(NGM[sheep,])
    antelope_ngm <- sum(NGM[antelope,])

    others_ngm <- sum(NGM[others,])

    cat("Humans:", human_ngm/sum_ngm*100, "%\n")
    cat("Pigs:", pigs_ngm/sum_ngm*100, "%\n")
    cat("Sheep:", sheep_ngm/sum_ngm*100, "%\n")
    cat("Antelope:", antelope_ngm/sum_ngm*100, "%\n")
    cat("Others:", others_ngm/sum_ngm*100, "%\n")
  }

