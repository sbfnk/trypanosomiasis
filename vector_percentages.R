source('mixing.R')
source('plot_prevs.R')
source('components.R')
source('flybites.R')

data <- read.csv('data/epidemic.csv')
data_names <- rownames(data)
M <- data$pos_tbg + data$pos_tbng
N <- data$N
gamma <- data$rec_rate
mu <- data$mortality
       
data_names_red <- rownames(data)[M>0]
M_red <- M[M>0]
N_red <- N[M>0]
gamma_red <- gamma[M>0]
mu_red <- mu[M>0]

b <- c(0.001357134,0.005492852,0.004488954,0.02236555,0.0075093,6.159601e-05,4.553958e-05,0.0003571679,0.0003146496,0.002079010,0.0001083927,0.0002977478,0.001189424,0.001924247,0.001844523,0.01061361,0.001110220,9.137268e-05,0.005420015,4.847543e-05,0.002821891,1.877051e-05,5.250066e-07,0.003070415,0.003472695,0.004951839,1.189856e-05,0.004059073,0.002832169,6.61021e-06,0.002165979,0.0004639785,0.001867161,0.0009626925,0.0003860022,5.796255e-06,0.0003765397,2.064087e-05)
NGM <- diag(N) %*% mixing(b) %*% diag(1/(mu + gamma))
sink('out/ddr_flybites.dat')
flybites(NGM)
sink()

b <- c(7.903457e-05,0.01377720,0.004046161,0.07116131,0.01040221,0.002677451,0,0.04035635,0,0.001035528,0,0.0004056032,0.0005139624,0.001469887,0.003485261,0.02375452,0,0,0.01279922,0,0.001523343,0,0,0.01265212,0.002697245,0,0,0.02161392,0.002498156,0,0.007036866,0,0.001413993,0.002120933,0.001304112,0,0,0)
NGM_red <- diag(N_red) %*% mixing(b) %*% diag(1/(mu_red + gamma_red))
sink('out/ddr_red_flybites.dat')
flybites(NGM_red, antelope=c(16))
sink()

b <- c(2.167887e-08,4.820159e-08,5.011142e-11,1.289012e-06,4.486956e-07,0.0003227288,3.037489e-05,1.787614e-07)
mixing_structure <- as.matrix(read.csv('data/mixing1.csv', header=F))
NGM <- diag(N) %*% mixing(b, mixing_structure) %*% diag(1/(mu + gamma))
sink('out/dd_flybites.dat')
flybites(NGM)
sink()

b <- c(1.055808e-09,9.06602e-10,1.203618e-06,1.415980e-06,4.146252e-11,0.0003398417,1.738282e-09,0.0001174262)
mixing_structure <- as.matrix(read.csv('data/mixing1.csv', header=F))
mixing_structure_red <- mixing_structure[M>0,][,M>0]
NGM_red <- diag(N_red) %*% mixing(b, mixing_structure_red) %*% diag(1/(mu_red + gamma_red))
sink('out/dd_red_flybites.dat')
flybites(NGM_red, antelope=c(16))
sink()

b <- c(9.82214e-09,1.666936e-05,1.156880e-05,6.367033e-05,2.939591e-05,6.413902e-06,9.022743e-06,2.409961e-05,1.538586e-06,6.684299e-06,1.199471e-06,4.389808e-08,2.652156e-06,3.59685e-06,1.440070e-05,1.135789e-05,3.065452e-05,4.926048e-06,1.343586e-05,4.104977e-07,9.694865e-06,4.627153e-05,2.516006e-07,2.423866e-06,8.490936e-06,1.152008e-05,6.259129e-08,2.687414e-06,1.334728e-05,1.405244e-07,6.052675e-06,1.363829e-07,5.669149e-06,2.439442e-06,4.021935e-06,5.919051e-07,6.577544e-07,3.559915e-08)
mixing_structure <- as.matrix(read.csv('data/mixing2.csv', header=F))
NGM <- diag(N) %*% mixing(b, mixing_structure) %*% diag(1/(mu + gamma))
sink('out/ddh_flybites.dat')
flybites(NGM)
sink()

b <- c(1.383193e-08,1.868983e-05,1.020734e-05,5.938432e-05,2.019028e-05,1.098337e-05,0,6.63383e-06,0,5.288411e-06,0,6.160857e-06,3.990223e-06,5.181757e-06,1.142659e-05,3.640035e-06,0,0,1.617036e-05,0,1.166462e-05,0,0,1.359035e-06,1.192522e-05,0,0,2.222662e-05,7.36959e-06,0,9.46223e-06,0,7.810611e-06,1.007531e-05,3.353963e-06)
mixing_structure <- as.matrix(read.csv('data/mixing2.csv', header=F))
mixing_structure_red <- mixing_structure[M>0,][,M>0]
NGM_red <- diag(N_red) %*% mixing(b, mixing_structure_red) %*% diag(1/(mu_red + gamma_red))
sink('out/ddh_red_flybites.dat')
flybites(NGM_red, antelope=c(16))
sink()
