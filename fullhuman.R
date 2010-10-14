source('mixing.R')
source('plot_prevs.R')
source('components.R')

data <- read.csv('data/epidemic.csv')
data_names <- rownames(data)
M <- data$pos_tbg
N <- data$N
gamma <- data$rec_rate
mu <- data$mortality
       
data_names_red <- rownames(data)[M>0]
M_red <- M[M>0]
N_red <- N[M>0]
gamma_red <- gamma[M>0]
mu_red <- mu[M>0]

b <- c(0.000988135,0.007115023,0.002818846,0.03659315,8.613598e-06,0.01530191,0.003335346,0.0006718254,0.06604545,0.004781628,0.0003120496,0.0003309162,0.0006305497,0.001195220,0.005333664,0.003569865,0.001549510,0.00777874,0.004583622,0.001350717,0.004938572,0.004325641,0.0005499413,0.0003719598,0.001507626,8.626667e-05,8.950453e-06,0.001724242,0.004259433,0.0002981199,0.001025018,0.06034878,0.008622904,0.005814368,0.001864966,0.0003477834,0.003275976,0.001114058)
NGM <- diag(N) %*% mixing(b) %*% diag(1/(mu + gamma))
pdf('out/ddr.pdf')
plot_prevs(pars=b, gamma=gamma, mu=mu, mixing_structure=NA, density=TRUE, N=N, M=M, labels=data_names)
pdf()
sink('out/ddr_components.dat')
components(NGM)
sink()

b <- c(0.000644347,0.00778643,0.005342955,0.007327648,0,0.006607401,0,0,0,0.00965748,0,0,0,0,0,0,0,0,0.01624581,0,0.003696261,0,0,0,0.006634353,0,0,0,0.008898947,0,0.06386972,0,0.005385333,0,0,0,0,0)
NGM_red <- diag(N_red) %*% mixing(b) %*% diag(1/(mu_red + gamma_red))
pdf('out/ddr_red.pdf')
plot_prevs(pars=b, gamma=gamma_red, mu=mu_red, mixing_structure=NA, density=TRUE, N=N_red, M=M_red, labels=data_names_red)
pdf()
sink('out/ddr_red_components.dat')
components(NGM_red, domestic=c(2:4), wild=c(5:12))
sink()

b <- c(4.300183e-07,0.0001094344,6.138935e-08,2.059766e-06,1.047227e-05,0.0001252614,1.438046e-05,1.945267e-06)
mixing_structure <- as.matrix(read.csv('data/mixing1.csv', header=F))
NGM <- diag(N) %*% mixing(b, mixing_structure) %*% diag(1/(mu + gamma))
pdf('out/dd.pdf')
plot_prevs(pars=b, gamma=gamma, mu=mu, mixing_structure=mixing_structure, density=TRUE, N=N, M=M, labels=data_names)
dev.off()
sink('out/dd_components.dat')
components(NGM)
sink()

b <- c(3.158982e-09,2.377476e-06,3.011942e-05,6.302696e-06,6.752572e-06,0.0002214773,4.736548e-06,0.0001445645)
mixing_structure <- as.matrix(read.csv('data/mixing1.csv', header=F))
mixing_structure_red <- mixing_structure[M>0,][,M>0]
NGM_red <- diag(N_red) %*% mixing(b, mixing_structure_red) %*% diag(1/(mu_red + gamma_red))
pdf('out/dd_red.pdf')
plot_prevs(pars=b, gamma=gamma_red, mu=mu_red, mixing_structure=mixing_structure_red, density=TRUE, N=N_red, M=M_red, labels=data_names_red)
dev.off()
sink('out/dd_red_components.dat')
components(NGM_red, domestic=c(2:4), wild=c(5:12))
sink()

b <- c(1.074891e-09,1.718487e-05,9.314082e-06,1.300090e-06,2.801727e-07,2.312986e-07,4.735239e-07,1.300974e-05,9.27968e-09,1.627528e-05,2.566608e-06,1.061708e-08,1.631640e-07,5.757445e-08,8.493155e-08,4.887062e-07,4.450815e-07,1.129888e-08,2.124374e-05,4.597474e-07,9.953295e-06,4.193987e-06,1.005658e-06,3.253240e-08,4.096546e-06,1.208186e-06,1.057232e-05,2.601784e-06,2.561739e-05,1.199461e-05,2.817132e-06,1.309010e-06,2.275126e-05,3.426837e-07,8.824426e-07,2.948902e-05,1.079002e-07,7.691737e-08)
mixing_structure <- as.matrix(read.csv('data/mixing2.csv', header=F))
NGM <- diag(N) %*% mixing(b, mixing_structure) %*% diag(1/(mu + gamma))
pdf('out/ddh.pdf')
plot_prevs(pars=b, gamma=gamma, mu=mu, mixing_structure=mixing_structure, density=TRUE, N=N, M=M, labels=data_names)
dev.off()
sink('out/ddh_components.dat')
components(NGM)
sink()

b <- c(1.699889e-09,1.783278e-05,8.097965e-06,1.451894e-07,0,8.067715e-05,0,0,0,1.605634e-05,0,0,0,0,0,0,0,0,2.193679e-05,0,8.477288e-06,0,0,0,6.120427e-06,0,0,0,7.898174e-06,0,8.08661e-06,0,8.34342e-06)
mixing_structure <- as.matrix(read.csv('data/mixing2.csv', header=F))
mixing_structure_red <- mixing_structure[M>0,][,M>0]
NGM_red <- diag(N_red) %*% mixing(b, mixing_structure_red) %*% diag(1/(mu_red + gamma_red))
pdf('out/ddh_red.pdf')
plot_prevs(pars=b, gamma=gamma_red, mu=mu_red, mixing_structure=mixing_structure_red, density=TRUE, N=N_red, M=M_red, labels=data_names_red)
dev.off()
sink('out/ddh_red_components.dat')
components(NGM_red, domestic=c(2:4), wild=c(5:12))
sink()

