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

NGM_ddr<- diag(N) %*% mixing(b_ddr) %*% diag(1/(mu + gamma))
pdf('out/ddr.pdf')
plot_prevs(pars=b_ddr, gamma=gamma, mu=mu, mixing_structure=NA, density=TRUE, N=N, M=M, labels=data_names)
pdf()
sink('out/ddr_components.dat')
components(NGM_ddr)
sink()
sink('out/ddr_singles.dat')
findres(NGM_ddr, depth = 1)
sink()

NGM_ddr_red <- diag(N_red) %*% mixing(b_ddr_red) %*% diag(1/(mu_red + gamma_red))
pdf('out/ddr_red.pdf')
plot_prevs(pars=b_ddr_red, gamma=gamma_red, mu=mu_red, mixing_structure=NA, density=TRUE, N=N_red, M=M_red, labels=data_names_red)
pdf()
sink('out/ddr_red_components.dat')
components(NGM_ddr_red, domestic=c(2:4), wild=c(5:12))
sink()
sink('out/ddr_red_singles.dat')
findres(NGM_ddr_red, depth = 1)
sink()

mixing_structure1<- as.matrix(read.csv('data/mixing1.csv', header=F))
NGM_dd<- diag(N) %*% mixing(b_dd, mixing_structure1) %*% diag(1/(mu + gamma))
pdf('out/dd.pdf')
plot_prevs(pars=b_dd, gamma=gamma, mu=mu, mixing_structure=mixing_structure1, density=TRUE, N=N, M=M, labels=data_names)
dev.off()
sink('out/dd_components.dat')
components(NGM_dd)
sink()
sink('out/dd_singles.dat')
findres(NGM_dd, depth = 1)
sink()

mixing_structure1_red <- mixing_structure1[M>0,][,M>0]
NGM_dd_red <- diag(N_red) %*% mixing(b_dd_red, mixing_structure1_red) %*% diag(1/(mu_red + gamma_red))
pdf('out/dd_red.pdf')
plot_prevs(pars=b_dd_red, gamma=gamma_red, mu=mu_red, mixing_structure=mixing_structure1_red, density=TRUE, N=N_red, M=M_red, labels=data_names_red)
dev.off()
sink('out/dd_red_components.dat')
components(NGM_dd_red, domestic=c(2:4), wild=c(5:12))
sink()
mixing_structure <- as.matrix(read.csv('data/mixing1.csv', header=F))
sink('out/dd_red_singles.dat')
findres(NGM_dd_red, depth = 1)
sink()

mixing_structure2 <- as.matrix(read.csv('data/mixing2.csv', header=F))
NGM_ddh <- diag(N) %*% mixing(b_ddh, mixing_structure2) %*% diag(1/(mu + gamma))
pdf('out/ddh.pdf')
plot_prevs(pars=b_ddh, gamma=gamma, mu=mu, mixing_structure=mixing_structure2, density=TRUE, N=N, M=M, labels=data_names)
dev.off()
sink('out/ddh_components.dat')
components(NGM_ddh)
sink()
sink('out/ddh_singles.dat')
findres(NGM_ddh, depth = 1)
sink()

mixing_structure2_red <- mixing_structure2[M>0,][,M>0]
NGM_ddh_red <- diag(N_red) %*% mixing(b_ddh_red, mixing_structure2_red) %*% diag(1/(mu_red + gamma_red))
pdf('out/ddh_red.pdf')
plot_prevs(pars=b_ddh_red, gamma=gamma_red, mu=mu_red, mixing_structure=mixing_structure2_red, density=TRUE, N=N_red, M=M_red, labels=data_names_red)
dev.off()
sink('out/ddh_red_components.dat')
components(NGM_ddh_red, domestic=c(2:4), wild=c(5:12))
sink()
sink('out/ddh_red_singles.dat')
findres(NGM_ddh_red, depth = 1)
sink()
