tryps <- read.csv(file="epidemic.csv",head=TRUE,sep=",")
l2=(tryps$pos_tbg+tryps$pos_tbng)/(tryps$N-(tryps$pos_tbg+tryps$pos_tbng)) * (tryps$mortality + tryps$rec_rate)
b=l2/sqrt(sum(tryps$N*l2*l2/(l2+tryps$mortality+tryps$rec_rate)));
R = (b*tryps$N %o% b) %*% diag(1/(tryps$mortality+tryps$rec_rate))
R0 = max(Re(eigen(R)$values))
cat ("R0 = ", R0, "\n")
unit=diag(1,38,38)
for (i in c(1:38)) {
  proj=matrix(0,38,38)
  proj[i,i]=1
  U=(unit-proj) %*% R
  V=proj %*% R
  usr=max(Re(eigen(U)$values))
  vsr=max(Re(eigen(V)$values))
  cat ("i = ", i, ", usr = ", usr, ", vsr = ", vsr, "\n")
}
