delta <- c(0.5,-0.1)
u <- matrix(
  c(1, 1, 1, 1, 1, 1,
    1, 1, 1, 0, 0, 0),
  ncol=2, nrow=6
)
psi <- diag(1, nrow=6, ncol=6)
Sigma <- matrix(0,ncol=6, nrow=6)

Sigma <- psi
for(p in 1:length(delta)){
  Sigma[1:6, 1:6] <- Sigma[1:6, 1:6] + delta[p]*u[,p, drop=F]%*%t(u[,p, drop=F])
}
Sigma

u[,1, drop=F] %*% t(u[,2, drop=F])
