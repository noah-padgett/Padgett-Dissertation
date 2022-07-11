
library(lavaan)

mod <- '

f =~ lam1*y1 + lam2*y2 + lam3*y3 + lam4*y4 + lam5*y5

y1 ~~ psi1*y1
y2 ~~ psi2*y2
y3 ~~ psi3*y3
y4 ~~ psi4*y4
y5 ~~ psi5*y5

omega := ((lam1 + lam2 + lam3 + lam4 + lam5)**2)/((lam1 + lam2 + lam3 + lam4 + lam5)**2 + psi1 + psi2 + psi3 + psi4 + psi5)
'

dat <- as.data.frame(t(sim.data$ystar))
colnames(dat) <- paste0("y",1:5)
fit <- sem(mod, data=dat, std.lv=T)
summary(fit, standardized=T)


mod2 <- '

f =~ lam1*y1 + lam2*y2 + lam3*y3 + lam4*y4 + lam5*y5

psi1 := 1 - lam1**2
psi2 := 1 - lam2**2
psi3 := 1 - lam3**2
psi4 := 1 - lam4**2
psi5 := 1 - lam5**2

omega := ((lam1 + lam2 + lam3 + lam4 + lam5)**2)/((lam1 + lam2 + lam3 + lam4 + lam5)**2 + psi1 + psi2 + psi3 + psi4 + psi5)
'
fit2 <- sem(mod2, data=sim.data$y, ordered=T, std.lv=T)
summary(fit2, standardized=T)


fit2 <- sem(mod2, data=sim.data$Ysampled, ordered=T, std.lv=T)
summary(fit2, standardized=T)
