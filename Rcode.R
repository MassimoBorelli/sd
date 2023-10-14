
############################################################
# import dataset and divide numbers by 1000 
# in order to avoid trouble with point/comma decimal symbol

url = "https://raw.githubusercontent.com/MassimoBorelli/sd/main/wanXiEta.csv"
wanXiEta = read.csv(url, header = TRUE, sep = ";")
attach(wanXiEta); xi_n = xi_n/1000; eta_n = eta_n/1000
############################################################


############################################################
# formulas (1) and (2), page 3
BlomXi = 2*qnorm((n-0.375)/(n+0.25))
BlomEta = 2*qnorm((0.75*n-0.125)/(n+0.25))
############################################################


############################################################
# residuals delta(n) and epsilon(n), page 4
residualityXi = DELTA = xi_n - BlomXi
residualityEta = EPSILON = eta_n - BlomEta
############################################################


############################################################
# Section 2 Figure 1 page 4
par(mfrow = c(2,2))
plot(n, xi_n, pch = 20, col = "blue", xlab = "", ylab = "")
points(n, BlomXi, col = "red", pch = 22)
legend( x="bottomright", 
        legend=c(expression(paste(xi, "(n)")),expression(paste("2",Phi^-1, "(.)"))),        
        col=c("blue","red"), 
        pch=c(20,22) )
plot(n, residualityXi, "h",
     main = expression(paste( delta, "(n) = ", xi, "(n)", " - 2" ,Phi^-1 , "(.)" ) )
     , xlab = "", ylab = "")
#par(mfrow = c(1,1))
#par(mfrow = c(1,2))
plot(n, eta_n, pch = 20, col = "blue", xlab = "", ylab = "")
points(n, BlomEta, col = "red", pch = 22)
legend( x="bottomright", 
        legend=c(expression(paste(eta, "(n)")),expression(paste("2",Phi^-1, "(.)"))),        
        col=c("blue","red"), 
        pch=c(20,22) )
plot(n, residualityEta, "h",
     main = expression(paste( epsilon, "(n) = ", eta, "(n)", " - 2" ,Phi^-1 , "(.)" ) )
     , xlab = "", ylab = "")
par(mfrow = c(1,1))
############################################################


############################################################
# summary lm, page 5
Y = n/log(EPSILON)
modelEpsilon = lm (Y ~ n)
coef(modelEpsilon)
summary(modelEpsilon)
############################################################


############################################################
# unpublished diagnostic
par(mfrow = c(2,2))
plot(modelEpsilon)
par(mfrow = c(1,1))
# it appears a not negligible curvature in residuals
# it could be better to introduce a curvature term
modelEpsilonSquared = lm (Y ~ n + I(n^2))
summary(modelEpsilonSquared)
par(mfrow = c(2,2))
plot(modelEpsilonSquared)
par(mfrow = c(1,1))
# therefore it would be better to adopt: 
ourImprovedNewEpsilon = exp(n/(-2.6065 -0.2609 * n + 0.0006 * n^2))
# instead of 'ourNewEpsilon'
############################################################


############################################################
# page 5, our new epsilon(n) and our new eta(n)
ourNewEpsilon = exp(n/(-2.8822 - 0.2308 * n ))
ourNewEta = BlomEta + ourNewEpsilon
############################################################


############################################################
# Section 2.1 Figure 2 page 5
par(mfrow = c(1,2))
plot(n, Y, pch = 24, col = "magenta", xlab = "n", ylab = "",
     main = expression(paste( "Y(n) = n / log(", epsilon, "(n))")))
abline(modelEpsilon, pch = 8, col = "purple",lwd = 2, lty = 4)
legend(x="bottomleft", legend="Y(n)", col="magenta", pch=24)
plot(n, eta_n, pch = 20, col = "blue", xlab = "", ylab = "")
points(n, BlomEta, col = "red", pch = 22)
points(n, ourNewEta, col = "orange", pch = 5)
legend( x="bottomright", 
        legend=c(expression(paste(eta, "(n)")),expression(paste("2",Phi^-1, "(.)")), 
                 expression(paste("2",Phi^-1, "(.) + ", epsilon, "(n)"))),        
        col=c("blue","red","orange"), 
        pch=c(20,22,5) )
par(mfrow = c(1,1))
############################################################


############################################################
# Section 2.1 page 5, limit of our new epsilon(n) when n tends to \infty
exp(1/coef(modelEpsilon)[[2]])
############################################################


############################################################
# Section 2.1 page 6, sup and inf
max(abs(eta_n - BlomEta)) < 0.58  # 0.5795677
min(abs(eta_n - BlomEta)) > 0.02  # 0.0296914
max(abs(eta_n - ourNewEta)) < 0.03  # 0.02972438
min(abs(eta_n - ourNewEta)) > 0.0001  # 0.0001589739
############################################################


############################################################
# Section 2.2 page 6, sup and inf
max(abs(xi_n - BlomXi)) < 0.051  # 0.0509116
min(abs(xi_n - BlomXi)) > 0.0002  # 0.000276987
############################################################


############################################################
# Section 2.2 page 6, derivative with central difference quotient
derivative = (EPSILON[3:49] - EPSILON[1:47])/2
plot(2:48,1/derivative, "l") # oscillations in right graph tail
############################################################


############################################################
# summary lm, page 6
logn = log(n)
(modelDelta = lm (residualityXi ~ logn))
coef(modelDelta)
summary(modelDelta)
############################################################


############################################################
# page 6, our new delta(n) 
# page 7, formula (8) our new xi(n)
ourNewDelta = -0.0626 + 0.0197 * logn
ourNewEta = BlomXi + ourNewDelta
############################################################






