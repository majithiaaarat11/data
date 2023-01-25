install.packages("POT")
install.packages("tseries")
install.packages("tseries")
install.packages("PerformanceAnalytics")
install.packages("ugarch")
install.packages("mixtools")
install.packages("ReIns")
install.packages("mclust")
install.packages('evir')
install.packages("ugarch")
install.packages("readxl")
library(mclust)
library(evir)
library(ReIns)
library(mixtools)
library(readxl)
library(tseries)
library(mixtools)
library(readxl)
library(tseries)
library(PerformanceAnalytics)
library(ugarch)
data <- read_excel("533960-605045.xlsx")

returns <- data[1:2]

#1. Variance?covariance approach based on the two series

mean_stock_1 <- mean(returns$Stock1)
variance_stock_1 <- var(returns$Stock1)
mean_stock_2 <- mean(returns$Stock2)
variance_stock_2 <- var(returns$Stock2)
covariance_matrix <- cov(returns)
covariance_stocks <- covariance_matrix[1,2]

w1 <- 0.6 #weight of the Stock 1
w2 <- 0.4 #weight of the Stock 2

portfolio <- w1*returns[,1] + w2*returns[,2] #Is this correct?
colnames(portfolio) <- "stocks"
mean_portfolio <- mean(portfolio$stocks)
variance_portfolio <- var(portfolio$stocks)
VaR99_1 <- mean_portfolio + sqrt(variance_portfolio)*qnorm(.99) #10.8263
#Results
mean_stock_1
variance_stock_1
mean_stock_2
variance_stock_2
covariance_stocks
mean_portfolio
variance_portfolio
VaR99_1
#2. Historical Simulation based on the portfolio loss returns

VaR99_2 <- quantile(portfolio$stocks, 0.99) # 11.9628
VaR99_2

#3. A normal mixture model mixing two normal distributions

jarque.bera.test(portfolio)
JB_statistic <- 445.81
p_value <- 2.2e-16


res <- normalmixEM(portfolio$stocks, k=2)
n = 100000
mixture_m <- res$mu
mixture_s <- res$sigma
mixture_p <- res$lambda
mixture_p

plot(res,2)
simGM=rnormmix(n, res$lambda, res$mu,
               res$sigma)
portfolio_sorted = sort(simGM)

VaR99_3 <- quantile(portfolio_sorted, 0.99) #12.3901
VaR99_3
#4. An EVT approach

threshold = hill(portfolio$stocks, p =0.99)
threshold_df = data.frame(unlist(threshold[["OrderStatistic"]]), unlist(threshold[["Alpha"]]))

k = 150
portfolio_sorted = sort(portfolio$stocks)
n = length(portfolio$stocks)
alpha = rep(0, 2000)

for (j in 1:k)
{ 
  for(i in 1:k)
  {
    alpha[j] = alpha[j] + (log(portfolio_sorted[n-i+1]) - log(portfolio_sorted[n-j]))
  }
  alpha[j] = (alpha[j]/j)
}

alpha = alpha^(-1)
plot(alpha)
alpha = alpha[150]
alpha
portfolio_sorted[n-k]
VaR99_4 <- portfolio_sorted[n-k]*(k/(n*0.01))^(1/alpha)
VaR99_4 #12.3267

#5. Fitting a GARCH model with normal innovations (QMLE) and then
#predict the next day VaR using the dynamic historical simulation method
spec = ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = TRUE));
estimate <- ugarchfit(spec, portfolio$stocks);
estimate 
# mu = -0.03007, omega = 0.128307, alpha1= 0.059622, beta1 = 0.934137
FTSE.dy<-estimate@fit$residuals
predicted_vol <- sqrt(0.128307+0.059622*(-0.500)^2+0.934137*(FTSE.dy[2000])^2)
VaR99_5_error <- quantile(FTSE.dy, p=0.99,method = "historical")
VaR99_5 <- -(-0.030007)+predicted_vol*VaR99_5_error




#6. Fitting the two series to a bivariate distribution with the marginals fol-
#low (different) Student?t distributions and the dependence is modeled by a reverse
#Clayton copula. Here the reverse means that (1 -F1(X1),1 -F2(X2)) follows
#the Clayton copula. In that way, (X1,X2) possesses upper tail dependence.
library("copula")
library("fitdistrplus")

mle1 <- fitdistr(returns$Stock1, "t")
mle1

mle2 <- fitdistr(returns$Stock2, "t")
mle2

tau<- cor(returns[,1], returns[,2], method = "kendall")
theta <- 2/((1-tau)-2)

hill <- 1999;
lamda <- matrix(0,hill,1);
sortedStock1 <- sort(returns$Stock1);
sortedStock2 <- sort(returns$Stock2);
n = 2000

for (k in 1:hill) {
    for (t in 1:n){
       if(returns$Stock1[t]>sortedStock1[n-k] && returns$Stock2[t]>sortedStock2[n-k]) {
       lamda[k] = lamda[k] + 1
}
    }
  lamda[k] = lamda[k]/k;
}

lamda
plot(lamda)

estimated_tail_dep <-  lamda[250]
calculated_tail_dep = 2^(-1/theta) 

m4 = fitCopula(claytonCopula(), data = cbind(returns$Stock1, returns$Stock2), method = 'itau')
summary(m4)

library(copula)
library(survival)
nb <- length(returns$Stock1)
rata <- cbind(returns$Stock1,returns$Stock2)
tauEst <- cor(returns$Stock1, returns$Stock2, method="kendall")
thetaClayton <- 2/(1-tauEst)

start=c(2,2,2)

myClayton <- mvdc(copula=claytonCopula(param=1, dim=2), c("t", "t"),
                  paramMargins=list(list(df=mle1$estimate['df']), list(df= mle1$estimate['df'])))
fitClayton <- fitMvdc(rata, myClayton,start=start)
fitClayton


corBLfit <- fitClayton@mvdc@copula@parameters
meanBfit <- fitClayton@estimate[1]
nb <- length(returns$Stock1)
sdBfit <- sqrt(nb/(nb-1))*fitClayton@estimate[2] #unbiased MLE of the sd
rateLfit <- fitClayton@estimate[3]

mcClayton <- mvdc(copula=claytonCopula(param=corBLfit, dim=2),
                  margins=c("t", "t"),
                  paramMargins=list(list(mle1$estimate['df']), list(df= mle2$estimate['df'] )))

nsim <- 100000
blSim <- rMvdc(nsim, mcClayton)
simPortfolio = w1*blSim[,1]+w2*blSim[,2];
VaR99_6 = quantile(simPortfolio,p= 0.99)




