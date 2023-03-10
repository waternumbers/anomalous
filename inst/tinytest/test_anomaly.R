## this code is adapted from that of tha anomaly package.
## Major changes are that:
## - ported from testhat to tinytest
## - multivariate tests are not used...
## - we try to second guess the defaults of the capa package
## - we have to allow for the differnt package output formats

load("anomaly_test_objects.RData")

## conversion functions
fc <- function(tmp,type){
    out <- collective_change(tmp)
    out <- out[out$segment==type,c("start","end")]
    return(out) 
}

fp <- function(tmp,type){
    out <- point_change(tmp)
    out <- out[out$segment==type,"index"]
    if(is.null(out)){ out <- integer(0) }
    return(out)
}


## ########################################################
## test_that("Example 1, 2, and 2a from vignettes",
set.seed(0)
x <- rnorm(5000)
x[401:500] <- rnorm(100, 4, 1)
x[1601:1800] <- rnorm(200, 0, 0.01)
x[3201:3500] <- rnorm(300, 0, 10)
x[c(1000, 2000, 3000, 4000)] <- rnorm(4, 0, 100)
x <- (x - median(x))/mad(x)

## Example 1
mu <- rep(0,length(x))
sigma <- rep(1,length(x))
beta <- 4*log(length(x))
betaP <- 3*log(length(x))

expect_silent({ res <- capa(x,mu,sigma,gaussMeanVar,beta,betaP,min_length=10) })
expect_equivalent( anomaly_paper_example_1_collective_anomalies[,c("start","end")],fc(res,"gaussMeanVar") )
expect_equal( anomaly_paper_example_1_point_anomalies$location, fp(res,"gaussPoint") )
  
	
## Example 2
beta <- 3*log(length(x))
expect_silent({ res <- capa(x,mu,sigma,gaussMean,beta,betaP,min_length=10) })
expect_equivalent( anomaly_paper_example_2_collective_anomalies[,c("start","end")],fc(res,"gaussMean") )
expect_equal( anomaly_paper_example_2_point_anomalies$location, fp(res,"gaussPoint") ) # this fails

## double checking
##all(fp(res,"gaussPoint") %in% anomaly_paper_example_2_point_anomalies$location) ## TRUE
##tmp <- sapply(anomaly_paper_example_2_point_anomalies$location, function(ii){ update(gaussPoint(betaP,ii),x[ii],mu[ii],sigma[ii])@cost - update(gaussFixed(0,ii),x[ii],mu[ii],sigma[ii])@cost }) ## should be negative for point anomalies
##anomaly_paper_example_2_point_anomalies$location[tmp>0]
##setdiff( anomaly_paper_example_2_point_anomalies$location,fp(res,"gaussPoint"))



## Example 2a
expect_silent({ res <- capa(1 + 2 * x, mu, sigma,gaussMean,beta,betaP,min_length=10) })
expect_equivalent( anomaly_paper_example_2_a_collective_anomalies$start,fc(res,"gaussMean")$start ) ## fail
expect_equal( anomaly_paper_example_2_a_point_anomalies$location, fp(res,"gaussPoint") ) ### fail


## ###############################
## test_that("Example 4 from vignettes",
data("machinetemp")
x <- (machinetemp$temperature - median(machinetemp$temperature)) / mad(machinetemp$temperature)

mu <- rep(0,length(x))
sigma <- rep(1,length(x))
beta <- 4*log(length(x))
betaP <- 3*log(length(x))

## Part 1
expect_silent({ res <- capa(x,mu,sigma,gaussMean,beta,betaP,min_length=2) })
expect_equivalent( anomaly_paper_example_4_1_collective_anomalies[,c("start","end")],fc(res,"gaussMean") ) ## fail
expect_equal( anomaly_paper_example_4_1_point_anomalies$location,fp(res,"gaussPoint") )

## Part 2
library(robustbase)
n <- length(x)
x.lagged <- matrix(c(x[1:(n - 1)], x[2:n]), n - 1, 2)
phi <- robustbase::covMcd(x.lagged, cor = TRUE)$cor[1,2]
inflated_penalty <- 3 * (1 + phi) / (1 - phi) * log(n)

expect_silent({ res <- capa(x, mu, sigma, gaussMean, beta = inflated_penalty, betaP = inflated_penalty) })
expect_equal( anomaly_paper_example_4_2_collective_anomalies[, c("start","end")], fc(res) )
expect_equal( anomaly_paper_example_4_2_point_anomalies$location, fp(res) )



## test_that("Example 5 from vignettes",
## {
## 	data(simulated)

## 	# Part 2
## 	{
## 		res <- capa(sim.data, type = "mean", min_seg_len = 2)

## 		expect_equal(anomaly_paper_example_5_1_result,res,ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_1_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_1_point_anomalies,point_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_1_point_plot,plot(res),ignore_attr=TRUE)
##     	}

## 	# Part 2
## 	{
## 		beta <- 2 * log(ncol(sim.data):1)
## 		beta[1] <- beta[1] + 3 * log(nrow(sim.data))
## 		res<-capa(sim.data, type= "mean", min_seg_len = 2,beta = beta)

## 		expect_equal(anomaly_paper_example_5_2_result,res)
## 		expect_equal(anomaly_paper_example_5_2_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_2_point_anomalies,point_anomalies(res),ignore_attr=TRUE)
## 		expect_equal(anomaly_paper_example_5_2_point_plot,plot(res),ignore_attr=TRUE)

## 	}
## })


## test_that("Example 6 from vignettes",
## {
## 	set.seed(0)
## 	x1 <- rnorm(500)
## 	x2 <- rnorm(500)
## 	x3 <- rnorm(500)
## 	x4 <- rnorm(500)
## 	x1[151:200] <- x1[151:200] + 2
## 	x2[171:200] <- x2[171:200] + 2
## 	x3[161:190] <- x3[161:190] - 3
## 	x1[351:390] <- x1[371:390] + 2
## 	x3[351:400] <- x3[351:400] - 3
## 	x4[371:400] <- x4[371:400] + 2
## 	x4[451] <- x4[451] * max(1, abs(1 / x4[451])) * 6
## 	x4[100] <- x4[100] * max(1, abs(1 / x4[100])) * 6
## 	x2[050] <- x2[050] * max(1, abs(1 / x2[050])) * 6
## 	x1 <- (x1 - median(x1))/mad(x1)
## 	x2 <- (x2 - median(x2))/mad(x2)
## 	x3 <- (x3 - median(x3))/mad(x3)
## 	x4 <- (x4 - median(x4))/mad(x4)
## 	x <- cbind(x1, x2, x3, x4)
## 	res <- capa(x, max_lag = 20, type = "mean")

## 	expect_equal(anomaly_paper_example_6_result,res,ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_6_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_6_point_anomalies,point_anomalies(res),ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_6_point_plot,plot(res),ignore_attr=TRUE)
## })

## test_that("Example 7 from vignettes",
## {
## 	data(simulated)
## 	res <- pass(sim.data, max_seg_len = 20, alpha = 3)
	
## 	expect_equal(anomaly_paper_example_7_result,res,ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_7_collective_anomalies,collective_anomalies(res),ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_7_point_plot,plot(res),ignore_attr=TRUE)
## })

## test_that("Example 8 from vignettes",
## {
## 	data(simulated)
## 	bard.res <- bard(sim.data)
## 	sampler.res <- sampler(bard.res, gamma = 1/3, num_draws = 1000)

## 	expect_equal(anomaly_paper_example_8_bard_result,bard.res,ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_8_sampler_result,sampler.res,ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_8_collective_anomalies,collective_anomalies(sampler.res),ignore_attr=TRUE)
## 	expect_equal(anomaly_paper_example_8_point_plot,plot(sampler.res),ignore_attr=TRUE)
## })







