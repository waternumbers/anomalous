# All printed values should be TRUE
rm(list=ls())
load("test_objects.RData")


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
idx <- anomaly_paper_example_1_point_anomalies$location
cost_background <- -2*dnorm(x[idx],0,1,log=TRUE)
cost_pa <- -2*dnorm(x[idx],0,abs(x[idx]),log=TRUE) +  anomaly_paper_example_1_result@beta_tilde ## point anomaly no gamma
all( cost_pa < cost_background )
gamma <- exp(-(1+anomaly_paper_example_1_result@beta_tilde))
cost_pa_pen <- log(2*pi) + log( gamma + x[idx]^2) +1 + anomaly_paper_example_1_result@beta_tilde
all( cost_pa_pen < cost_background )
	
## Example 2
idx <- anomaly_paper_example_2_point_anomalies$location
cost_background <- -2*dnorm(x[idx],0,1,log=TRUE)
cost_pa <- -2*dnorm(x[idx],0,abs(x[idx]),log=TRUE) +  anomaly_paper_example_2_result@beta_tilde
all( cost_pa < cost_background )
gamma <- exp(-(1+anomaly_paper_example_2_result@beta_tilde))
cost_pa_pen <- log(2*pi) + log( gamma + x[idx]^2) +1 + anomaly_paper_example_2_result@beta_tilde
all( cost_pa_pen < cost_background )


## Example 2a
idx <- anomaly_paper_example_2_a_point_anomalies$location
cost_pa <- -2*dnorm(x[idx],0,abs(x[idx]),log=TRUE) +  anomaly_paper_example_2_a_result@beta_tilde
cost_background <- -2*dnorm(x[idx],0,1,log=TRUE)
all( cost_pa < cost_background )
gamma <- exp(-(1+anomaly_paper_example_2_a_result@beta_tilde))
cost_pa_pen <- log(2*pi) + log( gamma + x[idx]^2) +1 + anomaly_paper_example_2_a_result@beta_tilde
all( cost_pa_pen < cost_background )



## ###############################
## test_that("Example 4 from vignettes",

## Part 1
idx <- anomaly_paper_example_4_1_point_anomalies$location
cost_pa <- -2*dnorm(x[idx],0,abs(x[idx]),log=TRUE) +  anomaly_paper_example_4_1_result@beta_tilde
cost_background <- -2*dnorm(x[idx],0,1,log=TRUE)
all( cost_pa < cost_background )
gamma <- exp(-(1+anomaly_paper_example_4_1_result@beta_tilde))
cost_pa_pen <- log(2*pi) + log( gamma + x[idx]^2) +1 + anomaly_paper_example_4_1_result@beta_tilde
all( cost_pa_pen < cost_background )

## Part 2
idx <- anomaly_paper_example_4_2_point_anomalies$location
cost_pa <- -2*dnorm(x[idx],0,abs(x[idx]),log=TRUE) +  anomaly_paper_example_4_2_result@beta_tilde
cost_background <- -2*dnorm(x[idx],0,1,log=TRUE)
all( cost_pa < cost_background )
gamma <- exp(-(1+anomaly_paper_example_4_2_result@beta_tilde))
cost_pa_pen <- log(2*pi) + log( gamma + x[idx]^2) +1 + anomaly_paper_example_4_2_result@beta_tilde
all( cost_pa_pen < cost_background )

