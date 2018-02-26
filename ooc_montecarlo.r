#
# ooc_montecarlo.r
#
rm(list=ls(all=TRUE))
#
library(MASS)
library(MCMCpack)
library(hitandrun)
library(ooc)
library(statar)
library(Matrix)
library(reshape)
library(ggplot2)
library(viridis)
library(foreach)
library(doParallel)
#
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)
#
set.seed(1985)
#
setwd("c:/")
#
#
# Utility functions plot
#
x <- seq(0, 2.0, by=0.01)
linearutility <- 1 - abs(x-0)
normalutility <- exp(-2 * ((x-0)^2))
quadraticutility <- 1 - 2*(x-0)^2
#
#
#pdf("Dropbox/OC_publicopinion/images/montecarlo_utilityfunctions.pdf", height=6, width=6)
#
plot(c(0,2), c(-1,1), type="n", bty="n", cex.lab=1.1, xlab="Distance from Ideal Point", ylab="Utility", main="Linear, Normal, and Quadratic Utility Functions")
lines(x,linearutility, lwd=2, lty=1)
lines(x,normalutility, lwd=2, lty=2)
lines(x,quadraticutility, lwd=2, lty=3)
#
dev.off()
#
#
# DEFINE MONTE CARLO FUNCTION
#
montecarlo.oc <- function(n=1200, q=20, ndim=2, utility.probs=c(0.33,0.33,0.33),
	missing=0.1, error.respondents=c(0.1,1), error.issues=c(2,1)){
#
# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%    BEGIN    %%%%%%%
# %%%%% Monte Carlo %%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%
#
heteroskedastic.respondents <- runif(n, error.respondents[1], error.respondents[2])
heteroskedastic.issues <- rgamma(q, error.issues[1], error.issues[2])
#
correlations <- runif(n, -0.1, 0.7)
knowledge <- xtile(correlations, 3)
#
# 1.) Generate respondent ideal points
#
mu <- rep(0,ndim)
Sigma <- list()
idealpoints <- matrix(NA, nrow=n, ncol=ndim)
for (i in 1:n){
	Sigma[[i]] <- matrix(1, nrow=ndim, ncol=ndim)
	Sigma[[i]][lower.tri(Sigma[[i]])] <- runif(sum(lower.tri(Sigma[[i]])), correlations[i]-0.2, correlations[i]+0.2)
	Sigma[[i]][upper.tri(Sigma[[i]])] <- t(Sigma[[i]])[upper.tri(t(Sigma[[i]]))]
	Sigma[[i]] <- as.matrix(nearPD(Sigma[[i]])$mat)
	idealpoints[i,] <- mvrnorm(1, mu=mu, Sigma=Sigma[[i]])
}
#
idealpoints[idealpoints > 2] <- 2
idealpoints[idealpoints < -2] <- -2
#
# 2.) Generate normal vectors
#
normalvectors <- hypersphere.sample(ndim, q)
#
for (j in 1:q){
if (normalvectors[j,1] < 0) normalvectors[j,] <- -1 * normalvectors[j,]
}
#
# 3.) Project respondents on normal vectors
#
respondent.projections <- idealpoints %*% t(normalvectors)
#
# 4.) Generate outcome locations
#
outcome.locations <- apply(respondent.projections, 2, function(x){sort(runif(5,min(x),max(x)))})
#
# 5.) Define utility functions
#
linear.utility <- function(idealpoint, choices){
    tmp <- 1 - (3*abs(idealpoint - choices))
    return(tmp)}
#
normal.utility <- function(idealpoint, choices){
    tmp <- 5 * exp(-1 * ((idealpoint - choices)^2))
    return(tmp)}
#
quad.utility <- function(idealpoint, choices){
    tmp <- 1 - 2 * (idealpoint - choices)^2
    return(tmp)}
#
# 6.) Calculate utility and response probabilities
#
utility.fun <- sample(c("linear","normal","quadratic"), n, replace=TRUE, prob=utility.probs)
#
systematic.utility <- total.utility <- probmat <- list()
for (j in 1:q){
systematic.utility[[j]] <- matrix(NA, nrow=n, ncol=5)
total.utility[[j]] <- matrix(NA, nrow=n, ncol=5)
probmat[[j]] <- matrix(NA, nrow=n, ncol=5)
}
#
for (i in 1:n){
for (j in 1:q){
#
	if(utility.fun[i]=="linear"){
	systematic.utility[[j]][i,] <- linear.utility(respondent.projections[i,j], outcome.locations[,j])
	}
#
	if(utility.fun[i]=="normal"){
	systematic.utility[[j]][i,] <- normal.utility(respondent.projections[i,j], outcome.locations[,j])
	}
#
	if(utility.fun[i]=="quadratic"){
	systematic.utility[[j]][i,] <- quad.utility(respondent.projections[i,j], outcome.locations[,j])
	}
#
total.utility[[j]][i,] <- systematic.utility[[j]][i,] / exp(heteroskedastic.respondents[i] * heteroskedastic.issues[j])
probmat[[j]][i,] <- pnorm(total.utility[[j]][i,]) / sum(pnorm(total.utility[[j]][i,]))
}}
#
for (j in 1:q){
probmat[[j]][is.na(probmat[[j]])] <- 1
}
#
# 7.) Generate perfect and error (probabilistic) responses
#
simulated.responses <- perfect.responses <- matrix(NA, nrow=n, ncol=q)
#
for (i in 1:n){
for (j in 1:q){
simulated.responses[i,j] <- sample(1:5, 1, prob=probmat[[j]][i,])
# Note: perfect voting could of course also be simulated using the maximum of total.utility
perfect.responses[i,j] <- which.max(systematic.utility[[j]][i,])
}}
#
# 8.) Insert missing data at random
#
if(missing > 0){
miss <- expand.grid(n=1:n, q=1:q)
missing.selected <- as.matrix(miss)[sample(1:nrow(miss), (n*q*missing), replace=FALSE),]
for (i in 1:nrow(missing.selected)){
simulated.responses[missing.selected[i,1], missing.selected[i,2]] <- NA
perfect.responses[missing.selected[i,1], missing.selected[i,2]] <- NA
}
}
#
# 9.) Organize output
#
correctvotes <- sum(diag(table(simulated.responses,perfect.responses)))
totalvotes <- sum(table(simulated.responses,perfect.responses))
#
return(list(simulated.responses = simulated.responses,
	    perfect.responses = perfect.responses,
	    idealpoints = idealpoints,
	    normalvectors = normalvectors,
	    heteroskedastic.respondents = heteroskedastic.respondents,
	    heteroskedastic.issues = heteroskedastic.issues,
	    correlations = correlations,
	    knowledge = knowledge,
	    error = (1 - (correctvotes / totalvotes))
	))
#
}
#
# %%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%     END    %%%%%%%
# %%%%% Monte Carlo %%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@  RUN ORDERED  @@@@@@@@@@@@@
# @@@@@@@@  OPTIMAL CLASSIFICATION  @@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#
# utility.probs = c(linear, normal, quadratic)
# error.respondents = normal(mu, sigma)
# error.issues = gamma(shape, rate)
#
# Simulations randomize:
#
#		1.) Respondent utility functions
#		2.) Missingness
#		3.) Error rate
#
#
ntrials <- 100
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
# Basic setup:
#sim <- montecarlo.oc(n=1500, q=40, ndim=2, utility.probs=runif(3,0,1), missing=0.1, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
#res <- ooc(sim$perfect.responses, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
#
#
#
# I. VARY AMOUNT OF MISSING DATA
#
fit.2dim.missing10 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=2, utility.probs=runif(3,0,1), missing=0.1, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.2dim.missing10, file="Dropbox/OC_publicopinion/analysis/fit.2dim.missing10.rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
fit.2dim.missing25 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=2, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.2dim.missing25, file="Dropbox/OC_publicopinion/analysis/fit.2dim.missing25.rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
fit.2dim.missing40 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=2, utility.probs=runif(3,0,1), missing=0.4, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.2dim.missing40, file="Dropbox/OC_publicopinion/analysis/fit.2dim.missing40.rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
fit.3dim.missing10 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=3, utility.probs=runif(3,0,1), missing=0.1, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=3, min=10, lop=0.0001, polarity=rep(1,3), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- 1
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.3dim.missing10, file="Dropbox/OC_publicopinion/analysis/fit.3dim.missing10.rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
fit.3dim.missing25 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=3, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=3, min=10, lop=0.0001, polarity=rep(1,3), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- 1
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.3dim.missing25, file="Dropbox/OC_publicopinion/analysis/fit.3dim.missing25.rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
fit.3dim.missing40 <- foreach(i=101:(100+ntrials), .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=3, utility.probs=runif(3,0,1), missing=0.4, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=3, min=10, lop=0.0001, polarity=rep(1,3), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- 1
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.3dim.missing40, file="Dropbox/OC_publicopinion/analysis/fit.3dim.missing40.rda")
#
# !!!!!!!!!!!!!!!!!!!!!!!!!
#
#
#
#
# II. INFORMAL TEST OF CONSISTENCY:
#    VARY NUMBER OF ISSUES (25, 50, 100)
# 
fit.2dim.issues25 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=25, ndim=2, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.2dim.issues25, file="Dropbox/OC_publicopinion/analysis/fit.2dim.issues25.rda")
#
# 
fit.2dim.issues50 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=50, ndim=2, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.2dim.issues50, file="Dropbox/OC_publicopinion/analysis/fit.2dim.issues50.rda")
#
# 
fit.2dim.issues100 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=100, ndim=2, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.2dim.issues100, file="Dropbox/OC_publicopinion/analysis/fit.2dim.issues100.rda")
#
#
fit.3dim.issues25 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=25, ndim=3, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=3, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.3dim.issues25, file="Dropbox/OC_publicopinion/analysis/fit.3dim.issues25.rda")
#
# 
fit.3dim.issues50 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=50, ndim=3, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=3, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.3dim.issues50, file="Dropbox/OC_publicopinion/analysis/fit.3dim.issues50.rda")
#
# 
fit.3dim.issues100 <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=100, ndim=3, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=3, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
normvecs <- res$issues.unique[,grepl("normVector", colnames(res$issues.unique)) & colnames(res$issues.unique)!="normVectorAngle2D"]
#
# Perform Procrustes rotation
rotation <- procrustes(x, sim$idealpoints)$R
rotated.idealpoints <- x %*% rotation
rotated.normalvectors <- normvecs %*% rotation
rotated.normalvectors[(rotated.normalvectors[,1] < 0),] <- -1 * rotated.normalvectors[(rotated.normalvectors[,1] < 0),]
#
truetheta <- apply(sim$normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
truetheta[sim$normalvectors[,2] < 0] <- -1 * truetheta[sim$normalvectors[,2] < 0]
#
theta <- apply(rotated.normalvectors, 1, function(x){ acos_d(sum(x * c(1,0))/(sqrt(sum(x * x)) * sqrt(sum(c(1,0) * c(1,0)))))})
theta[rotated.normalvectors[,2] < 0] <- -1 * theta[rotated.normalvectors[,2] < 0]
#
respondents.correct <- res$respondents[,"percent.correctScale"]
issues.PRE <- res$issues.unique[,"PREScale"]
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
cornormvec <- cor(as.vector(dist(rotated.normalvectors)), as.vector(dist(sim$normalvectors)))
cortheta <- cor(truetheta, theta)
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corconstraintresp <- cor(sim$correlations, respondents.correct)
corhetissues <- cor(sim$heteroskedastic.issues, issues.PRE)
#
error <- sim$error
c(corx, cornormvec, cortheta, corhetresp, corconstraintresp, corhetissues, error)
}
#
save(fit.3dim.issues100, file="Dropbox/OC_publicopinion/analysis/fit.3dim.issues100.rda")
#
#
#
#
# III. COMPARE FITS OF MODERATES AND EXTREMISTS
#    
fit.2dim.moderatesextremists <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=2, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
respondents.correct <- res$respondents[,"percent.correctScale"]
#
moderate.1dim <- rotated.idealpoints[,1] > quantile(sim$idealpoints[,1])["25%"] &
rotated.idealpoints[,1] < quantile(sim$idealpoints[,1])["75%"]
moderate.2dim <- rotated.idealpoints[,2] > quantile(sim$idealpoints[,2])["25%"] &
rotated.idealpoints[,2] < quantile(rotated.idealpoints[,2])["75%"]
moderate <- moderate.1dim & moderate.2dim
#
extreme.1dim <- rotated.idealpoints[,1] < quantile(sim$idealpoints[,1])["25%"] |
rotated.idealpoints[,1] > quantile(sim$idealpoints[,1])["75%"]
extreme.2dim <- rotated.idealpoints[,2] < quantile(sim$idealpoints[,2])["25%"] |
rotated.idealpoints[,2] > quantile(sim$idealpoints[,2])["75%"]
extreme <- extreme.1dim & extreme.2dim
#
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
corx.moderate <- cor(as.vector(dist(x[moderate,])), as.vector(dist(sim$idealpoints[moderate,])))
corx.extreme <- cor(as.vector(dist(x[extreme,])), as.vector(dist(sim$idealpoints[extreme,])))
#
#
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corhetresp.moderate <- cor(sim$heteroskedastic.respondents[moderate], respondents.correct[moderate])
corhetresp.extreme <- cor(sim$heteroskedastic.respondents[extreme], respondents.correct[extreme])
#
corconstraintresp <- cor(sim$correlations, respondents.correct)
corconstraintresp.moderate <- cor(sim$correlations[moderate], respondents.correct[moderate])
corconstraintresp.extreme <- cor(sim$correlations[extreme], respondents.correct[extreme])
#
error <- sim$error
c(corx, corx.moderate, corx.extreme, corhetresp, corhetresp.moderate, corhetresp.extreme, corconstraintresp, corconstraintresp.moderate, corconstraintresp.extreme)
}
#
save(fit.2dim.moderatesextremists, file="Dropbox/OC_publicopinion/analysis/fit.2dim.moderatesextremists.rda")
#
#
fit.3dim.moderatesextremists <- foreach(i=1:ntrials, .combine='rbind', .packages=c("ooc", "statar", "Matrix")) %dopar% {
set.seed(i)
sim <- montecarlo.oc(n=1500, q=40, ndim=3, utility.probs=runif(3,0,1), missing=0.25, error.respondents=sort(runif(2,0,0.5)), error.issues=c(runif(1,0,3), 0.5))
issuescales <- sim$simulated.responses
res <- ooc(issuescales, dims=2, min=10, lop=0.0001, polarity=rep(1,2), iter=25, nv.method="svm.reg", cost=1)
x <- res$respondents[,grepl("coord", colnames(res$respondents))]
respondents.correct <- res$respondents[,"percent.correctScale"]
#
moderate.1dim <- rotated.idealpoints[,1] > quantile(sim$idealpoints[,1])["25%"] &
rotated.idealpoints[,1] < quantile(sim$idealpoints[,1])["75%"]
moderate.2dim <- rotated.idealpoints[,2] > quantile(sim$idealpoints[,2])["25%"] &
rotated.idealpoints[,2] < quantile(rotated.idealpoints[,2])["75%"]
moderate.3dim <- rotated.idealpoints[,3] > quantile(sim$idealpoints[,3])["25%"] &
rotated.idealpoints[,3] < quantile(rotated.idealpoints[,3])["75%"]
moderate <- moderate.1dim & moderate.2dim & moderate.3dim
#
extreme.1dim <- rotated.idealpoints[,1] < quantile(sim$idealpoints[,1])["25%"] |
rotated.idealpoints[,1] > quantile(sim$idealpoints[,1])["75%"]
extreme.2dim <- rotated.idealpoints[,2] < quantile(sim$idealpoints[,2])["25%"] |
rotated.idealpoints[,2] > quantile(sim$idealpoints[,2])["75%"]
extreme.3dim <- rotated.idealpoints[,3] < quantile(sim$idealpoints[,3])["25%"] |
rotated.idealpoints[,3] > quantile(sim$idealpoints[,3])["75%"]
extreme <- extreme.1dim & extreme.2dim & extreme.3dim
#
#
corx <- cor(as.vector(dist(x)), as.vector(dist(sim$idealpoints)))
corx.moderate <- cor(as.vector(dist(x[moderate,])), as.vector(dist(sim$idealpoints[moderate,])))
corx.extreme <- cor(as.vector(dist(x[extreme,])), as.vector(dist(sim$idealpoints[extreme,])))
#
#
corhetresp <- cor(sim$heteroskedastic.respondents, respondents.correct)
corhetresp.moderate <- cor(sim$heteroskedastic.respondents[moderate], respondents.correct[moderate])
corhetresp.extreme <- cor(sim$heteroskedastic.respondents[extreme], respondents.correct[extreme])
#
corconstraintresp <- cor(sim$correlations, respondents.correct)
corconstraintresp.moderate <- cor(sim$correlations[moderate], respondents.correct[moderate])
corconstraintresp.extreme <- cor(sim$correlations[extreme], respondents.correct[extreme])
#
error <- sim$error
c(corx, corx.moderate, corx.extreme, corhetresp, corhetresp.moderate, corhetresp.extreme, corconstraintresp, corconstraintresp.moderate, corconstraintresp.extreme)
}
#
save(fit.3dim.moderatesextremists, file="Dropbox/OC_publicopinion/analysis/fit.3dim.moderatesextremists.rda")
#
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
stopCluster(cl)
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $$$$$$$$$$  PLOTS  $$$$$$$$$$$$$$$$
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#
load("Dropbox/OC_publicopinion/analysis/fit.2dim.missing10.rda")
load("Dropbox/OC_publicopinion/analysis/fit.2dim.missing25.rda")
load("Dropbox/OC_publicopinion/analysis/fit.2dim.missing40.rda")
load("Dropbox/OC_publicopinion/analysis/fit.3dim.missing10.rda")
load("Dropbox/OC_publicopinion/analysis/fit.3dim.missing25.rda")
load("Dropbox/OC_publicopinion/analysis/fit.3dim.missing40.rda")
#
load("Dropbox/OC_publicopinion/analysis/fit.2dim.issues25.rda")
load("Dropbox/OC_publicopinion/analysis/fit.2dim.issues50.rda")
load("Dropbox/OC_publicopinion/analysis/fit.2dim.issues100.rda")
load("Dropbox/OC_publicopinion/analysis/fit.3dim.issues25.rda")
load("Dropbox/OC_publicopinion/analysis/fit.3dim.issues50.rda")
load("Dropbox/OC_publicopinion/analysis/fit.3dim.issues100.rda")
#
load("Dropbox/OC_publicopinion/analysis/fit.2dim.moderatesextremists.rda")
load("Dropbox/OC_publicopinion/analysis/fit.2dim.moderatesextremists.rda")
#
#
#
# 10%
fit <- fit.2dim.missing10
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.twodim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.twodim$errorlevel <- factor(xtile(fit.twodim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.twodim$cor.errorissues <- abs(fit.twodim$cor.errorissues)
#
names(fit.twodim)[names(fit.twodim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.twodim)[names(fit.twodim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.twodim)[names(fit.twodim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
twodim.df <- melt(fit.twodim, id=c("errorlevel"))
#
#ggplot(twodim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
twodim.df$Missing <- "10% Missing"
twodim.df10 <- twodim.df
#
# 25%
fit <- fit.2dim.missing25
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.twodim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.twodim$errorlevel <- factor(xtile(fit.twodim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.twodim$cor.errorissues <- abs(fit.twodim$cor.errorissues)
#
names(fit.twodim)[names(fit.twodim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.twodim)[names(fit.twodim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.twodim)[names(fit.twodim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
twodim.df <- melt(fit.twodim, id=c("errorlevel"))
#
#ggplot(twodim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
twodim.df$Missing <- "25% Missing"
twodim.df25 <- twodim.df
#
# 40%
fit <- fit.2dim.missing40
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.twodim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.twodim$errorlevel <- factor(xtile(fit.twodim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.twodim$cor.errorissues <- abs(fit.twodim$cor.errorissues)
#
names(fit.twodim)[names(fit.twodim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.twodim)[names(fit.twodim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.twodim)[names(fit.twodim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
twodim.df <- melt(fit.twodim, id=c("errorlevel"))
#
#ggplot(twodim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
twodim.df$Missing <- "40% Missing"
twodim.df40 <- twodim.df
#
#
#
twodim.df <- rbind(twodim.df10, twodim.df25, twodim.df40)
#
# Create line breaks
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr <- Vectorize(swr)
twodim.df$variable <- swr(twodim.df$variable)
#
twodim.df$variable <- factor(twodim.df$variable, levels=c("True and recovered\nideal points", "True and recovered\nissue normal\nvectors", "Reversed issue\nerror variance and\nrecovered issue PRE"))
#
pdf("Dropbox/OC_publicopinion/images/ooc_montecarlo_twodim.pdf", height=9, width=8)
#
ggplot(twodim.df, aes(errorlevel, value, fill=errorlevel)) +
	facet_wrap(~Missing+variable) +
	geom_boxplot() +
	xlab("") +
	ylab("Correlation\n") +
	scale_fill_viridis(discrete=TRUE, begin=0.8, end=0.2) +
	ylim(0,1) +
	guides(fill=FALSE) +
	ggtitle("Two Dimensions\n") +
	theme(plot.title = element_text(hjust = 0.5))
#
dev.off()
#
# 3-dimensional graph
#
# 10%
fit <- fit.3dim.missing10
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.threedim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.threedim$errorlevel <- factor(xtile(fit.threedim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.threedim$cor.errorissues <- abs(fit.threedim$cor.errorissues)
#
names(fit.threedim)[names(fit.threedim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.threedim)[names(fit.threedim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.threedim)[names(fit.threedim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
threedim.df <- melt(fit.threedim, id=c("errorlevel"))
#
#ggplot(threedim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
threedim.df$Missing <- "10% Missing"
threedim.df10 <- threedim.df
#
# 25%
fit <- fit.3dim.missing25
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.threedim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.threedim$errorlevel <- factor(xtile(fit.threedim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.threedim$cor.errorissues <- abs(fit.threedim$cor.errorissues)
#
names(fit.threedim)[names(fit.threedim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.threedim)[names(fit.threedim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.threedim)[names(fit.threedim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
threedim.df <- melt(fit.threedim, id=c("errorlevel"))
#
#ggplot(threedim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
threedim.df$Missing <- "25% Missing"
threedim.df25 <- threedim.df
#
# 40%
fit <- fit.3dim.missing40
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.threedim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.threedim$errorlevel <- factor(xtile(fit.threedim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.threedim$cor.errorissues <- abs(fit.threedim$cor.errorissues)
#
names(fit.threedim)[names(fit.threedim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.threedim)[names(fit.threedim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.threedim)[names(fit.threedim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
threedim.df <- melt(fit.threedim, id=c("errorlevel"))
#
#ggplot(threedim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
threedim.df$Missing <- "40% Missing"
threedim.df40 <- threedim.df
#
#
#
threedim.df <- rbind(threedim.df10, threedim.df25, threedim.df40)
#
# Create line breaks
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr <- Vectorize(swr)
threedim.df$variable <- swr(threedim.df$variable)
#
threedim.df$variable <- factor(threedim.df$variable, levels=c("True and recovered\nideal points", "True and recovered\nissue normal\nvectors", "Reversed issue\nerror variance and\nrecovered issue PRE"))
#
pdf("Dropbox/OC_publicopinion/images/ooc_montecarlo_threedim.pdf", height=9, width=8)
#
ggplot(threedim.df, aes(errorlevel, value, fill=errorlevel)) +
	facet_wrap(~Missing+variable) +
	geom_boxplot() +
	xlab("") +
	ylab("Correlation\n") +
	scale_fill_viridis(discrete=TRUE, begin=0.8, end=0.2) +
	ylim(0,1) +
	theme(plot.title = element_text(hjust = 0.5)) +
	guides(fill=FALSE) +
	ggtitle("Three Dimensions\n") +
	theme(plot.title = element_text(hjust = 0.5))
#
dev.off()
#
#table(fit[,"errorlevel"], fit.twodim$errorlevel)
#
#
#
#
# 25 issues
fit <- fit.2dim.issues25
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.twodim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.twodim$errorlevel <- factor(xtile(fit.twodim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.twodim$cor.errorissues <- abs(fit.twodim$cor.errorissues)
#
names(fit.twodim)[names(fit.twodim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.twodim)[names(fit.twodim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.twodim)[names(fit.twodim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
twodim.df <- melt(fit.twodim, id=c("errorlevel"))
#
#ggplot(twodim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
twodim.df$Missing <- "25 Issues"
twodim.df25 <- twodim.df
#
# 50 issues
fit <- fit.2dim.issues50
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.twodim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.twodim$errorlevel <- factor(xtile(fit.twodim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.twodim$cor.errorissues <- abs(fit.twodim$cor.errorissues)
#
names(fit.twodim)[names(fit.twodim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.twodim)[names(fit.twodim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.twodim)[names(fit.twodim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
twodim.df <- melt(fit.twodim, id=c("errorlevel"))
#
#ggplot(twodim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
twodim.df$Missing <- "50 Issues"
twodim.df50 <- twodim.df
#
# 100 issues
fit <- fit.2dim.issues100
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.twodim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.twodim$errorlevel <- factor(xtile(fit.twodim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.twodim$cor.errorissues <- abs(fit.twodim$cor.errorissues)
#
names(fit.twodim)[names(fit.twodim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.twodim)[names(fit.twodim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.twodim)[names(fit.twodim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
twodim.df <- melt(fit.twodim, id=c("errorlevel"))
#
#ggplot(twodim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
twodim.df$Missing <- "100 Issues"
twodim.df100 <- twodim.df
#
#
#
twodim.df <- rbind(twodim.df25, twodim.df50, twodim.df100)
#
# Create line breaks
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr <- Vectorize(swr)
twodim.df$variable <- swr(twodim.df$variable)
#
twodim.df$variable <- factor(twodim.df$variable, levels=c("True and recovered\nideal points", "True and recovered\nissue normal\nvectors", "Reversed issue\nerror variance and\nrecovered issue PRE"))
#
pdf("Dropbox/OC_publicopinion/images/ooc_montecarlo_twodim_numberissues.pdf", height=9, width=8)
#
ggplot(twodim.df, aes(errorlevel, value, fill=errorlevel)) +
	facet_wrap(~Missing+variable) +
	geom_boxplot() +
	xlab("") +
	ylab("Correlation\n") +
	scale_fill_viridis(discrete=TRUE, begin=0.8, end=0.2) +
	ylim(0,1) +
	guides(fill=FALSE) +
	ggtitle("Two Dimensions\n") +
	theme(plot.title = element_text(hjust = 0.5))
#
dev.off()
#
# 3-dimensional graph
#
# 25 issues
fit <- fit.3dim.issues25
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.threedim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.threedim$errorlevel <- factor(xtile(fit.threedim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.threedim$cor.errorissues <- abs(fit.threedim$cor.errorissues)
#
names(fit.threedim)[names(fit.threedim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.threedim)[names(fit.threedim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.threedim)[names(fit.threedim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
threedim.df <- melt(fit.threedim, id=c("errorlevel"))
#
#ggplot(threedim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
threedim.df$Missing <- "25 Issues"
threedim.df25 <- threedim.df
#
# 50 issues
fit <- fit.3dim.issues50
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.threedim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.threedim$errorlevel <- factor(xtile(fit.threedim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.threedim$cor.errorissues <- abs(fit.threedim$cor.errorissues)
#
names(fit.threedim)[names(fit.threedim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.threedim)[names(fit.threedim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.threedim)[names(fit.threedim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
threedim.df <- melt(fit.threedim, id=c("errorlevel"))
#
#ggplot(threedim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
threedim.df$Missing <- "50 Issues"
threedim.df50 <- threedim.df
#
# 100 issues
fit <- fit.3dim.issues100
colnames(fit) <- c("cor.idealpoints", "cor.normvecs", "cor.theta", "cor.errorrespondents", "cor.constraintrespondents", "cor.errorissues", "errorlevel")
fit.threedim <- data.frame(fit[,c("cor.idealpoints","cor.normvecs","cor.errorissues","errorlevel")])
#
fit.threedim$errorlevel <- factor(xtile(fit.threedim$errorlevel, 3), labels=c("Low error", "Medium error", "High error"))
fit.threedim$cor.errorissues <- abs(fit.threedim$cor.errorissues)
#
names(fit.threedim)[names(fit.threedim)=="cor.idealpoints"] <- "True and recovered ideal points"
names(fit.threedim)[names(fit.threedim)=="cor.normvecs"] <- "True and recovered issue normal vectors"
names(fit.threedim)[names(fit.threedim)=="cor.errorissues"] <- "Reversed issue error variance and recovered issue PRE"
#
threedim.df <- melt(fit.threedim, id=c("errorlevel"))
#
#ggplot(threedim.df, aes(errorlevel, value)) + facet_wrap(~variable) + geom_boxplot()
#
threedim.df$Missing <- "100 Issues"
threedim.df100 <- threedim.df
#
#
#
threedim.df <- rbind(threedim.df25, threedim.df50, threedim.df100)
#
# Create line breaks
swr = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr <- Vectorize(swr)
threedim.df$variable <- swr(threedim.df$variable)
#
threedim.df$variable <- factor(threedim.df$variable, levels=c("True and recovered\nideal points", "True and recovered\nissue normal\nvectors", "Reversed issue\nerror variance and\nrecovered issue PRE"))
#
pdf("Dropbox/OC_publicopinion/images/ooc_montecarlo_threedim_numberissues.pdf", height=9, width=8)
#
ggplot(threedim.df, aes(errorlevel, value, fill=errorlevel)) +
	facet_wrap(~Missing+variable) +
	geom_boxplot() +
	xlab("") +
	ylab("Correlation\n") +
	scale_fill_viridis(discrete=TRUE, begin=0.8, end=0.2) +
	ylim(0,1) +
	theme(plot.title = element_text(hjust = 0.5)) +
	guides(fill=FALSE) +
	ggtitle("Three Dimensions\n") +
	theme(plot.title = element_text(hjust = 0.5))
#
dev.off()
#
#table(fit[,"errorlevel"], fit.twodim$errorlevel)
#
#