### R code from vignette source 'rococo.Rnw'

###################################################
### code chunk number 1: Init
###################################################
options(width=75)
set.seed(0)
library(rococo)
rococoVersion <- packageDescription("rococo")$Version
rococoDateRaw <- packageDescription("rococo")$Date
rococoDateYear <- as.numeric(substr(rococoDateRaw, 1, 4))
rococoDateMonth <- as.numeric(substr(rococoDateRaw, 6, 7))
rococoDateDay <- as.numeric(substr(rococoDateRaw, 9, 10))
rococoDate <- paste(month.name[rococoDateMonth], " ",
                     rococoDateDay, ", ",
                     rococoDateYear, sep="")


###################################################
### code chunk number 2: InstallRoCoCo (eval = FALSE)
###################################################
## install.packages("rococo")


###################################################
### code chunk number 3: LoadRoCoCo (eval = FALSE)
###################################################
## library(rococo)


###################################################
### code chunk number 4: OpenVignette (eval = FALSE)
###################################################
## vignette("rococo")


###################################################
### code chunk number 5: ShowHelp (eval = FALSE)
###################################################
## help(rococo)


###################################################
### code chunk number 6: ToyData
###################################################
x1 <- rnorm(15)
y1 <- 2 * x1 + rnorm(length(x1), sd=0.25)
plot(x1, y1)


###################################################
### code chunk number 7: SimpleRococo
###################################################
rococo(x1, y1)


###################################################
### code chunk number 8: SimpleRococoTest
###################################################
rococo.test(x1, y1, alternative="two.sided")


###################################################
### code chunk number 9: RoCoCoTestFormula
###################################################
data(iris)
plot(~ Sepal.Length + Petal.Length, iris)
rococo.test(~ Sepal.Length + Petal.Length, iris, alternative="two.sided")


###################################################
### code chunk number 10: FlatNoisyData
###################################################
x2 <- rnorm(15)
f2 <- function(x) ifelse(x > 0.8, x - 0.8, ifelse(x < -0.8, x + 0.8, 0))
y2 <- f2(x2) + rnorm(length(x2), sd=0.1)
plot(x2, y2)


###################################################
### code chunk number 11: FlatNoisyDefault
###################################################
rococo.test(x2, y2, alternative="greater")


###################################################
### code chunk number 12: FlatNoisyOther
###################################################
rococo.test(x2, y2, similarity="gauss", alternative="greater")
rococo.test(x2, y2, similarity=c("classical", "gauss"), alternative="greater")


###################################################
### code chunk number 13: FlatNoisySameR
###################################################
rococo.test(x2, y2, similarity="gauss", r=0.1, alternative="greater")


###################################################
### code chunk number 14: FlatNoisyDifferentR
###################################################
rococo.test(x2, y2, similarity=c("linear", "gauss"), r=c(0.05, 0.1),
            alternative="greater")


###################################################
### code chunk number 15: CheckIQR
###################################################
rococo.test(x2, y2, similarity=c("linear", "gauss"), r=0, alternative="greater")
IQR(x2) * 0.1
IQR(y2) * 0.1


###################################################
### code chunk number 16: CheckIQR
###################################################
rococo.test(x2, y2, similarity=c("linear", "gauss"), tnorm="prod",
            alternative="greater")


###################################################
### code chunk number 17: YagertNorm
###################################################
DrastictNorm <- function(x, y)
{
    if (x == 1) y
    else if (y == 1) x
    else 0
}
YagertNorm <- function(lambda)
{
    fun <- function(x, y)
    {
        if (lambda == 0)
            DrastictNorm(x, y)
        else if (is.infinite(lambda))
            min(x, y)
        else
            max(0, 1 - ((1 - x)^lambda + (1 - y)^lambda)^(1 / lambda))
    }

    attr(fun, "name") <- paste("Yager t-norm with lambda =", lambda)

    fun
}
rococo(x2, y2, tnorm=YagertNorm(0.5))
rococo.test(x2, y2, tnorm=YagertNorm(0.2))


###################################################
### code chunk number 18: LargeNoShuffles
###################################################
res <- rococo.test(x2, y2, numtests=100000, storeValues=TRUE)
res


###################################################
### code chunk number 19: ExactpValueExample
###################################################
rococo.test(x2[1:8], y2[1:8], exact=TRUE)


###################################################
### code chunk number 20: LargeNoShufflesPic
###################################################
hist(res@perm.gamma, breaks=100, probability=TRUE, xlab="gamma",
     main="Distribution of gamma for random shuffles")
plot(function(x) dnorm(x, mean=res@H0gamma.mu, sd=res@H0gamma.sd),
     min(res@perm.gamma), max(res@perm.gamma), col="red", lwd=2, add=TRUE)


###################################################
### code chunk number 21: LargeNoShufflespVal
###################################################
res@p.value.approx


###################################################
### code chunk number 22: SimpleGaussCorTest
###################################################
gauss.cor(x1, y1)
gauss.cor.test(x1, y1, alternative="two.sided")
gauss.cor.test(~ Sepal.Length + Petal.Length, iris, alternative="two.sided")


###################################################
### code chunk number 23: GetBibTeX (eval = FALSE)
###################################################
## toBibtex(citation("rococo"))


