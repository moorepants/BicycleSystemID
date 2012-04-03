#data <- read.csv('whipple-structured-results.csv')
data <- read.csv('whipple-structured-results-x0-K-freeones.csv')

data <- subset(data, MeanFit > 30)

library(nlme)

# fit all the constant coefficients
ls.a31 = lm(a31 ~ 1, data)
ls.a41 = lm(a41 ~ 1, data)
ls.b31 = lm(b31 ~ 1, data)
ls.b41 = lm(b41 ~ 1, data)

me.a31 = lme(a31 ~ 1, data, ~1 | Speed)
me.a41 = lme(a41 ~ 1, data, ~1 | Speed)
me.b31 = lme(b31 ~ 1, data, ~1 | Speed)
me.b41 = lme(b41 ~ 1, data, ~1 | Speed)

# fit the linear coefficients
ls.a33 = lm(a33 ~ ActualSpeed, data)
ls.a34 = lm(a34 ~ ActualSpeed, data)
ls.a43 = lm(a43 ~ ActualSpeed, data)
ls.a44 = lm(a44 ~ ActualSpeed, data)

me.a33 = lme(a33 ~ ActualSpeed, data, ~1 | Speed)
me.a34 = lme(a34 ~ ActualSpeed, data, ~1 | Speed)
me.a43 = lme(a43 ~ ActualSpeed, data, ~1 | Speed)
me.a44 = lme(a44 ~ ActualSpeed, data, ~1 | Speed)

# fit the quadratic coefficients
ls.a32 = lm(a32 ~ I(ActualSpeed^2), data)
ls.a42 = lm(a42 ~ I(ActualSpeed^2), data)

me.a32 = lme(a32 ~ I(ActualSpeed^2), data, ~1 | Speed)
me.a42 = lme(a42 ~ I(ActualSpeed^2), data, ~1 | Speed)

# plot the fitted curves
x <- seq(min(data$ActualSpeed), max(data$ActualSpeed), num=200)
png('whipple-structured-x0-K-freeones.png', width=1000, height=600)
layout(matrix(1:10, 2, 5))

plot(a31 ~ ActualSpeed, data)
abline(ls.a31)
abline(a=me.a31$coef$fixed, b=0, col="red")

plot(a41 ~ ActualSpeed, data)
abline(ls.a41)
abline(a=me.a41$coef$fixed, b=0, col="red")

plot(a32 ~ ActualSpeed, data)
y <- ls.a32$coef %*% rbind(1, x^2)
lines(x, y)
y <- me.a32$coef$fixed %*% rbind(1, x^2)
lines(x, y, col="red")

plot(a42 ~ ActualSpeed, data)
y <- ls.a42$coef %*% rbind(1, x^2)
lines(x, y)
y <- me.a42$coef$fixed %*% rbind(1, x^2)
lines(x, y, col="red")

plot(a33 ~ ActualSpeed, data)
abline(ls.a33)
abline(me.a33$coef$fixed, col="red")
plot(a43 ~ ActualSpeed, data)
abline(ls.a43)
abline(me.a43$coef$fixed, col="red")

plot(a34 ~ ActualSpeed, data)
abline(ls.a34)
abline(me.a34$coef$fixed, col="red")
plot(a44 ~ ActualSpeed, data)
abline(ls.a44)
abline(me.a44$coef$fixed, col="red")

plot(b31 ~ ActualSpeed, data)
abline(ls.b31)
abline(a=me.b31$coef$fixed, b=0, col="red")
plot(b41 ~ ActualSpeed, data)
abline(ls.b41)
abline(a=me.b41$coef$fixed, b=0, col="red")

dev.off()

png('../plots/empirical-root-loci-basic-regression.png')
speed <- seq(0, 10, 0.1)
plot(0, 0, xlim=c(0, 10), ylim=c(-10, 10))
for (v in speed) {
  A <- rbind(c(0, 0, 1, 0), c(0, 0, 0, 1),
    c(ls.a31$coef, ls.a32$coef %*% rbind(1, v^2), ls.a33$coef %*% rbind(1, v),
      ls.a34$coef %*% rbind(1, v)),
    c(ls.a41$coef, ls.a42$coef %*% rbind(1, v^2), ls.a43$coef %*% rbind(1, v),
      ls.a44$coef %*% rbind(1, v)))
  e <- eigen(A, only.values=TRUE)
  points(c(v, v, v, v), Re(e$values), col='black')
  points(c(v, v, v, v), Im(e$values), col='blue')
}
dev.off()

png('../plots/empirical-root-loci-mixed-effects.png')
plot(0, 0, xlim=c(0, 10), ylim=c(-10, 10))
for (v in speed) {
  A <- rbind(c(0, 0, 1, 0), c(0, 0, 0, 1),
    c(me.a31$coef$fixed, me.a32$coef$fixed %*% rbind(1, v^2), me.a33$coef$fixed %*% rbind(1, v),
      me.a34$coef$fixed %*% rbind(1, v)),
    c(me.a41$coef$fixed, me.a42$coef$fixed %*% rbind(1, v^2), me.a43$coef$fixed %*% rbind(1, v),
      me.a44$coef$fixed %*% rbind(1, v)))
  e <- eigen(A, only.values=TRUE)
  points(c(v, v, v, v), Re(e$values), col='black')
  points(c(v, v, v, v), Im(e$values), col='blue')
}
dev.off()
