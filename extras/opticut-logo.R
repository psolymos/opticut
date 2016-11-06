library(opticut)

y <- t(matrix(c(
    1,	1,	1,	1,	1,	1,	1,	1,	1,	1,
    1,	1,	0,	0,	1,	1,	1,	0,	0,	1,
    1,	0,	1,	1,	0,	1,	0,	1,	1,	1,
    1,	0,	1,	1,	0,	1,	0,	1,	1,	1,
    1,	0,	1,	1,	0,	1,	0,	1,	1,	1,
    1,	1,	0,	0,	1,	1,	1,	0,	0,	1,
    1,	1,	1,	1,	1,	1,	1,	1,	1,	1),
    7, 10, byrow=TRUE))
colnames(y) <- c("O","P","T","I","C","U","t")
yy <- 1-y*0.9
x <- opticut(yy ~ 1, strata=c("I","N","D","i","C","A","T","O","R","S"))
plot(x, sort=FALSE, cut=-Inf, xlab="", ylab="", show_I=FALSE, show_S=FALSE)

library(animation)
zzz <- seq(0.05, 0.95, 0.15)
saveGIF({
    for (i in c(zzz, rev(zzz))) {
        yy <- 1-y*i
        x <- opticut(yy ~ 1, strata=c("I","N","D","i","C","A","T","O","R","S"))
        plot(x, sort=FALSE, cut=-Inf,
            xlab="", ylab="", show_I=FALSE, show_S=FALSE)
    }
}, ani.width=600, ani.height=400, interval=0.2)
