#========================================================#
# STAT 6021: Linear Models for Data Science              |
# R tutorial for learning unit 8                         |
# R functions                                            |
#========================================================#
#                                                        |
# -------------------------------------------------------+
# R function for ridge regression                        |
# -------------------------------------------------------+

ridge.reg.coefficients <- function(y.vect, X0.mat, plot=TRUE, grid.size=25, grid.st=0.001, grid.fn=0.5) {
	# Collect parameters
	n <- dim(X0.mat)[1]
	k <- dim(X0.mat)[2]
	p <- k + 1
	# Unit-length scaling
	y.bar <- mean(y.vect)
	y.cent <- y.vect - y.bar
	SS.T <- sum(y.cent^2)
	y.vect.scl <- y.vect / sqrt(SS.T)
	X.mat.scl <- matrix(data=NA, nrow=n, ncol=k)
	x.bar.list <- numeric(length=k)
	css.list <- numeric(length=k)
	for (j in 1:k) {
		x.bar.list[j] <- mean(X0.mat[,j])
		xj.cent <- X0.mat[,j] - x.bar.list[j]
		css.list[j] <- sum(xj.cent^2)
		X.mat.scl[,j] <- xj.cent / sqrt(css.list[j])
	}
	# Calculate ridge trace diagram
	ridge.k.grid <- exp(seq(from=log(grid.st), to=log(grid.fn), length.out=grid.size))
	b.hat.R.scl.list <- matrix(data=0, nrow=k, ncol=grid.size)
	y.vect.aug <- rbind(y.vect.scl, matrix(data=0, nrow=k, ncol=1))
	for (iVAL in 1:grid.size) {
		ridge.k.val <- ridge.k.grid[iVAL]
		X.mat.aug <- rbind(X.mat.scl, sqrt(ridge.k.val)*diag(k))
		XpX.mat.aug <- t(X.mat.aug) %*% X.mat.aug
		Xpy.mat.aug <- t(X.mat.aug) %*% y.vect.aug
		XpX.inv.aug <- solve(XpX.mat.aug)
		b.hat.R.scl.list[,iVAL] <- XpX.inv.aug %*% Xpy.mat.aug
	}
	if (plot) {
		plot(ridge.k.grid, rep(x=0, times=grid.size), pch=3, cex=1, ylim=c(min(b.hat.R.scl.list), max(b.hat.R.scl.list)), xlab="ridge constant, k", ylab="fitted ridge regression coefficient", main = "Ridge trace diagram")
		abline(h=0, lty=1, lwd=1)
		for (j in 1:k) {
			lines(ridge.k.grid, b.hat.R.scl.list[j,], type="l", lty=1, lwd=3)
		}
	}
	# Convert to the original scale and calculate MS.Res and R2.
	X.mat <- as.matrix(cbind(rep(x=1, times=n), X0.mat))
	b.hat.R.list <- matrix(data=0, nrow=p, ncol=grid.size)
	SS.Res.list <- numeric(length=grid.size)
	R2.list <- numeric(length=grid.size)
	for (iVAL in 1:grid.size) {
		b.hat.R.list[1,iVAL] <- y.bar
		for (j in 1:k) {
			b.hat.R.list[j+1,iVAL] <- b.hat.R.scl.list[j,iVAL] / sqrt(css.list[j] / SS.T)
			b.hat.R.list[1,iVAL] <- b.hat.R.list[1,iVAL] - b.hat.R.list[j+1,iVAL]*x.bar.list[j]
		}
		SS.Res.list[iVAL] <- sum((y.vect - X.mat %*% b.hat.R.list[,iVAL])^2)
		R2.list[iVAL] <- 1 - SS.Res.list[iVAL] / SS.T
	}
	MS.Res.list <- SS.Res.list / (n-p)
	out.list <- list(ridge.k.grid=ridge.k.grid, b.hat.R.list=b.hat.R.list, MS.Res.list=MS.Res.list, R2.list=R2.list)
	return(out.list)
}

# -------------------------------------------------------+
# R function for principle components regression         |
# -------------------------------------------------------+

prin.comp.coefficients <- function(y.vect, X0.mat) {
	# Collect parameters
	n <- dim(X0.mat)[1]
	k <- dim(X0.mat)[2]
	p <- k + 1
	# Unit-length scaling
	y.bar <- mean(y.vect)
	y.cent <- y.vect - y.bar
	SS.T <- sum(y.cent^2)
	y.vect.scl <- y.vect / sqrt(SS.T)
	X.mat.scl <- matrix(data=NA, nrow=n, ncol=k)
	x.bar.list <- numeric(length=k)
	css.list <- numeric(length=k)
	for (j in 1:k) {
		x.bar.list[j] <- mean(X0.mat[,j])
		xj.cent <- X0.mat[,j] - x.bar.list[j]
		css.list[j] <- sum(xj.cent^2)
		X.mat.scl[,j] <- xj.cent / sqrt(css.list[j])
	}
	# Calculate principal components and coefficient estimates
	XpX.mat.scl <- t(X.mat.scl) %*% X.mat.scl
	eig.out <- eigen(XpX.mat.scl)
	Lambda.mat <- diag(eig.out$values)
	Lambda.inv <- diag(1/eig.out$values)
	T.mat <- eig.out$vectors
	Z.mat <- X.mat.scl %*% T.mat
	Zpy.mat <- t(Z.mat) %*% y.vect.scl
	a.hat.scl <- Lambda.inv %*% Zpy.mat[1:j,]
	a.hat.PC.scl.list <- matrix(data=0, nrow=k, ncol=k)
	b.hat.PC.scl.list <- matrix(data=0, nrow=k, ncol=k)
	for (j in 1:k) {
		a.hat.PC.scl.list[1:j,j] <- a.hat.scl[1:j,]
		b.hat.PC.scl.list[,j] <- T.mat %*% a.hat.PC.scl.list[,j]
	}
	# Convert to the original scale and calculate MS.Res and R2.
	X.mat <- as.matrix(cbind(rep(x=1, times=n), X0.mat))
	grid.size <- dim(b.hat.PC.scl.list)[2]
	b.hat.PC.list <- matrix(data=0, nrow=p, ncol=grid.size)
	SS.Res.list <- numeric(length=grid.size)
	R2.list <- numeric(length=grid.size)
	for (iVAL in 1:grid.size) {
		b.hat.PC.list[1,iVAL] <- y.bar
		for (j in 1:k) {
			b.hat.PC.list[j+1,iVAL] <- b.hat.PC.scl.list[j,iVAL] / sqrt(css.list[j] / SS.T)
			b.hat.PC.list[1,iVAL] <- b.hat.PC.list[1,iVAL] - b.hat.PC.list[j+1,iVAL]*x.bar.list[j]
		}
		SS.Res.list[iVAL] <- sum((y.vect - X.mat %*% b.hat.PC.list[,iVAL])^2)
		R2.list[iVAL] <- 1 - SS.Res.list[iVAL] / SS.T
	}
	MS.Res.list <- SS.Res.list / (n-p)
	out.list <- list(b.hat.PC.list=b.hat.PC.list, MS.Res.list=MS.Res.list, R2.list=R2.list)
	return(out.list)
}
