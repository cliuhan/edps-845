#' @title Equipercentile Equating Analytic Standard Error
#' 
#' @description This function is used to calcuate the analytic standard error
#' of equipercentile equating without smoothing under single- or random-
#' groups design.
#' 
#' @details The delta method is a commonly used to compute standard errors.
#' 
#' @param x and y are frequency tables.
#' @examples
#' 
#' Get the analytic standard errors. 
#' 
#' library(equate)
#' rx <- as.freqtab(ACTmath[, 1:2])
#' ry <- as.freqtab(ACTmath[, c(1, 3)])
#' mod <- equip(x, y)
#' se <- mod$se
#' 
#' @export
equip <- function (x, y, ly = min(scales(y)), ky = max(scales(y))) {
	yscale <- scales(y)
	yn <- sum(y)
	if (!equate:::is.freqtab(x)) {
		prank <- sort(unique(x))
		xscale <- yscale
		xn <- 0
	}
	else {
		prank <- round(px(x), 10)
		xscale <- scales(x)
		xn <- sum(x)
	}
	sn <- length(yscale)
	yinc <- round(diff(yscale), 8)
	yincl <- c(yinc[1]/2, yinc/2)
	yinch <- c(yinc/2, yinc[sn - 1]/2)
	yx <- numeric(length = length(prank))
	fy <- round(fx(y), 10)
	xnone <- prank == 0
	xone <- prank == 1
	xbot <- sum(xnone) + 1
	xtop <- sum(!xone)
	yxi <- xbot:xtop
	xyone <- which(xscale[xone] > (ky + yinch[sn])) + xtop
	yx[xnone] <- ly - yincl[1]
	yx[xone] <- ky + yinch[sn]
	yx[xyone] <- xscale[xyone]
	if (any(yx == 0)) {
		yu <- sapply(yxi, function(i) sum(fy <= prank[i]) + 1)
		yu2 <- yu - 1
		yu2[yu2 > 0] <- fy[yu2]
		g0 <- fy[yu] - yu2
		yx[yxi] <- yscale[yu] - yincl[yu] + ((prank[yxi] - yu2)/g0) * 
			(yinch + yincl)[yu]
		if (any(y == 0)) {
			yxi <- (xbot + sum(prank[!xnone] <= min(fy))):xtop
			yl <- sapply(yxi, function(i) sum(fy < prank[i]))
			yl2 <- fy[yl + 1]
			yxtemp <- yscale[yl] + yincl[yu] + ((prank[yxi] - 
					fy[yl])/(yl2 - fy[yl])) * (yinch + yincl)[yu]
			yx[yxi] <- (yx[yxi] + yxtemp)/2
		}
	}
	# Get se
	se <- numeric(length = length(prank))
	se[xnone] <- NA
	se[yxi] <- sqrt((1/(g0^2)) * (prank[yxi] * (1 - prank[yxi]) * (xn + yn)/(xn * yn) - 
			(fy[yu] - prank[yxi]) * (prank[yxi] - yu2)/(yn * g0)))
	
	if (xn) {
		out <- list(yx = yx, se = se)
	}
	else out <- list(yx = yx[match(x, prank)], se = se[match(x, prank)])
	return(out)
}
