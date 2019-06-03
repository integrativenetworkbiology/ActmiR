

InfermiRactivity <- function(miRNA, miRexpr, expr, target, cutoff) {
	act <- list()
	sd <- list()
	ww <- list()
	weights_dff <- c()
	gg <- intersect(target, rownames(expr))
	if(length(which( miRexpr[miRNA,]==min(miRexpr)))>3){
		basesample <- which(miRexpr[miRNA,]==min(miRexpr))
	} else {
		basesample <- order(miRexpr[miRNA,])[1:3]
	}
	bb <-rowMeans(expr[gg, basesample])
	deg <- apply(expr[gg,], 2, function(x) bb-x)
	act[[1]] <- apply(deg, 2, function(x) { tmp <- summary(lm(x~bb))$coefficients
                                                if(dim(tmp)[1]<2) { NA } else { tmp[2,1]}})
	cc <- cor(act[[1]],t(expr[gg,]))
	weights <- as.numeric(cc)
	weights[which(weights>(-1*cutoff))] <- cutoff
	weights <- abs(weights)
	ww[[1]] <- weights
	iterations <- 20
	for( s in 2:iterations ) {
        	deg <- apply(expr[gg,], 2, function(x) bb-x)
	        act[[s]] <- apply(deg, 2, function(x) { tmp <- summary(lm(x~bb, weights=ww[[s-1]]))$coefficients
        	                                        if(dim(tmp)[1]<2) { NA } else { tmp[2,1]}})
	        cc <- cor(act[[s]],t(expr[gg,]))
	        weights <- as.numeric(cc)
	        weights[which(weights>(-1*cutoff))] <- cutoff
	        weights <- abs(weights)
	        ww[[s]] <- weights
	        tmp <- (ww[[s]]-ww[[s-1]])
	        tmp <- tmp[-which(tmp==0|is.na(tmp))]
	        if(length(tmp)==0) { weights_dff[s-1]  <- 0 } else {weights_dff[s-1] <- mean((tmp^2)^(1/2), na.rm=T)}
	        if(s==iterations) {
	                activity_wls <- act[[s]]
	        } else if(weights_dff[s-1] < 0.005) {
	                activity_wls <- act[[s]]
	                break
        	} 
	}
	activity_wls
}



