######################################################################
# copyright (c) 2010-09-30, Katalin Csillery
#
#     These R functions are free and distributed in the hope that they
#     will be useful, but without any warranty; without even the
#     implied warranty of merchantability or fitness for a particular
#     purpose.
# 
#     These functions might be distributed in the future as part of an
#     R package for population genetic simulations and data analysis
# 
######################################################################

relat <- function(data, id=1, li=FALSE, wang=TRUE, queller=FALSE, ritland=FALSE, lynch=TRUE){
  ## id is 1 (ID's are in col 1) OR zero
  if(!(id %in% c(0,1))) print("id needs to be 0 (no id) or 1 (id in first column)")

  numind <- dim(data)[1]
  mycols <- dim(data)[2]
  if(id==1) id.pairs <- apply(make.index.pairs(1:numind), 1, function(a)
       paste(data[a[1], 1], data[a[2], 1], sep="-"))
  if(id==0) id.pairs <- apply(make.index.pairs(1:numind), 1, function(a)
       paste(a[1], a[2], sep="-"))

  data <- make.pairs(as.matrix(data[,(1+id):mycols]))
  n <- dim(data[[1]])[1]
  af <- lapply(data, function(a) as.numeric(table(a)/(4*n)))
  an <- lapply(data, function(a) as.numeric(names(table(a))))

  if(li) r.li <- Li(data, n, af, an) else r.li <- NULL
  if(wang) r.wang <- Wang(data, n, af, an) else r.wang <- NULL
  if(queller) r.queller <- Queller(data, n, af, an) else r.queller <- NULL
  if(ritland) r.ritland <- Ritland(data, n, af, an) else r.ritland <- NULL
  if(lynch) r.lynch <- Lynch(data, n, af, an) else r.lynch <- NULL
  
  res <- as.data.frame(cbind(r.li, r.wang, r.queller, r.ritland, r.lynch))
  names(res) <- c("li", "wang", "queller", "ritland", "lynch")[c(li, wang, queller, ritland, lynch)]
  rownames(res) <- id.pairs
  return(res)
}

Li <- function(gen, n, af, an){ # the Li estimator

    loci <- length(af)
    sims <- forLi(gen, n, af, an)
    
    relat <- function(a,sxy){ 
        s0 <- sum(a^2*(2-a))
        (sxy-s0)/(1-s0)
    }
    
    r1 <- lapply(af,relat,sxy=1)
    r2 <- lapply(af,relat,sxy=0.75)
    r3 <- lapply(af,relat,sxy=0.5)
    r4 <- lapply(af,relat,sxy=0)
    
    for (i in 1:loci){
        tmp <- sims[,i]
        tmp[tmp==1] <- r1[[i]]
        tmp[tmp==0.75] <- r2[[i]]
        tmp[tmp==0.5] <- r3[[i]]
        tmp[tmp==0] <- r4[[i]]
        sims[,i] <- tmp
    }
    
    r <- apply(sims, 1, mean)          
    
    return(r)
}# end of li function

## ############################################################
Lynch <- function(gen, n, af, an){ ## Lynch and Ritland estimator
    loci <- length(af)
    mat <- forLynch(gen, n, af, an)
    mat1 <- mat[[1]]
    mat2 <- mat[[2]]
    rm(mat)
    relat <- function(v){ # the Lynch and Ritland estimator
        (v[1]*(v[3] + v[7]) + v[2]*(v[6] + v[5]) - 4*v[1]*v[2]) / ((1 + v[4])*(v[1] + v[2]) - 4*v[1]*v[2])
    }
    weight <- function(v){
        ((1 + v[4])*(v[1] + v[2]) - 4*v[1]*v[2])/(2*v[1]*v[2])
    }
    ## for mat1 (reference individual x)
    r1 <- matrix(apply(mat1, 1, relat), ncol=loci)
    w1 <- matrix(apply(mat1, 1, weight), ncol=loci)
    W1 <- apply(w1, 1, sum)
    r1 <- r1*w1
    r1 <- apply(r1, 1, sum)/W1
    ## for mat2 (reference individual y)
    r2 <- matrix(apply(mat2, 1, relat), ncol=loci)
    w2 <- matrix(apply(mat2, 1, weight), ncol=loci)
    W2 <- apply(w2, 1, sum)
    r2 <- r2*w2
    r2 <- apply(r2, 1, sum)/W2
    ## and the average
    r <- rbind(c(r1), c(r2))
    r <- apply(r, 2, mean)
    return(r)
}#end of Lynch function
  
## ###################################################
Ritland <- function(gen, n, af, an){ ## Ritland estimator

    fun <- function(a, all){
        b <- apply(matrix(all, ncol=1), 1, function(x) match(a, unique(sort(x)), nomatch=0))
        apply(b, 2, function(x) ((x[1]+x[3])*(x[2]+x[4]))/4)
    }
    relat <- function(s, p, na){
        sum((s-p^2)/((na-1)*p))
    }
    loci <- length(af)
    pairs <- dim(gen[[1]])[1]
    num.all <- unlist(lapply(af, length))
    w <- sum(unlist(lapply(af, length))-1)

    S <- list()
    for(i in 1:loci){
        S[[i]] <- apply(gen[[i]], 1, fun, an[[i]])
    }
    
    r <- matrix(nrow=pairs, ncol=loci)
    for(i in 1:loci) r[,i] <- apply(S[[i]], 2, relat, p=af[[i]], na=num.all[i])*(num.all[i]-1)  
    apply(r, 1, function(a) sum(a)/w)*2
}# end of Ritland function

## ###################################################
Wang <- function(gen, n, af, an){ ## Wang estimator
    
    sims <- forWang(gen, n, af, an)
    loci <- length(af)

    relat <- function(a,p1,p2,p3,p4){ # the Wang estimator
        a2 <- sum(a^2)
        a3 <- sum(a^3)
        a4 <- sum(a^4)
        b <- 2*a2^2 - a4
        c <- a2 - 2*a2^2 +a4
        d <- 4*(a3 - a4)
        e <- 2*(a2 - 3*a3 + 2*a4)
        f <- 4*(a2 - a2^2 - 2*a3 +2*a4)
        g <- (1 - 7*a2 + 4*a2^2 + 10*a3 - 8*a4)
        V <- (1 - b)^2*(e^2*f + d*g^2) - (1 - b)*(e*f - d*g)^2 + 2*c*d*f*(1 - b)*(g + e) + c^2*d*f*(d + f)
        fi <- (d*f*( (e+g)*(1-b)+c*(d+f) )*(p1-1) + d*(1-b)*(g*(1-b-d)+f*(c+e))*p3 + f*(1-b)*(e*(1-b-f)+d*(c+g))*p2)/V
        delta <- (c*d*f*(e+g)*(p1+1-2*b) + ((1-b)*(f*e^2+d*g^2)-(e*f-d*g)^2)*(p1-b) + c*(d*g-e*f)*(d*p3-f*p2) - c^2*d*f * (p3+p2-d-f) - c*(1-b)*(d*g*p3+e*f*p2))/V
        r <- 0.5*fi + delta
        return(r)
    }
    
    weight <- function(a){ # locus specific weights
        a2 <- sum(a^2)
        a3 <- sum(a^3)
        2*a2 - a3 # returns U
    }
    chw <- function(a){
        w <- lapply(a, weight)
        1/(sum(1/unlist(w))*unlist(w))
    }
    w <- chw(af)
    
    r1 <- lapply(af,relat,p1=1,p2=0,p3=0,p4=0)
    r2 <- lapply(af,relat,p1=0,p2=1,p3=0,p4=0)
    r3 <- lapply(af,relat,p1=0,p2=0,p3=1,p4=0)
    r4 <- lapply(af,relat,p1=0,p2=0,p3=0,p4=1)
    
    for (i in 1:loci){
        tmp <- sims[,i]
        tmp[tmp==1] <- r1[[i]]
        tmp[tmp==0.75] <- r2[[i]]
        tmp[tmp==0.5] <- r3[[i]]
        tmp[tmp==0] <- r4[[i]]
        sims[,i] <- tmp*w[i]
    }
    
    return(apply(sims, 1, sum))
    
}# end of Wang function

Queller <- function(gen, n, af, an){ # the Queller and Goodnight estimator
    
    loci <- length(af)
    mat <- forQueller(gen, n, af, an)
    mat1 <- mat[[1]]
    mat2 <- mat[[2]]
    rm(mat)
    relat.denom <- function(v) 0.5*(v[6] + v[5] + v[3] + v[7]) - v[1] - v[2]
    relat.num <- function(v) 1 + v[4] - v[1] - v[2]
    
    ## for mat1 (reference individual x)
    r1.d <- matrix(apply(mat1, 1, relat.denom), ncol=loci)
    r1.n <- matrix(apply(mat1, 1, relat.num), ncol=loci)
    r1 <- apply(r1.d, 1, sum)/apply(r1.n, 1, sum)
    ## for mat2 (reference individual y)
    r2.d <- matrix(apply(mat2, 1, relat.denom), ncol=loci)
    r2.n <- matrix(apply(mat2, 1, relat.num), ncol=loci)
    r2 <- apply(r2.d, 1, sum)/apply(r2.n, 1, sum)
    ## and the average of those:
    r <- rbind(c(r1), c(r2))
    r <- apply(r, 2, mean)
    return(r)
    
}#end of Queller function

## ####################################################################

## misc: functions for internal use by some of the function above

## ####################################################################

make.pairs <- function(x){

  ind <- dim(x)[1]
  loci <- dim(x)[2]/2
  if(loci==1) y <- matrix(x, ncol=1)
  else y <- rbind(x[,(2*(1:loci)-1)], x[,(2*(1:loci))])
  
  fun <- function(a){
    a <- cbind(a[1:ind], a[(ind+1):(2*ind)])
    as.matrix(cbind(a[rep(1:(ind-1), times=(ind-1):1),1],
          a[sequence(c((ind-1):1))+rep(1:(ind-1), times=(ind-1):1),1],
          a[rep(1:(ind-1), times=(ind-1):1),2],
          a[sequence(c((ind-1):1))+rep(1:(ind-1), times=(ind-1):1),2]))
  }
  l <- list()
  for(i in 1:loci) l[[i]] <- fun(y[,i])
  return(l)
} # end of make.pairs

make.index.pairs <- function(x){

  ind <- length(x)
  
  fun <- function(a){
    a <- cbind(a[1:ind], a[(ind+1):(2*ind)])
    as.matrix(cbind(a[rep(1:(ind-1), times=(ind-1):1),1],
          a[sequence(c((ind-1):1))+rep(1:(ind-1), times=(ind-1):1),1]))
  }
  fun(x)
  
} # end of make.index.pairs

forLi <- function(gen, n, af, an){ # gen is a list
    sim <- matrix(NA, n, length(af))
    for(i in 1:length(af)){
        sim[,i] <- (4-(abs(sign(gen[[i]][,1]-gen[[i]][,2]))+abs(sign(gen[[i]][,3]-gen[[i]][,4])))*(abs(sign(gen[[i]][,1]-gen[[i]][,4]))+abs(sign(gen[[i]][,3]-gen[[i]][,2]))))/4
    }
    return(sim)
}# end of forL

forLynch <- function(gen, n, af, an){ # gen is a list
  sim1 <- matrix(NA, n*length(af), 7)
  sim2 <- matrix(NA, n*length(af), 7)
  for(i in 1:length(af)){
    s1 <- list(abs(sign(gen[[i]][,3]-gen[[i]][,2])), abs(sign(gen[[i]][,1]-gen[[i]][,3])), abs(sign(gen[[i]][,1]-gen[[i]][,4])), abs(sign(gen[[i]][,1]-gen[[i]][,2])), abs(sign(gen[[i]][,3]-gen[[i]][,4])))
    
    s2 <- list(abs(sign(gen[[i]][,1]-gen[[i]][,4])), abs(sign(gen[[i]][,2]-gen[[i]][,4])), abs(sign(gen[[i]][,2]-gen[[i]][,3])), abs(sign(gen[[i]][,1]-gen[[i]][,2])), abs(sign(gen[[i]][,3]-gen[[i]][,4])))
    
    sim1[((i-1)*n+1):(i*n),3:7] <- trunc(cos(matrix(unlist(s1),n,5)*(2/pi)))
    sim1[((i-1)*n+1):(i*n),1:2] <- cbind(af[[i]][match(gen[[i]][,1],an[[i]])],af[[i]][match(gen[[i]][,3],an[[i]])])
    sim2[((i-1)*n+1):(i*n),3:7] <- trunc(cos(matrix(unlist(s2),n,5)*(2/pi)))
    sim2[((i-1)*n+1):(i*n),1:2] <- cbind(af[[i]][match(gen[[i]][,2],an[[i]])],af[[i]][match(gen[[i]][,4],an[[i]])])
  }
  return(list(sim1, sim2))
}

forWang <- function(gen, n, af, an){ # gen is a list
  sim <- matrix(NA, n, length(af))
  for(i in 1:length(af)){
    sim[,i] <- (4-(abs(sign(gen[[i]][,1]-gen[[i]][,2]))+abs(sign(gen[[i]][,3]-gen[[i]][,4])))*(abs(sign(gen[[i]][,1]-gen[[i]][,4]))+abs(sign(gen[[i]][,3]-gen[[i]][,2]))))/4
  }
  return(sim)
}# end of forWang

forQueller <- function(gen, n, af, an){ # gen is a list
  sim1 <- matrix(NA, n*length(af), 7)
  sim2 <- matrix(NA, n*length(af), 7)
  for(i in 1:length(af)){
    s1 <- list(abs(sign(gen[[i]][,2]-gen[[i]][,3])), abs(sign(gen[[i]][,1]-gen[[i]][,3])), abs(sign(gen[[i]][,1]-gen[[i]][,4])), abs(sign(gen[[i]][,1]-gen[[i]][,2])), abs(sign(gen[[i]][,3]-gen[[i]][,4])))

    s2 <- list(abs(sign(gen[[i]][,1]-gen[[i]][,4])), abs(sign(gen[[i]][,2]-gen[[i]][,4])), abs(sign(gen[[i]][,2]-gen[[i]][,3])), abs(sign(gen[[i]][,1]-gen[[i]][,2])), abs(sign(gen[[i]][,3]-gen[[i]][,4])))
    
    sim1[((i-1)*n+1):(i*n),3:7] <- trunc(cos(matrix(unlist(s1),n,5)*(2/pi)))
    sim1[((i-1)*n+1):(i*n),1:2] <- cbind(af[[i]][match(gen[[i]][,1],an[[i]])],af[[i]][match(gen[[i]][,3],an[[i]])])
    sim2[((i-1)*n+1):(i*n),3:7] <- trunc(cos(matrix(unlist(s2),n,5)*(2/pi)))
    sim2[((i-1)*n+1):(i*n),1:2] <- cbind(af[[i]][match(gen[[i]][,2],an[[i]])],af[[i]][match(gen[[i]][,4],an[[i]])])
  }
  return(list(sim1, sim2))
} # end of of forQueller 
