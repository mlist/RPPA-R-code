#######################################################
#######################################################
# estimateGamma
# This function estimates the "gamma" parameter in the variable
# slope normalization model for RPPAs
# INPUTS
#    	Xhat: A matrix with samples in the rows and proteins in the columns
#			  to be used.  This matrix should have already been adjusted for
#			  for column effects (the column median should have already been subtracted off) 
#		method: The regression method that is used to estimate the "gammas".  The default
#				is "pca" which uses perpendicular regression.  There is also a method
# 				that just uses regular least squares regression.  This first estimates 
#				the ratio and then estimates the inverse, but yields similar results to
#				perpendicular regression. 


estimateGamma <- function(Xhat,method="pca") {
  
  #####################
  # use perpendicular regression
  if(method=="pca") {
    nCol <- ncol(Xhat)
    gamma <- matrix(0,nrow=nCol,ncol=nCol)
    means <- apply(Xhat,2,mean)
    
    for (i in 1:(nCol-1)) {
      for (j in (i+1):nCol) {
        r <- cor(Xhat[,i],Xhat[,j],use="complete.obs")
        a <- Xhat[,i]
        n <- length(a)
        tt <- r*sqrt((n-2)/(1-r^2))
        chk <- pt(tt,n-2,lower.tail=F)
        if (chk < .05) {
          eig <- eigen(var(cbind(Xhat[,i],Xhat[,j]),na.rm=T))
          tmp <- (-1)*eig$vectors[1,2]/eig$vectors[2,2]
          gamma[i,j] <- tmp
        }
      }
    }
    gamma[gamma<=0] <- 1
    upper <- upper.tri(gamma)
    ind <- which(upper,arr.ind=T)
    
    design <- matrix(0,ncol=nCol,nrow=nrow(ind))
    for(i in 1:nrow(ind)) {
      design[i,ind[i,1]] <- -1
      design[i,ind[i,2]] <- 1
    }
    
    loggamma <- log(gamma[upper])
    
    newrow <- rep((1/nCol),nCol)
    nonsingular <- rbind(newrow,design)
    lestimateMean <- qr.solve(nonsingular,c(0,loggamma))
    
    estimate1 <- exp(lestimateMean)
    val  <- estimate1
    
  }
  ####################
  # use least squares
  if(method!="pca") {
    nCol <- ncol(Xhat)
    gamma <- matrix(0,nrow=nCol,ncol=nCol)
    means <- apply(Xhat,2,mean,na.rm=T)
    
    for (i in 1:(nCol-1)) {
      for (j in (i+1):nCol) {
        gamma[i,j] <- 
          sum((Xhat[,i]-means[i])*(Xhat[,j]-means[j]),na.rm=T)/sum((Xhat[,i]-means[i])^2,na.rm=T)
      }
    }
    gamma[gamma<0] <- 1
    upper <- upper.tri(gamma)
    ind <- which(upper,arr.ind=T)
    
    design <- matrix(0,ncol=nCol,nrow=nrow(ind))
    for(i in 1:nrow(ind)) {
      design[i,ind[i,1]] <- -1
      design[i,ind[i,2]] <- 1
    }
    
    loggamma <- log(gamma[upper])
    
    newrow <- rep((1/nCol),nCol)
    nonsingular <- rbind(newrow,design)
    lestimateMean <- qr.solve(nonsingular,c(0,loggamma))
    
    estimate1 <- exp(lestimateMean)
    
    gamma <- matrix(0,nrow=nCol,ncol=nCol)
    means <- apply(Xhat,2,mean,na.rm=T)
    for (i in 1:(nCol-1)) {
      for (j in (i+1):nCol) {
        gamma[i,j] <- 
          sum((Xhat[,j]-means[j])*(Xhat[,i]-means[i]),na.rm=T)/sum((Xhat[,j]-means[j])^2,na.rm=T)
      }
    }
    gamma[gamma<0] <- 1
    upper <- upper.tri(gamma)
    ind <- which(upper,arr.ind=T)
    
    design <- matrix(0,ncol=nCol,nrow=nrow(ind))
    for(i in 1:nrow(ind)) {
      design[i,ind[i,1]] <- 1
      design[i,ind[i,2]] <- -1
    }
    
    loggamma <- log(gamma[upper])
    
    newrow <- rep((1/nCol),nCol)
    nonsingular <- rbind(newrow,design)
    lestimateMean <- qr.solve(nonsingular,c(0,loggamma))
    
    estimate2 <- exp(lestimateMean)
    
    val <- (estimate1 + estimate2)/2	
    
  }
  val
}
#######################################################
#######################################################
# pair score function
# This function computes the percent inconsistant when clustering
# the replicates separately.
# 
# INPUTS
#  		data1: A matrix with samples in the rows and proteins in the columns
#			   to be clustered (replicate 1).  
# 		data2: A matrix similar to data1 that will also be clustered (replicate 2)
# 		distmet: The distance metric to be used.  The default is pearson correlation coefficient
# 		cmet: The linkage method.  Default is average linkage
#		k: Number of groups. Default is 8

pairscore <- function(data1,data2,distmet="pearson",cmet="average",k=8) {
  antibody2 <- dimnames(data2)[[2]]
  antibody1 <- dimnames(data1)[[2]]
  if (sum(antibody1 != antibody2) != 0) { stop("Names don't match") }
  
  hc1 <- hclust(distanceMatrix(data1, distmet), cmet)
  cutter1 <- cutree(hc1,k=k)		
  names(cutter1) <- antibody1
  sc1 <- matrix(NA,ncol=length(antibody1),nrow=length(antibody1))
  for (i in 1:(length(antibody1)-1)) {
    for (j in (i+1):length(antibody1)) {
      sc1[i,j] <- (cutter1[antibody1[i]]==cutter1[antibody1[j]])+0
    }
  }
  hc2 <- hclust(distanceMatrix(data2, distmet), cmet)
  cutter2 <- cutree(hc2,k=k)
  names(cutter2) <- antibody2
  sc2 <- matrix(NA,ncol=length(antibody2),nrow=length(antibody2))
  for (i in 1:(length(antibody2)-1)) {
    for (j in (i+1):length(antibody2)) {
      sc2[i,j] <- (cutter2[antibody2[i]]==cutter2[antibody2[j]])+0
    }
  }
  sc1 <- sc1[upper.tri(sc1)]
  sc2 <- sc2[upper.tri(sc2)]
  tmp <- abs(sc1-sc2)
  sum(tmp)/length(tmp)
}
