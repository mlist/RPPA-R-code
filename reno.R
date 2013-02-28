#############################################################################################################
############################################ Library ########################################################
#############################################################################################################

library(cobs) # require package "SparseM"
library(Hmisc)
library(quantreg)
library(msm) # This library is only used for generating simulation data

#############################################################################################################
############################################# General function ##############################################
#############################################################################################################

lp = function(p) {logb(p/(1-p), 2)}
ea = function(a) {2^(a)/(1+2^(a))}
foo = function(x, A, B, alpha=0, beta=10/(B-A)) {
  A + (B-A)*ea(alpha+beta*x)
}

#############################################################################################################
################################################### Tabus ###################################################
#############################################################################################################

#### Tabus method
#### Hu. et al edition

dilutionFitrep = function(data, nrep, verbose=FALSE, trace=FALSE) {
  
  # We assume the data is in the form of a matrix, with each row
  # representing a dilution series (increasing concentrations)
  # for the same protein in a different sample. Try to fit a
  # common logistic curve, placing all series appropriately
  # along it.
  
  nsamp = nrow(data)
  ndilut = ncol(data)
  balance = (1:ndilut)-(1+ndilut)/2
  
  # make a first guess at the extreme levels we might see
  mini = apply(data, 1, min)
  maxi = apply(data, 1, max)
  A = min(mini)
  B = max(maxi)
  
  # transform the data to where it should be linear
  lindata = lp((data-A)/(B-A))
  
  # initial estimate of the logistic slope across the dilution steps
  beta = mean(apply(lindata, 1, function(x){
    y = x[!is.na(x) & !is.infinite(x)]
    max(y)-min(y)
  }))/(ndilut-1)
  
  # for each sample, estimate the offset
  passer = apply(lindata, 1, function(x, beta, balance) {
    median(x/beta - balance)
  }, beta, balance)
  
  
  passermat=matrix(passer,ncol=nrep,byrow=T)
  passeravg=apply(passermat,1,mean)
  passer=rep(passeravg,each=nrep)
  
  
  # put our current guess at the x and y values into vectors
  
  xval = rep(balance, nsamp) + rep(passer, each=ndilut)
  yval = as.vector(t(data))
  ox = order(xval)
  A2 = median(yval[ox][1:100])
  B2 = median(rev(yval[ox])[1:100])
  
  dum = data.frame(xval=xval, yval=yval)
  nls.model = nls(yval ~ alpha + beta*ea(gamma*xval), data=dum,
                  start=list(alpha=A, beta=B-A, gamma=beta),
                  control=nls.control(maxiter=10000,minFactor=1/2^50))
  
  cf = coef(nls.model)
  p.alpha = cf['alpha']
  p.beta  = cf['beta']
  p.gamma = cf['gamma']
  
  pass2 = rep(NA, (nrow(data)/nrep))
  warn2  = rep('', (nrow(data)/nrep))
  names(warn2) = (dimnames(data)[[1]])[seq(1,nrow(data),by=nrep)]
  
  for (i in 1:(nrow(data)/nrep)) {
    if(verbose) print(i)
    dum = data.frame(Y=as.vector(data[(nrep*(i-1)+1):(nrep*i),]),B=rep(balance,each=nrep))
    dum$p.alpha = p.alpha
    dum$p.beta = p.beta
    dum$p.gamma = p.gamma
    
    tmp = try(nls(Y ~ p.alpha + p.beta*ea(p.gamma*(B+X)), data=dum,
                  start=list(X=passer[nrep*(i-1)+1]),trace=trace,control=nls.control(maxiter=10000,minFactor=1/2^50)))
    
    if (is(tmp, 'try-error')) {
      warning('nls error')
      warn2[i] = paste(warn2[i], 'nls-error', collapse=' ')
      tmp = try(nls(Y ~ p.alpha + p.beta*ea(p.gamma*(B+X)), data=dum,
                    start=list(X=passer[nrep*(i-1)+1]),trace=trace, algorithm='plinear',
                    control=nls.control(maxiter=10000,minFactor=1/2^50)))
      
      if (is(tmp, 'try-error')) {
        warning('unavoidable nls error')
        warn2[i] = paste(warn2[i], 'unavoidable nls-error')
        pass2[i] = passer[nrep*(i-1)+1]
      } else {
        pass2[i] = coef(tmp)['X']
      }
    } else {
      pass2[i] = coef(tmp)
    }
  }
  
  list(A=A, B=B, A2=A2, B2=B2, passer=passer, pass2=pass2,
       p.alpha=p.alpha, p.beta=p.beta, p.gamma=p.gamma, balance=balance,
       warn=warn2, fit=nls.model)
}

#############################################################################################################
####################################### Nonparametric method ################################################
#############################################################################################################

##########functions to estimate x using the logistic approach#########

lp = function(p) {
  p=ifelse(p<0.0000001,0.0000001,p)
  logb(p/(1-p), 2)}

ea = function(a) {2^(a)/(1+2^(a))}

foo = function(x, A, B, alpha=0, beta=10/(B-A)) {
  A + (B-A)*ea(alpha+beta*x)}

############functions to estimate x using the nonparametric approach############

fspline=function(xvec, aknot, acoef){
  
  nknot=length(aknot)
  aknotnew=c(aknot[1], aknot[1],aknot, aknot[nknot], aknot[nknot])
  ncoef=length(acoef)
  
  xvec[xvec<(aknot[1]+(1e-8))]=aknot[1]+(1e-8)
  xvec[xvec > (aknot[nknot]-(1e-8))]=aknot[nknot]-(1e-8)
  
  if (ncoef == nknot + 2)
  {acoefnew=acoef[-ncoef]
   a=spline.des(aknotnew, xvec,ord=3)
   fvalvec= (a$design)%*%acoefnew
  }
  
  return(fvalvec)
}

getfityrep=function(x,x0mat,fity,nrep,data){
  fityrep=fity[seq(1,length(fity),by=nrep)]
  x0matsub=x0mat[seq(1,dim(data)[1],by=nrep),]
  x0vec=as.vector(x0matsub)
  
  if (x<min(x0vec)) youtput=fityrep[x0vec==min(x0vec)]
  if (x>max(x0vec)) youtput=fityrep[x0vec==max(x0vec)]
  
  if (sum(x0vec==x)>0) {
    youtput=fityrep[x0vec==x]
    if (length(youtput)>1) youtput=youtput[1]     }
  if (sum(x0vec==x)==0 && x>min(x0vec) && x<max(x0vec)){
    
    x0vecsub=x0vec[(x0vec-x)>0]
    fitysub=fityrep[(x0vec-x)>0]
    xright=x0vecsub[(x0vecsub-x)==min(x0vecsub-x)]
    yright=fitysub[(x0vecsub-x)==min(x0vecsub-x)]
    
    if (length(xright)>1) {xright=xright[1]; yright=yright[1]}
    x0vecsub=x0vec[(x0vec-x)<0]
    fitysub=fityrep[(x0vec-x)<0]
    xleft=x0vecsub[(x0vecsub-x)==max(x0vecsub-x)]
    yleft=fitysub[(x0vecsub-x)==max(x0vecsub-x)]
    
    if (length(xleft)>1) {xleft=xleft[1]; yleft=yleft[1]}
    
    if (length(xright)>0 && length(xleft)>0){
      lineslope=(yright-yleft)/(xright-xleft)
      lineint=yleft-lineslope*xleft
      youtput=lineint+lineslope*x
    }
  }
  
  if (length(youtput)>1) youtput=youtput[1]
  return(youtput)
}

getnewxrep=function(nrep,fity,xi,datamat,gridsize,windowsize,miny,disminyy,lowxwt,data,x0mat,balance){
  
  datavec=as.vector(datamat)
  allwt=rep(1,length(datavec))
  if (sum((datavec-miny)<=disminyy)>0){
    allwt[(datavec-miny)<=disminyy]=lowxwt
  }
  
  x0vec=as.vector(x0mat[seq(1,dim(data)[1],by=nrep),])
  xicandvec=seq(max(xi-windowsize,sort(x0vec)[6]+1),min(xi+windowsize,-sort(-x0vec)[6]-1),length=gridsize)
  l1dis=NULL
  
  for (j in 1:length(xicandvec)){
    x0dilute=xicandvec[j]+balance
    fitydilute=NULL
    for (k in 1:length(balance)){
      outputk=getfityrep(x0dilute[k],x0mat,fity,nrep,data)
      fitydilute=c(fitydilute,outputk)
    }
    
    fitydilute=rep(fitydilute,each=nrep)
    l1dis=c(l1dis, sum(abs(datavec-fitydilute)*allwt) )
  }
  
  
  output=xicandvec[l1dis==min(l1dis,na.rm=T)]
  len.out=length(output)
  
  if (len.out>1) { loc=as.integer(len.out/2)
                   output=output[loc] }
  minl1=min(l1dis,na.rm=T)
  return(list(output=output,minl1=minl1))
}

getnonpest=function(data,nrep,randW=FALSE){
  
  #####the first step#####
  
  ###obtain the initial value x###
  
  nsamp = nrow(data)
  ndilut = ncol(data)
  balance = (1:ndilut)-(1+ndilut)/2
  
  # make a first guess at the extreme levels we might see
  mini = apply(data, 1, min)
  maxi = apply(data, 1, max)
  A = min(mini)
  B = max(maxi)
  
  # transform the data to where it should be linear
  
  lindata = lp((data-A)/(B-A))
  
  # initial estimate of the logistic slope across the dilution steps
  
  beta = mean(apply(lindata, 1, function(x){
    y = x[!is.na(x) & !is.infinite(x)]
    max(y)-min(y)
  }))/(ndilut-1)
  
  
  # for each sample, estimate the offset
  
  passer = apply(lindata, 1, function(x, beta, balance) {
    median(x/beta - balance)
  }, beta, balance)
  
  passermat=matrix(passer,ncol=nrep,byrow=T)
  passeravg=apply(passermat,1,mean)
  passeravg=passeravg-median(passeravg)
  passeravg=rep(passeravg,each=nrep)
  
  ###the second step###
  ###doing 5 iterations###
  x0new=passeravg
  
  if(randW){weights=rexp(length(as.vector(data)))}else
  {weights=rep(1,length(as.vector(data)))}
  for (ii in 1:5){
    
    numsearchpt=ii*20
    x0mat=NULL
    for (i in 1:nrow(data)) x0mat=rbind(x0mat,balance+x0new[i])
    
    if (ii==1) {a=cobs(as.vector(x0mat),as.vector(data),w=weights, nknots=40,constraint="increase",lambda=-1,print.warn = FALSE, print.mesg = FALSE)
                lambda1=a$lambda}
    if (ii >1)   a=cobs(as.vector(x0mat),as.vector(data),w=weights, nknots=40,constraint="increase",lambda=lambda1,print.warn = FALSE, print.mesg = FALSE)
    fity=fitted(a)
    
    aknot=a$knots
    acoef=a$coef
    
    if (ii<5) windowsize=2
    if (ii==5) windowsize=1
    
    x0old=x0new
    x0new=NULL
    
    for (i in 1:(dim(data)[1]/nrep)){
      Y = as.vector(data[(nrep*(i-1)+1):(nrep*i),])
      Steps = rep(balance,each=nrep)
      
      nlrqres=nlrq(Y ~ fspline(X+Steps,aknot,acoef),start = list(X = x0old[(nrep*(i-1)+1)]), trace = FALSE, control=nlrq.control(maxiter=20,eps=1e-03))
      
      tmp = try(nlrqres, TRUE)
      if (length(coef(tmp))==0) x0new=c(x0new, x0old[(nrep*(i-1)+1)])
      if (length(coef(tmp))>0 && abs(coef(tmp)-x0old[(nrep*(i-1)+1)])>2){
        temp=getnewxrep(nrep,fity,x0old[nrep*(i-1)+1],data[(nrep*(i-1)+1):(nrep*i),],numsearchpt,windowsize,miny=500,disminyy=2000,lowxwt=1,data,x0mat,balance)
        x0new=c(x0new,temp$output)
      }
      if (length(coef(tmp))>0 && abs(coef(tmp)-x0old[(nrep*(i-1)+1)])<=2) x0new = c(x0new, coef(tmp))
    }
    x0new=x0new-median(x0new)
    x0new1=x0new
    x0new=rep(x0new,each=nrep)
    
  }
  x0mat=NULL
  for (i in 1:nrow(data)) x0mat=rbind(x0mat,balance+x0new[i])
  a=cobs(as.vector(x0mat),as.vector(data),nknots=40,constraint="increase",lambda=lambda1,print.warn = FALSE, print.mesg = FALSE)
  return(list(x0new=x0new1,lambda1=lambda1,fitted=matrix(fity,nsamp,ndilut),residuals=matrix(data-fity,nsamp,ndilut),fit=a))
}


#############################################################################################################
###################################### Spline Coefficients ##################################################
#############################################################################################################

# To obtain the coefficients of an existing spline in a matrix format.
# Each row represents the coefficents in an certain interval in an decreasing order, i.e. "x^2", "x", "1".

coef_for_spline<-function(fit){
  
  knot=fit$knots;
  
  coefficient<-matrix(0,(length(knot)-1),3)
  colnames(coefficient)=c("x^2","x","1")
  
  for(ii in 1:(length(knot)-1)){
    
    x=(3*knot[ii]+knot[ii+1])/4;
    y=(2*knot[ii]+2*knot[ii+1])/4;
    z=(knot[ii]+3*knot[ii+1])/4;
    
    fx=predict(fit,x)[,2];
    fy=predict(fit,y)[,2];
    fz=predict(fit,z)[,2];
    
    a=-(fy*x - fz*x - fx*y + fz*y + fx*z - fy*z)/((-x + y)*(y - z)*(-x + z))
    b=-(-fy*x^2 + fz*x^2 + fx*y^2 - fz*y^2 - fx*z^2 +
          fy*z^2)/((x - y)*(x - z)*(y - z))
    d=-(-fz*x^2*y + fz*x*y^2 + fy*x^2*z - fx*y^2*z - fy*x*z^2 +
          fx*y*z^2)/((y - z)*(x^2 - x*y - x*z + y*z))
    
    coefficient[ii,]=t(t(c(a,b,d)))
  }
  
  return(coefficient)
}

#############################################################################################################
########################################### Derivatives #####################################################
#############################################################################################################

# To obtain the derivative of an existing spline at one particular point.

deriv_for_spline<-function(xx,fit,coefficient){
  
  knot=fit$knots;
  xx=t(t(xx));
  
  xx.tf=t(apply(xx, 1, function(x,knot){x>knot}, knot));
  xx.interval=apply(xx.tf,1,sum);
  
  deri=2*coefficient[xx.interval,1]*xx+coefficient[xx.interval,2]
  
  return(as.vector(deri))
}

#############################################################################################################
############################################ Prediction #####################################################
#############################################################################################################

# To obtain the predictation of an existing spline at the boundaries in a matrix format.

pred_for_spline<-function(fit){
  
  knot=fit$knots
  rank=1:length(knot)
  
  e=10^(-3)
  
  prediction=predict(fit,knot[1:length(knot)-1])
  prediction=rbind(prediction,predict(fit,knot[length(knot)]-e))
  prediction[,2]=round(prediction[,2],3)
  
  # removing the flat region.
  # lowerbound will denote the start interval, after removing the flat region.
  
  lowerbound=max(rank[prediction[,2]==min(prediction[,2])])
  upperbound=min(rank[prediction[,2]==max(prediction[,2])])
  
  prediction=prediction[lowerbound:upperbound,]
  
  list(prediction=prediction,lowerbound=lowerbound)
}

#############################################################################################################
############################################ Inverse ########################################################
#############################################################################################################

# To obtain the inverse of a given spline at any y. If y is outside the range of spline, we shrink it back
# to the range.

inverse_for_spline<-function(yy,fit,coefficient,prediction,lowerbound){
  
  e=10^(-3)
  yy=round(yy,3)
  
  yy[yy<=prediction[1,2]]=prediction[1,2]+e
  yy[yy>=(prediction[dim(prediction)[1],2])]=(prediction[dim(prediction)[1],2])-e
  yy=t(t(yy))
  
  yy.tf=t(apply(yy, 1, function(x,knot){x>prediction[,2]}, prediction));
  yy.interval=apply(yy.tf,1,sum);
  
  yy[yy==prediction[yy.interval,2]]=yy[yy==prediction[yy.interval,2]]+e;
  
  yy.interval=yy.interval+lowerbound-1;
  
  a=coefficient[yy.interval,1]
  b=coefficient[yy.interval,2]
  c=coefficient[yy.interval,3]
  
  # Here, we have to take the sign into consideration.
  # If a>0, we need to pick up the larger one to make the whole curve incresing
  # If a<0, we need to pick up the samller one to make the curve increasing
  # In general, we have to pick up the one with -b+sqrt(~)
  inve<-(-b+sqrt(apply(b^2-4*a*c+4*a*yy,1, function(x) {max(0,x)})))/2/a
  return(as.vector(inve))
}

#############################################################################################################
######################################## Initial Estimates ##################################################
#############################################################################################################

# To obtain the initial estimates

initial_x<-function(data,nrep=1){
  
  ####### First Guess #######
  ####### Treat them individually #######
  
  nsamp = nrow(data)
  ndilut = ncol(data)
  ntotal = nsamp*ndilut
  balance = (1:ndilut)-(1+ndilut)/2
  
  mini = apply(data, 1, min)
  maxi = apply(data, 1, max)
  A = min(mini)
  B = max(maxi)
  
  lindata = lp((data-A)/(B-A))
  
  beta = mean(apply(lindata, 1, function(x){
    y = x[!is.na(x) & !is.infinite(x)]
    max(y)-min(y)
  }))/(ndilut-1)
  
  passer = apply(lindata, 1, function(x, beta, balance) {
    median(x/beta - balance)
  }, beta, balance)
  
  passermat=matrix(passer,ncol=nrep,byrow=T)
  passeravg=apply(passermat,1,mean)
  passeravg=passeravg-median(passeravg)
  
  x.new=rep(passeravg,each=nrep)
  x.hat<-matrix(rep(passeravg,each=ndilut*nrep)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
  initial<-x.new
  initial=matrix(initial,nrep,(nsamp/nrep))
  
  list(initial=initial,x.hat=x.hat,x.new=x.new)
}

#############################################################################################################
################################### Update estimates within loops ###########################################
#############################################################################################################

# For lambda2=infinity

update_x_infinity<-function(xi,yi,fit,coefficient,prediction,lowerbound,balance,nrep,x.old,i,fity,data,numsearchpt,windowsize,x.hat){
  
  weights<-abs(deriv_for_spline(xi,fit,coefficient));
  value<-inverse_for_spline(yi,fit,coefficient,prediction,lowerbound);
  value<-value-rep(balance,nrep);
  
  if(max(round(weights,4))!=0){
    rqfit<-rq(value~1,weights=weights)
    c1<-rqfit$coefficients[1]
    x.update=c(rep(c1,nrep))
  } else  {   x.update=rep(median(xi),nrep) }
  
  if( max(abs(x.update-x.old[(nrep*(i-1)+1:nrep)]))>2 )
  {
    temp=getnewxrep(nrep,fity,x.old[nrep*(i-1)+1],data[(nrep*(i-1)+1):(nrep*i),],numsearchpt,windowsize,miny=500,disminyy=2000,lowxwt=1,data,x.hat,balance)
    x.update=rep(temp$output,nrep)
  }
  return(x.update)
}

# For finite lam2

update_x_lam2<-function(xi,yi,fit,XX,coefficient,prediction,lowerbound,balance,nrep,lambda2,x.old,i,fity,data,numsearchpt,windowsize,x.hat,ndilut){
  
  weights<-abs(deriv_for_spline(xi,fit,coefficient));
  value<-inverse_for_spline(yi,fit,coefficient,prediction,lowerbound);
  value<-value-rep(balance,nrep);
  
  if(max(round(weights,4))!=0){
    yy<-c(value,0,0,0)
    weights<-c(weights,lambda2,lambda2,lambda2)
    rqfit<-rq(yy~XX-1,weights=weights)
    x.update=c(rqfit$coefficients[1],rqfit$coefficients[2],rqfit$coefficients[3])
  } else  {   x.update=rep(median(xi),nrep)}
  
  if( max(abs(x.update-x.old[(nrep*(i-1)+1:nrep)]))>2 )
  {
    temp1=getnewxrep(nrep,fity,x.old[nrep*(i-1)+1],data[(nrep*(i-1)+1):(nrep*i),],numsearchpt,windowsize,miny=500,disminyy=2000,lowxwt=1,data,x.hat,balance)
    x_nrep3=rep(temp1$output,nrep)
    x_nrep1=NULL;
    
    for(ww in 1:nrep)
    {
      temp2=getnewxrep(nrep=1,fity,x.old[nrep*(i-1)+ww],data[nrep*(i-1)+ww,],numsearchpt,windowsize,miny=500,disminyy=2000,lowxwt=1,data,x.hat,balance)
      x_nrep1=c(x_nrep1,temp2$output)
    }
    
    loss_nrep3=sum(weights[1:(ndilut*nrep)]*abs((rep(x_nrep3,each=ndilut)-value)))
    loss_nrep1=sum(weights[1:(ndilut*nrep)]*abs((rep(x_nrep1,each=ndilut)-value)))+lambda2*(abs(x_nrep1[1]-x_nrep1[2])+abs(x_nrep1[1]-x_nrep1[3])+abs(x_nrep1[3]-x_nrep1[2]))
    
    if(loss_nrep3<loss_nrep1) x.update=x_nrep3 else {x.update=x_nrep1}
  }
  return(x.update)
}

#############################################################################################################
#################################### choose x from different iterations #####################################
#############################################################################################################

# After certain number of iterations, we track the whole iteration results and search for the best one in the
# sense of minimizing the prediction error.

final_x<-function(track.xi,yi,fit,ndilut,nrep,balance){
  
  Object.f<-NULL;
  
  for (ii in 1:4) {
    xi = matrix(rep(track.xi[,ii],ndilut),nrep,ndilut) + matrix(rep(balance,each=nrep),nrep,ndilut)
    xi = as.vector(xi)
    yi = as.vector(yi)
    Object.f = c(Object.f,sum(abs(track.xi[c(1,1,2),ii]-track.xi[c(2,3,3),ii])) + sum(abs(predict(fit,xi)[,2]-yi)))
  }
  
  pickid=order(Object.f)[1]
  return(track.xi[,pickid])
}

#############################################################################################################
################################################# Functions #################################################
#############################################################################################################

lp = function(p) {
  p=ifelse(p<0.0000001,0.0000001,p)
  logb(p/(1-p), 2)}
ea = function(a) {2^(a)/(1+2^(a))}
foo = function(x, A, B, alpha=0, beta=10/(B-A)) {
  A + (B-A)*ea(alpha+beta*x)
}

getfityrep=function(x,x0mat,fity,nrep,data){
  fityrep=fity[seq(1,length(fity),by=nrep)]
  x0matsub=x0mat[seq(1,dim(data)[1],by=nrep),]
  x0vec=as.vector(x0matsub)
  
  if (x<min(x0vec)) youtput=fityrep[x0vec==min(x0vec)]
  if (x>max(x0vec)) youtput=fityrep[x0vec==max(x0vec)]
  
  if (sum(x0vec==x)>0) {
    youtput=fityrep[x0vec==x]
    if (length(youtput)>1) youtput=youtput[1]     }
  if (sum(x0vec==x)==0 && x>min(x0vec) && x<max(x0vec)){
    
    x0vecsub=x0vec[(x0vec-x)>0]
    fitysub=fityrep[(x0vec-x)>0]
    xright=x0vecsub[(x0vecsub-x)==min(x0vecsub-x)]
    yright=fitysub[(x0vecsub-x)==min(x0vecsub-x)]
    
    if (length(xright)>1) {xright=xright[1]; yright=yright[1]}
    x0vecsub=x0vec[(x0vec-x)<0]
    fitysub=fityrep[(x0vec-x)<0]
    xleft=x0vecsub[(x0vecsub-x)==max(x0vecsub-x)]
    yleft=fitysub[(x0vecsub-x)==max(x0vecsub-x)]
    
    if (length(xleft)>1) {xleft=xleft[1]; yleft=yleft[1]}
    
    if (length(xright)>0 && length(xleft)>0){
      lineslope=(yright-yleft)/(xright-xleft)
      lineint=yleft-lineslope*xleft
      youtput=lineint+lineslope*x
    }
  }
  
  if (length(youtput)>1) youtput=youtput[1]
  return(youtput)
}

getnewxrep=function(nrep,fity,xi,datamat,gridsize,windowsize,miny,disminyy,lowxwt,data,x0mat,balance){
  
  datavec=as.vector(datamat)
  allwt=rep(1,length(datavec))
  if (sum((datavec-miny)<=disminyy)>0){
    allwt[(datavec-miny)<=disminyy]=lowxwt
  }
  
  x0vec=as.vector(x0mat[seq(1,dim(data)[1],by=nrep),])
  xicandvec=seq(max(xi-windowsize,sort(x0vec)[6]+1),min(xi+windowsize,-sort(-x0vec)[6]-1),length=gridsize)
  l1dis=NULL
  
  for (j in 1:length(xicandvec)){
    x0dilute=xicandvec[j]+balance
    fitydilute=NULL
    for (k in 1:length(balance)){
      outputk=getfityrep(x0dilute[k],x0mat,fity,nrep,data)
      fitydilute=c(fitydilute,outputk)
    }
    
    fitydilute=rep(fitydilute,each=nrep)
    l1dis=c(l1dis, sum(abs(datavec-fitydilute)*allwt) )
  }
  
  
  output=xicandvec[l1dis==min(l1dis,na.rm=T)]
  len.out=length(output)
  
  if (len.out>1) { loc=as.integer(len.out/2)
                   output=output[loc] }
  minl1=min(l1dis,na.rm=T)
  return(list(output=output,minl1=minl1))
}

##############################################################################################################
######################################### Full function for RERO #############################################
##############################################################################################################

# The whole procedure was used to deal with the situation that there are three replicates
# In order to pick up a reasonable estimates, we fix step.max=8.
# The changes after 4 or 5 iterations are tiny.

getnonpapprox=function(data,nrep=3,lambda2="infinity",CVindex="NULL",step.max=8){
  
  nsamp = nrow(data);
  ndilut = ncol(data);
  ntotal = nsamp*ndilut;
  balance = (1:ndilut)-(1+ndilut)/2;
  step = 1;
  
  x.ini = initial_x(data,nrep=3)
  x.hat = x.ini$x.hat
  x.new = x.ini$x.new
  
  Track.x=matrix(NA,nsamp,4) # we track the last 4 steps.
  
  #### lam2=Inifity
  if(lambda2=="infinity") {
    while(step<=step.max)
    {
      numsearchpt=step*20
      
      if (step<step.max) windowsize=2
      if (step==step.max) windowsize=1
      
      if(step==1) {fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=-1,print.warn = FALSE, print.mesg = FALSE);lambda=fit$lambda;}
      if(step!=1) {fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE);}
      fity=fitted(fit);
      
      x.old=x.new;
      x.new<-NULL;
      
      coefficient=coef_for_spline(fit);
      junk.pred=pred_for_spline(fit);
      prediction=junk.pred$prediction;
      lowerbound=junk.pred$lowerbound;
      
      for(i in 1:(nsamp/nrep))
      {
        xi<-as.vector(t(x.hat[nrep*(i-1)+1:nrep,]));   # x_j+d_l are read
        yi<-as.vector(t(data[nrep*(i-1)+1:nrep,]));
        x.new<-c(x.new,update_x_infinity(xi,yi,fit,coefficient,prediction,lowerbound,balance,nrep,x.old,i,fity,data,numsearchpt,windowsize,x.hat))
      }
      
      x.new=x.new-median(x.new)
      x.hat<-matrix(rep(x.new,each=ndilut)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
      if(step>=step.max-3) Track.x[,(4+step-step.max)]=x.new
      step=step+1
    }
    
    fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE)
    
    x.final<-NULL;
    for(jj in 1:(nsamp/nrep)){
      track.xi=Track.x[nrep*(jj-1)+1:nrep,]
      yi=data[nrep*(jj-1)+1:nrep,]
      x.final<-c(x.final,final_x(track.xi,yi,fit,ndilut,nrep,balance))
    }
    x.final<-x.final-median(x.final)
    x.hat.final<-matrix(rep(x.final,each=ndilut)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
    fit=cobs(as.vector(x.hat.final),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE)
    
  } else
    
    #### lam2: Finite
  {
    # Define the design matrix "XX" first
    x1<-c(rep(1,ndilut),rep(0,2*ndilut));
    x2<-c(rep(0,ndilut),rep(1,ndilut),rep(0,ndilut));
    x3<-c(rep(0,2*ndilut),rep(1,ndilut));
    X<-cbind(x1,x2,x3);
    XX=rbind(X,c(1,-1,0),c(0,1,-1),c(1,0,-1));
    
    if(CVindex=="NULL"){
      while(step<=step.max)
      {
        numsearchpt=step*20
        
        if (step<step.max) windowsize=2
        if (step==step.max) windowsize=1
        
        if(step==1) {fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=-1,print.warn = FALSE, print.mesg = FALSE);lambda=fit$lambda;}
        if(step!=1) {fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE);}
        fity=fitted(fit);
        
        x.old=x.new;
        x.new<-NULL;
        
        coefficient=coef_for_spline(fit);
        junk.pred=pred_for_spline(fit);
        prediction=junk.pred$prediction;
        lowerbound=junk.pred$lowerbound;
        
        for(i in 1:(nsamp/nrep))
        {
          xi<-as.vector(t(x.hat[nrep*(i-1)+1:nrep,]));   # x_j+d_l are read
          yi<-as.vector(t(data[nrep*(i-1)+1:nrep,]));
          x.new<-c(x.new,update_x_lam2(xi,yi,fit,XX,coefficient,prediction,lowerbound,balance,nrep,lambda2,x.old,i,fity,data,numsearchpt,windowsize,x.hat,ndilut))
        }
        
        x.new=x.new-median(x.new)
        x.hat<-matrix(rep(x.new,each=ndilut)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
        if(step>=step.max-3) Track.x[,(4+step-step.max)]=x.new
        step=step+1
        
      }
      
      fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE)
      x.final<-NULL;
      for(jj in 1:(nsamp/nrep)){
        track.xi=Track.x[nrep*(jj-1)+1:nrep,]
        yi=data[nrep*(jj-1)+1:nrep,]
        x.final<-c(x.final,final_x(track.xi,yi,fit,ndilut,nrep,balance))
      }
      x.final<-x.final-median(x.final)
      x.hat.final<-matrix(rep(x.final,each=ndilut)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
      fit=cobs(as.vector(x.hat.final),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE)
      
    } else
    {
      #### lam2: Cross validation
      
      ndilut=ndilut-1;
      data_rec<-data
      x.hat_rec<-x.hat
      balance_rec<-balance
      data=data[,-CVindex]
      x.hat=x.hat[,-CVindex]
      balance=balance[-CVindex]
      
      x1<-c(rep(1,ndilut),rep(0,2*ndilut));
      x2<-c(rep(0,ndilut),rep(1,ndilut),rep(0,ndilut));
      x3<-c(rep(0,2*ndilut),rep(1,ndilut));
      X<-cbind(x1,x2,x3);
      XX=rbind(X,c(1,-1,0),c(0,1,-1),c(1,0,-1));
      
      while(step<=step.max)
      {
        numsearchpt=step*20
        
        if (step<step.max) windowsize=2
        if (step==step.max) windowsize=1
        
        if(step==1) {fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=-1,print.warn = FALSE, print.mesg = FALSE);lambda=fit$lambda;}
        if(step!=1) {fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE);}
        fity=fitted(fit);
        
        x.old=x.new;
        x.new<-NULL;
        
        coefficient=coef_for_spline(fit);
        junk.pred=pred_for_spline(fit);
        prediction=junk.pred$prediction;
        lowerbound=junk.pred$lowerbound;
        
        for(i in 1:(nsamp/nrep))
        {
          xi<-as.vector(t(x.hat[nrep*(i-1)+1:nrep,]));   # x_j+d_l are read
          yi<-as.vector(t(data[nrep*(i-1)+1:nrep,]));
          x.new<-c(x.new,update_x_lam2(xi,yi,fit,XX,coefficient,prediction,lowerbound,balance,nrep,lambda2,x.old,i,fity,data,numsearchpt,windowsize,x.hat,ndilut))
        }
        
        x.new=x.new-median(x.new)
        x.hat<-matrix(rep(x.new,each=ndilut)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
        if(step>=step.max-3) Track.x[,(4+step-step.max)]=x.new
        step=step+1
        
      }
      
      fit=cobs(as.vector(x.hat),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE)
      x.final<-NULL;
      for(jj in 1:(nsamp/nrep)){
        track.xi=Track.x[nrep*(jj-1)+1:nrep,]
        yi=data[nrep*(jj-1)+1:nrep,]
        x.final<-c(x.final,final_x(track.xi,yi,fit,ndilut,nrep,balance))
      }
      x.final<-x.final-median(x.final)
      x.hat.final<-matrix(rep(x.final,each=ndilut)+rep(balance,nsamp),nrow=nsamp,ncol=ndilut,byrow=T)
      fit=cobs(as.vector(x.hat.final),as.vector(data),nknots=40,constraint='increase',lambda=lambda,print.warn = FALSE, print.mesg = FALSE)
      
      fity=fitted(fit)
      est<-NULL
      for(kk in 1:nsamp)
        est<-c(est,getfityrep(x.final[kk]+balance_rec[CVindex],x.hat.final,fity,nrep,data))
    }}
  
  if(CVindex=="NULL") list(fit=fit,x.new=round(x.final,4)) else
    list(fit=fit,x.new=round(x.final,4),est=est)
}

# The end