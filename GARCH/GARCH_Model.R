setwd("C:/Users/think/Desktop")
data=read.csv("data.csv",header=T)

spx=data$LAST_PRICE[1:1001]
v2=0
spxreturn=rev(log(spx[-length(spx)]/spx[-1]))

#b[1]=sigma1 b[2]=sigma b[3]=alpha b[4]=beta
#likelihood function
garch11=function(R,b){
  v2[1]=b[1]*b[1];
  for(i in 2:length(R))
  v2[i]=(1-b[3]-b[4])*b[2]*b[2]+b[3]*R[i-1]*R[i-1]+b[4]*v2[i-1];
  return(log(v2)+R*R/v2+log(2*pi))
}

#quasi-Maximum likelihood estimation
initialvalue=c(sqrt(sum(spxreturn[1:20]*spxreturn[1:20])/20),sqrt(sum(spxreturn*spxreturn)/1000),0.1,0.8)
result=optim(initialvalue,fn=function(b){0.5*sum(garch11(spxreturn,b))}, method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"))

alpha=result$par[3]
beta=result$par[4]
sigma1=result$par[1]
sigma=result$par[2]

k=24

forecast=function(R,b){
	temp=0
	v2[1]=b[1]*b[1];
 	 for(i in 2:1001)
 	 v2[i]=(1-b[3]-b[4])*b[2]*b[2]+b[3]*R[i-1]*R[i-1]+b[4]*v2[i-1];

	sigma_f=v2[1001]
	for(i in 1:24){
	temp=temp+(b[3]+b[4])^(i-1)*(sigma_f-b[2]^2)
	}

	k*b[2]^2+temp
	attach(mtcars)
	par(mfrow=c(2,1))
	plot(v2,type='l',xaxt='n',yaxt='n',xlab='',ylab='',lwd=2,main='volatility')
	plot(R,type='l',xaxt='n',yaxt='n',xlab='',ylab='',col="grey",main='return')
	
}



b=c(sigma1,sigma,alpha,beta)

forecast(spxreturn,b)
#sqrt(forecast(spxreturn,b)*252/24)
