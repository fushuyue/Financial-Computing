setwd("C:/Users/think/Desktop")
data=read.csv("hw3 data.csv",header=T)

bs=function(which=c("call","put"),s0,k,div,r,t,vol){

	d1=(log(s0/k)+0.5*(r-div)*t)/(vol*sqrt(t))+0.5*vol*sqrt(t)
	d2=d1-vol*sqrt(t)
	if(which=="call")
		price=s0*exp(-div*t)*pnorm(d1)-k*exp(-r*t)*pnorm(d2)
	if(which=="put")
		price=k*exp(-r*t)*pnorm(-d2)-s0*exp(-div*t)*pnorm(-d1)
	price
}


EuropeanOptionImpliedVolatility=function(which=c("call","put"),price,s0,k,div,r,t,guess){

	if(which=="call")
	{
		c=bs("call",s0,k,div,r,t,guess)
		if(c<price)
		{
			for(i in 1:10)
			{
				guess2=guess+0.1*i
				c2=bs("call",s0,k,div,r,t,guess2)
				if(c2>price) break
			}

			for(j in 1:10000)
			{
				temp=(guess2+guess)/2
				if(abs((bs("call",s0,k,div,r,t,temp)-price))<0.000001) break
				if(bs("call",s0,k,div,r,t,temp)<price) guess=temp
				else guess2=temp			
			}
		}
		else{

			for(i in 1:10)
			{
				guess2=guess-0.1*i
				c2=bs("call",s0,k,div,r,t,guess2)
				if(c2<price) break
			}

			for(j in 1:10000)
			{
				temp=(guess2+guess)/2
				if(abs((bs("call",s0,k,div,r,t,temp)-price))<0.000001) break
				if(bs("call",s0,k,div,r,t,temp)>price) guess=temp
				else guess2=temp			
			}

		}
	}
	else{

		c=bs("put",s0,k,div,r,t,guess)
		if(c<price)
		{
			for(i in 1:10)
			{
				guess2=guess+0.1*i
				c2=bs("put",s0,k,div,r,t,guess2)
				if(c2>price) break
			}

			for(j in 1:10000)
			{
				temp=(guess2+guess)/2
				if(abs((bs("put",s0,k,div,r,t,temp)-price))<0.000001) break
				if(bs("put",s0,k,div,r,t,temp)<price) guess=temp
				else guess2=temp			
			}
		}
		else{

			for(i in 1:10)
			{
				guess2=guess-0.1*i
				c2=bs("put",s0,k,div,r,t,guess2)
				if(c2<price) break
			}

			for(j in 1:10000)
			{
				temp=(guess2+guess)/2
				if(abs((bs("put",s0,k,div,r,t,temp)-price))<0.000001) break
				if(bs("put",s0,k,div,r,t,temp)>price) guess=temp
				else guess2=temp			
			}

		}
	}

	(guess2+guess)/2
}


spx=data$SPX.Index;
indu=data$INDU.Index;
vix=data$VIX.Index;
vxd=data$VXD.Index;
v_spx_call=EuropeanOptionImpliedVolatility("call",((49.50+50.10)/2),spx[1],1865,0.0222,0.0025,(25/252),0.25)
v_spx_put=EuropeanOptionImpliedVolatility("put",((55.50+56.10)/2),spx[1],1865,0.0222,0.0025,(25/252),0.25)
v_djx_call=EuropeanOptionImpliedVolatility("call",((3.85+4.05)/2),(indu[1]/100),160,0.0248,0.0025,(25/252),0.25)
v_djx_put=EuropeanOptionImpliedVolatility("put",((4.70+4.90)/2),(indu[1]/100),160,0.0248,0.0025,(25/252),0.25)
v_spx_call
v_spx_put
v_djx_call
v_djx_put