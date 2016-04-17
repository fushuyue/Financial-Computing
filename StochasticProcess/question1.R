# simulate gaussian rv with correlation 0.5 
simulate_x = function(){

	cov = matrix(c(1,0.5,0.5,1),nrow=2)
	A = t(chol(cov))
	e = matrix(rnorm(2),nrow=2)
	x = A %*% rnorm(2)
	x
}

simulate_path=function(){
	T = 5*365-1
	x1 = vector();x2 = vector()
	for(i in 1:T){
		x = simulate_x()
		x1 = c(x1,x[1])
		x2 = c(x2,x[2])
	}
	x = data.frame(x1,x2)
	x
}

Revenue = function(){

	T = 5*365
	delta_t = 1
	h = 0.4
	Y = 200*20
	n = 50
	x=simulate_path()
	x1 = x$x1
	x2 = x$x2
	k = 40000
	# initialize
	s = vector()
	r = vector()
	payoff = vector()
	PE = vector()
	PG = vector()
	PE[1] = 20
	PG[1] = 24
	s[1] = 20-h*24
	r[1] = s[1]*Y
	payoff[1] = max((k-r[1]),0)
	for(i in 1:(T-1)){
			
		PE[i+1] = PE[i]+0.75*(32-PE[i])*delta_t+0.5*PE[i]*x1[i]
		PG[i+1] = PG[i]+0.15*(40-PG[i])*delta_t+0.6*PG[i]*x2[i]
		s[i+1] = PE[i+1]-h*PG[i+1]
		r[i+1] = max((s[i+1]*Y),0)
		payoff[i+1] = max((k-r[i+1]),0)
	}

	df = data.frame(s,r,payoff)
	df
}

keke = function(n){
	total = vector()
	not_operate = vector()
	less_than_40K = vector()
	payoff = vector()
	for(i in 1:n){
		df = Revenue()
		if(i%%100==0)
		print(c('i=',i))
		total[i] = sum(df$r)
		not_operate[i] = length(which(df$s<0))
		less_than_40K[i] = length(which(df$r<40000))
		payoff[i] = mean(df$payoff)
	}

	df = data.frame(total,not_operate,less_than_40K,payoff)
	df
}

df = keke(5000)
hist = qplot(x=df$total,data=df,geom='histogram',main='the histogram of the total revenue',xlab='total revenue',fill='firebrick',bins=100,xlim=c(1e+08,3.5e+08))
density = qplot(x=df$total,data=df,geom='density',color='firebircks1',main='the density of the total revenue',xlab='total revenue',xlim=c(1e+08,3.5e+08))
multiplot(hist,density)






