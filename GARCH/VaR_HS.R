data=read.csv("Desktop/q1.csv",header=T)
p1=data$Price.of.Stock.1[2:1001]
p2=data$Price.of.Stock.2[2:1001]
r1=data$Return.of.Stock.1[2:1001]
r2=data$Return.of.Stock.2[2:1001]

var=vector()
PL=1000*p1[1000]*(exp(r1)-1)+1000*p2[1000]*(exp(r2)-1)

-quantile(PL,0.05)

-mean(PL[PL<quantile(PL,0.05)])

