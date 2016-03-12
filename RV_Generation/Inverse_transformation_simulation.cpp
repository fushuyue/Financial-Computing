//
//  main.cpp
//  111
//
//  Created by fushuyue on 2/27/16.
//  Copyright © 2016 Monte Carlo Simulation. All rights reserved.
//

#include <iostream>
#include <cmath>

const double k=100;
const double t=(double)1/52;
const double s0=100;
const double q=0;
const double sigma=0.3;
const double r=0.005;
const double discont=exp(-r*t);
#define a0 2.50662823884
#define a1 -18.61500062529
#define a2 41.39119773534
#define a3 -25.44106049637
#define b0 -8.47351093090
#define b1 23.08336743743
#define b2 -21.06224101826
#define b3  3.13082909833
#define c0 0.3374754822726147
#define c1 0.9761690190917186
#define c2 0.1607979714918209
#define c3 0.0276438810333863
#define c4 0.0038405729373609
#define c5 0.0003951896511919
#define c6 0.0000321767881768
#define c7 0.0000002888167364
#define c8 0.0000003960315187
#define pi 3.141592653589793238

double get_uniform(){
	return (double)(random()/(((double)pow(2,31)-1)));
}

double min(double x,double y){
	if(x<y) return x;
	else return y;
	
}

double cdf(double x) {
	double v[16];double v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,c;
	v1 =1.253314137315500;
	v2 =0.6556795424187985;
	v3 =0.4213692292880545;
	v4 =0.3045902987101033;
	v5 =0.2366523829135607;
	v6 =0.1928081047153158;
	v7 =0.1623776608968675;
	v8 =0.1401041834530502;
	v9 =0.1231319632579329;
	v10 =0.1097872825783083;
	v11 = 0.09902859647173193;
	v12 = 0.09017567550106468;
	v13 = 0.08276628650136917;
	v14 = 0.0764757610162485;
	v15 =0.07106958053885211;
	c = 0.918938533204672;
	v[0]=0;
	v[1]=v1;v[2]=v2;v[3]=v3;v[4]=v4;v[5]=v5;v[6]=v6;v[7]=v7;v[15]=v15;
	v[8]=v8;v[9]=v9;v[10]=v10;v[11]=v11;v[12]=v12;v[13]=v13;v[14]=v14;
	double j,z,h,a,b,q,s,y,temp;
	temp=std::abs(x);
	j=min(temp+0.5,14);
	j=int(j);
	z=std::abs(j);
	h=std::abs(x)-z;
	a=0;
	a=v[int(j+1)];
	b=z*a-1;q=1;s=a+h*b;
	for(int ii=2;ii<=24-j;ii=ii*2)
	{
		a=(a+z*b)/ii;
		b=(b+z*a)/(ii+1);
		q=q*h*h;
		s=s+q*(a+b*h);
	}
	temp=-0.5*x*x-c;
	y=s*exp(temp);
	if(x>0) y=1-y;
	return y;
}


double inverse(){
	
	double u=get_uniform();
	
	double y=u-0.5;
	double x;
	if(std::abs(y)<0.42)
	{
		double r=y*y;
		x=y *(((a3 *r + a2) *r + a1) *r + a0)/((((b3*r+b2)*r+b1) *r+b0)*r+1);
	}
	else {double r=u;
	if(y>0)	r=1-u;
	r=log(-log(r));
	x=c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
	if(y<0) x=-x;
	
	 x=x+(u-cdf(x))/(exp(-0.5*x*x)/sqrt(2*pi));}
	return x;
	
}


double get_st()
{
	double temp;
	temp=s0*exp(t*(r-q-0.5*0.3*0.3)+0.3*sqrt(t)*inverse());
	
	return (double)temp;
}

double max(double a, double b){
	if(a>=b)
		return a;
	else
		return b;
}

double get_price()
{
	double st=get_st();
	if(st<=k) return discont*(-st+k);
	else return 0;
}

int main(int argc, const char * argv[]){
	
	double start_s=clock();
	long n;
	// number of times for simulation
	sscanf(argv[1],"%ld",&n);
	double *p=new double[n];
	for(long i=0;i<n;i++)
		p[i]=get_price();
	double sum=0;
	double sum1=0;
	for(long i=0;i<n;i++)
	{
		sum=sum+p[i];
		sum1=sum1+p[i]*p[i];
	}
	double se;
	se=sqrt((sum1-n*(sum/n)*(sum/n))/(n-1))/sqrt(n);
	std::cout<<n<<"	"<<sum/n<<"	"<<std::abs(sum/n-1.6547)<<"	"<<se;
	double stop_s=clock();
	double time=(stop_s-start_s)/double(CLOCKS_PER_SEC);
	std::cout<<"	"<<time<<std::endl;
	
}

