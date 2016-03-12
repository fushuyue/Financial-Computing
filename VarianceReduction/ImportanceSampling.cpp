#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

#define m1 2147483647
#define m2 2145483479
#define a12 63308
#define a13 -183326
#define a21 86098
#define a23 -539608
#define q12 33921
#define q13 11714
#define q21 24919
#define q23 3976
#define r12 12979
#define r13 2883
#define r21 7417
#define r23 2071
#define Invmp1 4.656612873077393e-10;
const double k=1900;

const double t=(double)1/12;
const double s0=2000;
const double q=0.02;
const double sigma=0.3;
const double r=0.005;
const double discont=exp(-r*t);

int x10, x11, x12, x20, x21, x22;

//seeds x10, x11, x12 should be between 1 and 2147483646, x20, x21, x22 should be between 1 and 2145483478

int Random()
{
	int h,p12,p13,p21,p23;
	h = x10/q13;
	p13 = -a13*(x10-h*q13)-h*r13;
	h = x11/q12; p12 = a12*(x11-h*q12)-h*r12;
	if(p13<0) p13 = p13+m1;
	if(p12<0) p12 = p12+m1;
	x10 = x11; x11 = x12; x12 = p12-p13;
	if(x12<0) x12 = x12+m1;
	h = x20/q23;
	p23 = -a23*(x20-h*q23)-h*r23;
	h = x22/q21;
	
	p21 = a21*(x22-h*q21)-h*r21;
	if(p23<0) p23 = p23+m2;
	if(p21<0) p21 = p21+m2;
	x20 =x21; x21 =x22; x22 =p21 - p23; if(x22 <0) x22 =x22 +m2;
	if (x12<x22) return (x12-x22+m1); else return (x12-x22);
}

double get_rand() {
	int Z;
	Z=Random();
	if(Z==0)
		Z=m1;
	return (double)(4.656612873077393e-10*Z);
}

double get_gaussian()
{
	return (double)(sqrt(-2.0*log(get_rand()))*cos(6.283185307999998*get_rand()));
}

double get_st(double mu)
{
	double temp;
	temp=s0*exp(mu+sigma*sqrt(t)*get_gaussian());
	return (double)temp;
}

double get_price(double mu)
{
	double st=get_st(mu);
	if(st>k) return discont*(st-k);
	else return 0;
}

double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a=fabs(z);
	double t = 1.0/(1.0+a*p);
	double b = c2*exp((-z)*(z/2.0));
	double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
	n = 1.0-b*n;
	if ( z < 0.0 ) n = 1.0 - n;
	return n;
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
									   const double& K,       // strike (exercise) price,
									   const double& r,       // interest rate
									   const double& sigma,   // volatility
									   const double& time,
									   const double& dividend){  // time to maturity
	double time_sqrt = sqrt(time);
	double d1 = (log(S/K)+(r-dividend)*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
	double d2 = d1-(sigma*time_sqrt);
	return S*exp(-dividend*time)*N(d1) - K*exp(-r*time)*N(d2);
};

double max(double x){
	if(x>k) return x;
	else return 0;
}

int main(int argc, const char * argv[]){
	double start_s=clock();
	
	x10=3;x11=123;x12=12314;x20=76;x21=9012;x22=24889;
	long n=1000000;
	double sum1=0;
	double sum2=0;
	double mu=(r-q-0.5*sigma*sigma)*t;
	double se,selectSe=1,selectMu=0;
	
	// choose the best mu
	for(int ii=0;ii<200;ii++)
	{
		double muHead=(-0.5+ii*0.001);
		for(long i=0;i<n;i++)
		{
			double st=get_st(muHead);
			double expectation=discont*max(st)*exp(-(mu*mu-muHead*muHead+2*(mu-muHead)*log(st/s0)));
			sum1+=expectation;
			sum2+=expectation*expectation;
		}
		se=sqrt((sum2-n*(sum1/n)*(sum1/n))/(n-1))/sqrt(n);
		if(se<selectSe)
		{
			selectSe=se;
			selectMu=muHead;
		}
		std::cout<<"se="<<se<<"	mu="<<muHead<<"\n";
		sum1=sum2=0;
	}
	
	for(long i=0;i<n;i++)
	{	double muHead=-0.5;
		double st=get_st(muHead);std::cout<<st<<"  ";
		double expectation=discont*max(st)*exp(-(mu*mu-muHead*muHead+2*(mu-muHead)*log(st/s0))/(2*sigma*sigma*t));
		sum1+=expectation;
		sum2+=expectation*expectation;
	}
	se=sqrt((sum2-n*(sum1/n)*(sum1/n))/(n-1))/sqrt(n);
	
	
	//std::cout<<"choose se="<<selectSe<<"	mu="<<selectMu<<"\n"<<"\n";
	
	std::cout<<"se="<<mu<<"	mu="<<-0.0550433<<"\n";
	
	//estimates, absolute error
	//std::cout<<n<<"	"<<std::setprecision(8)<<(double)discont*sum1/n<<"	"<<std::abs(discont*sum1/n-126.93)<<"	"<<se;
	
	
	double stop_s=clock();
	double time=(stop_s-start_s)/double(CLOCKS_PER_SEC);
	//std::cout<<"	"<<time<<std::endl;
	
}
