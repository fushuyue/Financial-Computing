#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
// constant number for simulate uniform random variables
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
int x10, x11, x12, x20, x21, x22;

// constant number from option contract
const double k=1900;
const double t=(double)1/12;
const double s0=2000;
const double q=0.02;
const double sigma=0.3;
const double r=0.005;
const double discont=exp(-r*t);

// simulate random number
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

// // simulate uniform random variables
double get_rand() {
	int Z;
	Z=Random();
	if(Z==0)
		Z=m1;
	return (double)(4.656612873077393e-10*Z);
}

// box-muller method
double get_gaussian()
{
	return (double)(sqrt(-2.0*log(get_rand()))*cos(6.283185307999998*get_rand()));
}

// geometric brownian motion model for stock prices
double get_st()
{
	double temp;
	temp=s0*exp(t*(r-q-0.5*sigma*sigma)+sigma*sqrt(t)*get_gaussian());
	return (double)temp;
}

// get option prices from stock price
// this is an asset or nothing call 
double get_price()
{
	double st=get_st();
	if(st>k) return discont*(st);
	else return 0;
}

// cdf of gaussian distribution
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

// bs formular
double option_price_call_black_scholes(const double& S,       // spot (underlying) price
									   const double& K,       // strike (exercise) price,
									   const double& r,       // interest rate
									   const double& sigma,   // volatility
									   const double& time,
									   const double& dividend){  // time to maturity
	double time_sqrt = sqrt(time);
	double d1 = (log(S/K)+(r-dividend)*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
	double d2 = d1-(sigma*time_sqrt);
	//return S*exp(-dividend*time)*N(d1) - K*exp(-r*time)*N(d2);
	return S*exp(-dividend*time)*N(d1);
};

double max(double x){
	if(x>k) return x;
	else return 0;
}


int main(int argc, const char * argv[]){
	double start_s=clock();
	
	// calculate the optimal weight for control variates
	double b=0;
	double xBar=0,yBar=0;
	std::vector <double> x;
	std::vector <double> y;
	x10=3;x11=123;x12=12314;x20=76;x21=9012;x22=24889;
	
	// 
	for(int j=1;j<=1000;j++)
	{
		double temp=get_st();
		x.push_back(temp);
		xBar=(xBar*(j-1)+temp)/j;
		double temp1=max(temp);
		y.push_back(temp1);
		yBar=(yBar*(j-1)+temp1)/j;
	}
	
	double temp1=0,temp2=0;
	
	for(int j=1;j<=1000;j++)
	{
		double xTemp=x[j];
		double yTemp=y[j];
		
		temp1=temp1+(xTemp-xBar)*(yTemp-yBar);
		temp2=(xTemp-xBar)*(xTemp-xBar)+temp2;
	}
	
	b=(double)temp1/temp2;
	

	long n=atoi(argv[1]);
	double sum1=0;
	double sum2=0;
	
	//std::cout<<"for the first call price= "<<std::setprecision(8)<<option_price_call_black_scholes(s0, k, r, sigma, t,q)<<"\n";
	//std::cout<<"for the second call price= "<<std::setprecision(8)<<option_price_call_black_scholes(s0, 2200, r, sigma, t,q)<<"\n";
	
	// simulate stock prices and calculate option price
	// Y = Ci+b(E[X]-Xi)
	// E[Y]=E[Ci]

	for(long i=0;i<n;i++)
	{
		double st=get_st();
		//double expectation=(max(st)+b*(s0*exp(t*(r-q))-st))*discont;
		double expectation=max(st)*discont;
		sum1+=expectation;
		sum2+=expectation*expectation;
	}
	
	// standard deviation
	double se;
	se=sqrt((sum2-n*(sum1/n)*(sum1/n))/(n-1))/sqrt(n);
	
	//estimates, absolute error
	std::cout<<n<<"	"<<std::setprecision(8)<<(double)sum1/n<<"	"<<std::abs(sum1/n-1463.0597)<<"	"<<se;
	
	
	double stop_s=clock();
	double time=(stop_s-start_s)/double(CLOCKS_PER_SEC);
	std::cout<<"	"<<time<<std::endl;

}
