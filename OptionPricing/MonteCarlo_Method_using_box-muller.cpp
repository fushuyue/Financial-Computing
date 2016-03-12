// box-muller
#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
using namespace std;
const double k=100;
const double t=(double)1/52;
const double s0=100;
const double q=0;
const double sigma=0.3;
const double r=0.005;
const double discont=exp(-r*t);
double start_s=clock();
double get_rand()
{
	return (double)(random()/(((double)pow(2,31)-1)));
}
//box muller
double get_gaussian()
{
	return (double)(sqrt(-2.0*log(get_rand()))*cos(6.283185307999998*get_rand()));
}
double get_st()
{
	double temp;
	temp=exp(t*(r-q-0.5*sigma*sigma)+sigma*sqrt(t)*get_gaussian());
	return (double)s0*temp;
}

double get_price()
{
	double st=get_st();
	if(st<=k) return discont*(-st+k);
	else return 0;
}

int main(int argc, const char * argv[]) {
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
	cout<<n<<"	"<<sum/n<<"	"<<std::abs(sum/n-1.6547)<<"	"<<se;
	double stop_s=clock();
	double time=(stop_s-start_s)/double(CLOCKS_PER_SEC);
	cout<<"	"<<time<<endl;
	
}