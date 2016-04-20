//
//  rv.cpp
//  IE525Project
//
//  Created by fushuyue on 4/19/16.
//  Copyright © 2016 Monte Carlo Simulation. All rights reserved.
//
#include "interface.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
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

double min(double a,double b)
{
	if (a<=b)
		return a;
	else
		return b;
};

double max(double a,double b)
{
	if (a<=b)
		return b;
	else
		return a;
};

int Random()
{
	int h, p12, p13, p21, p23;
	/* Component 1 */
	h = x10/q13; p13 = -a13*(x10-h*q13)-h*r13;
	h = x11/q12; p12 = a12*(x11-h*q12)-h*r12;
	if(p13<0) p13 = p13+m1; if(p12<0) p12 = p12+m1;
	x10 = x11; x11 = x12; x12 = p12-p13; if(x12<0) x12 = x12+m1;
	/* Component 2 */
	h = x20/q23; p23 = -a23*(x20-h*q23)-h*r23;
	h = x22/q21; p21 = a21*(x22-h*q21)-h*r21;
	if(p23<0) p23 = p23+m2; if(p21<0) p21 = p21+m2;
	/* Combination */
	if (x12<x22) return (x12-x22+m1); else return (x12-x22);
}
double Uniform()
{
	return (double)(rand()/(((double)pow(2,31)-1)));
}

double gaussian_box_muller()
{
	return (double)(sqrt(-2.0*log(Uniform()))*cos(6.283185307999998*Uniform()));
}


double Normalcdf(const double& z) {
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

double inverse(double u){
	
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



