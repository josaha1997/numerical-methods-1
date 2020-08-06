#include<iostream>
#include<cmath>
#include<stdlib.h>
#include<iomanip>
#include"poly.hpp"

using namespace std;

polynomial::polynomial(int degree)
{
	if(degree<=0)
	{
		cout<<"Invalid degree";
		exit(0);
	}
	this->degree=degree;
	poly=new float[degree+1];
}
void polynomial::acceptPoly()
{
	cout<<"Enter the coefficients of the polynomial:\n";
	for(int i=degree;i>=0;--i)
	{
		cout<<"x^"<<i<<":";
		cin>>poly[i];
	}
}
float polynomial::evaluatePoly(float x)const
{
	float answer=0.0;
	for(int i=0;i<=degree;++i)
		answer+=pow(x,i)*poly[i];
	return answer;
}

void polynomial::displayPoly()const
{
	for(int i=0;i<degree;++i)
		cout<<"("<<poly[i]<<")x^"<<i<<" + ";
	cout<<"("<<poly[degree]<<")x^"<<degree<<"\n";
}
float polynomial::findRoot(char token)const
{
	cout<<"\nEnter the interval for root of the function[a,b]:";
	float a,b;
	cin>>a>>b;
	displayPoly();
	cout<<"\na="<<a<<"\tb="<<b;
	float fa=evaluatePoly(a),fb=evaluatePoly(b),tol=0.00001;
	if(abs(fa)<=tol)
		return a;
	if(abs(fb)<=tol)
		return b;
		
	if(fa*fb>0)
	{
		cout<<"Given Interval is wrong!\n";
		exit(0);
	}

	cout<<"\n"<<setw(15)<<"a"<<setw(15)<<"b"<<setw(15)<<"c"<<setw(15)<<"f(a)"<<setw(15)<<"f(b)"<<setw(15)<<"f(c)"<<setw(15)<<"Update\n";
	float c=0.0;
	if(token=='b')
		c=(a+b)/2;
	if(token=='r')
		c=(a*fb-b*fa)/(fb-fa);
	while(abs(evaluatePoly(c))>=tol)
	{
		cout<<setprecision(8);
		cout<<fixed;
		cout<<setw(15)<<a<<setw(15)<<b<<setw(15)<<c<<setw(15)<<fa<<setw(15)<<fb<<setw(15)<<evaluatePoly(c)<<setw(15);
		
		if(fa*evaluatePoly(c)<0)
		{
			b=c;
			fb=evaluatePoly(b);
			cout<<"b=c\n";
		}
		if(fb*evaluatePoly(c)<0)
		{
			a=c;
			fa=evaluatePoly(a);
			cout<<"a=c\n";
		}
		if(token=='b')
			c=(a+b)/2;
		if(token=='r')
			c=(a*fb-b*fa)/(fb-fa);
		
	}
	return c;
}
polynomial polynomial::derivative()const 
{
	polynomial derivative(degree-1);
	for(int i=0;i<=derivative.degree;++i)
		derivative.poly[i]=poly[i+1]*(i+1);
	return derivative;
}
void polynomial::operator=(const polynomial &obj)
{
	degree=obj.degree;
	delete[]poly;
	poly=new float[degree+1];
	for(int i=0;i<=degree;++i)
		poly[i]=obj.poly[i];
}
float polynomial::newtonRaphson(float initialpoint)const
{
	float x=initialpoint,x1=0.0;
	polynomial derivative;
	derivative=this->derivative();
	cout<<"Points\n"<<x<<"\n";
	x1=x-(evaluatePoly(x)/derivative.evaluatePoly(x));
	while(abs(x1-x)>=0.000001)
	{
		cout<<x1<<"\n";
		x=x1;
		x1=x-(evaluatePoly(x)/derivative.evaluatePoly(x));
	}
	return x1;	
}
float polynomial::secant(float initialpoint1,float initialpoint2)const
{
	float x=initialpoint1,x1=initialpoint2,x2;
	cout<<"Points\n"<<x<<"\n"<<x1<<"\n";
	x2=(x*evaluatePoly(x1)-x1*evaluatePoly(x))/(evaluatePoly(x1)-evaluatePoly(x));
	while(abs(x2-x1)>=0.000001)
	{
		cout<<x2<<"\n";
		x=x1;
		x1=x2;
		x2=(x*evaluatePoly(x1)-x1*evaluatePoly(x))/(evaluatePoly(x1)-evaluatePoly(x));
	}
	return x2;	
}


