#include<iostream>
#include"poly.hpp"
#include<cmath>
using namespace std;
float evaluateFunc(float x)
{
	float func=3*x*pow(2,-x);
	return func;
}

Matrix acceptCoordinates(int degree)
{
	Matrix points(degree+1,2);
	float a,b;
	cout<<"Enter the interval to interpolate the function [a,b]:";
	cin>>a>>b;
	float difference= (b-a)/degree;
	points.matrix[0][0]=a;
	points.matrix[0][1]=evaluateFunc(a);
	for(int i=1;i<points.row;++i)
	{
		a+=difference;
		points.matrix[i][0]=a;
		points.matrix[i][1]=evaluateFunc(a);	
	}
	points.display();
	return points;	
}
float newtondifference(const Matrix &dd,const Matrix &points,float x)
{
	float result;
	result=dd.matrix[0][0];
	cout<<result<<"+";
	float prod=1.0;
	for(int i=1;i<=points.row-1;++i)
	{
		prod=1.0;
		for(int j=0;j<i;++j)
		{
			prod*=(x-points.matrix[j][0]);
			cout<<"(x-"<<points.matrix[j][0]<<")";
		}
		cout<<"("<<dd.matrix[i][i]<<")"<<"+";		
		result+=dd.matrix[i][i]*prod;
	}
	cout<<"\b";
	return result;
}
Matrix differencetable(const Matrix &points)
{
	Matrix dd(points.row,points.row);
	for(int i=0;i<dd.column;++i)
	{
		for(int j=0;j<dd.row;++j)
		{
			if(i==0)
				dd.matrix[j][0]=points.matrix[j][1];	
			else if(j>=i && i!=0)
				dd.matrix[j][i]=(dd.matrix[j][i-1]-dd.matrix[j-1][i-1])/(points.matrix[j][0]-points.matrix[j-i][0]);
			else
				dd.matrix[j][i]=0;
		}
	}
	dd.display();
	return dd;
}
int main()
{
	int degree;
	//cout<<"\nEnter the degree of the polynomial:";
	//cin>>degree;
	polynomial P(degree),derivative;
	//P.acceptPoly();
	//P.displayPoly();
	//cout<<"Root(bisection method)="<<P.findRoot('b')<<"\n";
	//cout<<"Root(regula falsie)="<<P.findRoot('r')<<"\n";
	//cout<<"Root(newton raphson)="<<P.newtonRaphson(1)<<"\n";
	//cout<<"Root(secant)="<<P.secant(1,2)<<"\n";
	cout<<"\nEnter the degree of the interpolated polynomial:";
	cin>>degree;
//	polynomial interpolate(degree);
	Matrix points=acceptCoordinates(degree);
	Matrix dd=differencetable(points);
	cout<<"Interpolated F(5):"<<newtondifference(dd,points,5);
	cout<<"\n"<<evaluateFunc(5);
	return 0;
}
