#include"matrix.hpp"
class polynomial
{
	public:	int degree;
			float *poly;
			polynomial(int degree=1);
			void acceptPoly();
			float evaluatePoly(float x)const;
			void displayPoly()const;
			float findRoot(char token)const;
			float newtonRaphson(float initialpoint)const;
			polynomial derivative()const;
			void operator=(const polynomial &obj);
			float secant(float initialpoint1,float initialpoint2)const;
		
};
