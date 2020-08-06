#include <fstream>
#include <string>
#include <iostream>

class Array 
{
	private:double* array;
	public:
		Array(double* array) : array(array) { }
		double& operator[](unsigned int index) 
		{
		    return array[index];
		}

};
class Matrix
{
	public:
			unsigned int row, column;
			double **matrix;
			//Default Constructor
			Matrix(unsigned int row=1, unsigned int column=1);
			//Copy Constructor
			Matrix(const Matrix &obj,unsigned int row=0,unsigned int column=0); 
			void dynamicMemory(unsigned int row,unsigned int column);
			//Read from File 
			int readMatrix(char *,int);
			//Write into file
			void writeMatrix(char *,std::string);
			//Arithmatic operations on matrix
			Matrix add(const Matrix &obj);
			Matrix subtract(const Matrix &obj);
			Matrix multiply(const Matrix &obj); 

			//Operator overloading
			Matrix operator +(const Matrix &obj) const;
			Matrix operator -(const Matrix &obj) const;
			Matrix operator *(const Matrix &obj) const;
			void operator =(const Matrix &obj);
			bool operator ==(const Matrix &obj) const;
			Array operator[](int rno);
			friend Matrix operator *(int scalar,const Matrix &obj);
			Matrix operator /(double scalar) const;

			//Gaussian Elemination
			Matrix gaussianElemination()const;
			Matrix upperTriangularMatrix()const;
			Matrix gaussJacobi()const;
			Matrix gaussSeidal()const;
			double powerMethod(const Matrix &obj)const;
			
			//Properties of matrix
			Matrix inverse()const;
			Matrix rowMerge(const Matrix &obj)const;
			Matrix scalarMultiply(double);
			Matrix transpose();
			Matrix findDiagonalyDominant()const;
			bool isIdentity();
			double trace();
			bool isOrthogonal();
			bool isSquare()const;
			bool isSymmetric();
			bool isNull(); 
			bool isDiagonal(); 
			bool isDiagonalyDominant()const;
			void display()const;
			//Destructor
			~Matrix();
			Matrix getSubMatrix(unsigned introwIndex, unsigned int columnIndex, unsigned int size);
			double determinant();
			//friend ifstream& operator >> (ifstream& in,const Matrix &obj);
			Matrix operator/(const Matrix &obj);
};

