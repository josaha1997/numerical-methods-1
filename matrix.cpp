#include "matrix.hpp"
#include <fstream>
#include <stdlib.h>
#include <cmath>
using namespace std;

//Copy Constructor
Matrix ::Matrix(const Matrix &obj,unsigned int row,unsigned int column)
{
	if(row==0)
		row=obj.row;
	if(column==0)
		column=obj.column;
    //Dynamic memory allocation (Deep copy)
    dynamicMemory(row,column);
	for(unsigned short int i=0;i<row;i++)
		for(unsigned short int j=0;j<column;j++)
			matrix[i][j]=obj.matrix[i][j];
}

//Default Constructor
Matrix ::Matrix(unsigned int row, unsigned int column)
{
    //Dynamic memory allocation
	dynamicMemory(row,column);
}

//Destructor
Matrix ::~Matrix()
{
	//Deleting dynamic memory requested
	for(int i=0;i<row;++i)
		delete[] matrix[i];
    delete[] matrix;
}

void Matrix::dynamicMemory(unsigned int row,unsigned int column)
{
	this->row=row;
	this->column=column;
	matrix = new double *[row];
    if(matrix==NULL)
    	exit(0);
    for (unsigned int i = 0; i < row; i++)
	{
		 matrix[i] = new double[column];
		 if(matrix[i]==NULL)
		 	exit(0);
	}
     
}

//Reading from File
int Matrix::readMatrix(char *filename,int pos)
{
	double tempVar;
	ifstream inFile(filename);
	//Check if file opened
    if (!inFile)
    {
        cout << "Unable to open file!!";
        exit(0);
    }
    //Taking file pointer to the position given
	inFile.seekg(pos, ios::beg);
    inFile >> row;
    inFile >> column;
    //Dynamic memory allocation
	dynamicMemory(row,column);
    //Writing matrix into file
    for (int i = 0; i < row; i++)
        for (int j = 0; j < column; j++)
    		if (inFile >> tempVar)
          		matrix[i][j] = tempVar;
    //Returning current position 
    int position=inFile.tellg();
    inFile.close();
    return position;
}

//Writing into file
void Matrix::writeMatrix(char *filename,std::string text)
{
	
	ofstream outFile(filename, std::ios::app);
	//Check if file opened
    if (!outFile)
    {
        cout << "Unable to open file!!";
        exit(0);
    }
    //Writing matrix into file
    outFile<<text<<"\n";
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        	outFile << matrix[i][j] << "\t";	
        outFile<<"\n";
    }
    outFile.close();
}

//Addition of two matrix
Matrix Matrix::add(const Matrix &obj)
{
	Matrix ans(row,column);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
		{
			ans.matrix[i][j]=matrix[i][j]+obj.matrix[i][j];
		}
	}
	return ans;
}

//Subtraction of two matrix
Matrix Matrix::subtract(const Matrix &obj)
{
	Matrix answer(row,column);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
		{
			answer.matrix[i][j]=matrix[i][j]-obj.matrix[i][j];
		}
	}
	return answer;
}

//Multiplication of two matrix
Matrix Matrix::multiply(const Matrix &obj)
{
	Matrix answer(row,column);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<column;j++)
		{
			for(int k=0;k<column;k++)
			answer.matrix[i][j]+=matrix[i][k]*obj.matrix[k][j];
		}
	}
	return answer;
}

//Operator overloading + to add two matrices
Matrix Matrix::operator +(const Matrix &rightOperand) const
{
	if(row!=rightOperand.row || column!=rightOperand.column)
		exit(0);
	Matrix result(row,column);
	for(int i=0;i<row;i++)
		for(int j=0;j<column;j++)
			result.matrix[i][j]=matrix[i][j]+rightOperand.matrix[i][j];
	return result;
}

//Operator overloading - to subtract two matrices
Matrix Matrix::operator -(const Matrix &rightOperand) const
{
	if(row!=rightOperand.row || column!=rightOperand.column)
		exit(0);
	Matrix result(row,column);
	for(int i=0;i<row;i++)
		for(int j=0;j<column;j++)
			result.matrix[i][j]=matrix[i][j]-rightOperand.matrix[i][j];
	return result;
}

//Operator overloading * to multiply two matrices
Matrix Matrix::operator *(const Matrix &rightOperand) const
{
	if(column!=rightOperand.row)
		exit(0);
	Matrix result(row,rightOperand.column);
	for(int i=0;i<result.row;i++)
		for(int j=0;j<result.column;j++)
			for(int k=0;k<column;k++)
				result.matrix[i][j]+=matrix[i][k]*rightOperand.matrix[k][j];
	return result;
}

//Overloading * for scalar multiplication
Matrix operator *(int scalar,const Matrix &matrix)
{
	Matrix result(matrix.row,matrix.column);
	for(int i=0;i<result.row;++i)
		for(int j=0;j<result.column;++j)
			result.matrix[i][j]=scalar*matrix.matrix[i][j];
	return result;
}

//Divide matrix by scalar value
Matrix Matrix:: operator /(double scalar)const
{
	Matrix result(row,column);
	for(int i=0;i<row;++i)
		for(int j=0;j<column;++j)
			result.matrix[i][j]=matrix[i][j]/scalar;
	return result;
}

//Operator overloading = to reassign object
void Matrix::operator =(const Matrix &rightOperand)
{
	if(this!=&rightOperand)
	{
		for(int i=0;i<row;++i)
			delete[] matrix[i];
    	delete[] matrix;
    	matrix=0;
		dynamicMemory(rightOperand.row,rightOperand.column);
		for(int i=0;i<row;++i)
			for(int j=0;j<column;++j)
				matrix[i][j]=rightOperand.matrix[i][j];
	}	
}

//Operator overlading == to check equality of two matrices
bool Matrix::operator ==(const Matrix &rightOperand) const
{
	if(row==rightOperand.row && column==rightOperand.column)
	{
		for(unsigned short int i=0;i<row;++i)
			for(unsigned short int j=0;j<column;++j)
				if(matrix[i][j]!=rightOperand.matrix[i][j])
					return false;
		return true;
	}
	else
		return false;
}
Array Matrix::operator[](int rno)
{
	return Array(matrix[rno]);
}

//Scalar Multiplication of two matrix
Matrix Matrix::scalarMultiply(double scalar_value)
{
	Matrix answer=*this;
	for(int i=0;i<row;i++)
		for(int j=0;j<column;j++)
			answer.matrix[i][j]*=scalar_value;
	return answer;
}

//Transpose of matrix
Matrix Matrix::transpose()
{
	Matrix result(column,row);
	for(int i=0;i<result.row;++i)
		for(int j=0;j<result.column;++j)
			result.matrix[i][j]=matrix[j][i];
	return result;
}

//Gauss
Matrix Matrix::upperTriangularMatrix()const
{
	Matrix upperT=*this;
	for(int i=0;i<upperT.row;++i)
	{
		double factorD=upperT.matrix[i][i];
		for(int j=i;j<upperT.row;++j)
		{
			double factorM=upperT.matrix[j][i];
			for(int k=0;k<upperT.column;++k)
			{
				if(j==i)	
					upperT.matrix[i][k]/=factorD;
				if(j>i)
					upperT.matrix[j][k]-=upperT.matrix[i][k]*factorM;
			}
		}
	} 
	return upperT;
}
Matrix Matrix::rowMerge(const Matrix &rightMatrix)const
{
	Matrix result(row,column+rightMatrix.column);
	for(int i=0;i<row;++i)
	{
		for(int j=0;j<column;++j)
			result[i][j]=matrix[i][j];
		int j=0;
		for(int k=column;k<result.column;++k)
			result[i][k]=rightMatrix.matrix[i][j++];
	}	
	return result;
}
Matrix Matrix::inverse()const
{

	Matrix inverse=upperTriangularMatrix();

	Matrix inverse1;
	inverse1=inverse.getSubMatrix(0,row,row);

	inverse=inverse.getSubMatrix(0,0,row);

	inverse1=inverse1.transpose();
	inverse=inverse.transpose();
	inverse=inverse.rowMerge(inverse1);

	inverse=inverse.upperTriangularMatrix();

	inverse=inverse.getSubMatrix(0,row,row);
		inverse.display();
	return inverse;
	
}
double Matrix::trace()
{
	double trace=0.0;
	for(unsigned short int i=0;i<row;++i)
		trace+=matrix[i][i];
	return trace;
}

bool Matrix::isOrthogonal()
{
	Matrix transpose=this->transpose();
	Matrix result=transpose*(*this);
	return(result.isIdentity());
}
Matrix Matrix::gaussianElemination()const
{
	double temp;
	Matrix solution(row,1);
	for(int i=0;i<row;++i)
		solution.matrix[i][0]=0.0;
	Matrix upperT=upperTriangularMatrix();
	solution.matrix[row-1][0]=upperT.matrix[row-1][column-1];
	cout<<solution.matrix[row-1][0];
	for(int i=row-2;i>=0;--i)
	{
		temp=0.0;
		for(int j=column-2;j>=0;--j)
			temp-=upperT.matrix[i][j]*solution.matrix[j][0];
		solution.matrix[i][0]=upperT.matrix[i][column-1]+temp;
	}
	return solution;	
}

//Gauss Jacobi
Matrix Matrix:: gaussJacobi()const
{

	Matrix sqmatrix(*this,row,column-1);
	Matrix position=sqmatrix.findDiagonalyDominant();
	if(!position.isDiagonalyDominant())
	{
		cout<<"Not a diagonally dominant matrix";
		exit(0);
	}
	int iteration=1;
	//cout<<iteration;
	Matrix solution(15,position.column);
	for(int i=0;i<position.column;++i)
	{	
		//cout<<"position"<<position.matrix[0][i];
		solution.matrix[0][i]=0.0;
		//cout<<"\n\n->"<<solution.matrix[0][i];
	}
	while(iteration<15)
	{
		for(int i=0;i<row;++i)
		{
			double factor=matrix[i][int(position.matrix[0][i])];
			//cout<<factor<<"\n";
			solution.matrix[iteration][int(position.matrix[0][i])]=matrix[i][column-1]/factor;
			//cout<<solution.matrix[iteration][int(position.matrix[0][i])]<<"\n";
			for(int j=0;j<column-1;++j)
			{
				if(j!=position.matrix[0][i])
					solution.matrix[iteration][int(position.matrix[0][i])]-=(matrix[i][j]/factor)*solution.matrix[iteration-1][j];
			}
			cout<<solution.matrix[iteration][int(position.matrix[0][i])]<<"\t";
		}
		++iteration;
		cout<<"\n";
	}
	return solution;
}

//Gauss Seidal
Matrix Matrix:: gaussSeidal()const
{

	Matrix sqmatrix(*this,row,column-1);
	Matrix position=sqmatrix.findDiagonalyDominant();
	if(!position.isDiagonalyDominant())
	{
		cout<<"Not a diagonally dominant matrix";
		exit(0);
	}
	int iteration=1;
	int k=0;
	//cout<<iteration;
	Matrix solution(15,position.column);
	for(int i=0;i<position.column;++i)
	{	
		//cout<<"position"<<position.matrix[0][i];
		solution.matrix[0][i]=0.0;
		//cout<<"\n\n->"<<solution.matrix[0][i];
	}
	while(iteration<15)
	{
		for(int i=0;i<row;++i)
		{
			double factor=matrix[i][int(position.matrix[0][i])];
			//cout<<factor<<"\n";
			solution.matrix[iteration][int(position.matrix[0][i])]=matrix[i][column-1]/factor;
			//cout<<solution.matrix[iteration][int(position.matrix[0][i])]<<"\n";
			k=position.matrix[0][i];
			for(int j=0;j<column-1;++j)
			{
				if(j!=k)
				{
					if(k==0 || j>k)
						solution.matrix[iteration][int(position.matrix[0][i])]-=(matrix[i][j]/factor)*solution.matrix[iteration-1][j];
					if(k!=0 && j<k)
						solution.matrix[iteration][int(position.matrix[0][i])]-=(matrix[i][j]/factor)*solution.matrix[iteration][j];
					
							
				}
			}
			cout<<solution.matrix[iteration][int(position.matrix[0][i])]<<"\t";
		}
		++iteration;
		cout<<"\n";
	}
	return solution;
}

double Matrix::powerMethod(const Matrix &x0)const
{
	Matrix sqmatrix(*this,row,column-1);
	int k=0;
	double eigenmax;
	x0.display();
	Matrix Ax;
	Ax=sqmatrix*x0;
	Ax.display();
	Matrix x1;
	while(k<=10)
	{
		eigenmax=Ax[0][0];
		for(int i=1;i<Ax.row;++i)
			if(abs(eigenmax)<abs(Ax[i][0]))
				eigenmax=Ax[i][0];
		x1=(Ax)/eigenmax;
		x1.display();
		Ax=sqmatrix*x1;
		k++;
	}
	for(int i=1;i<Ax.row;++i)
			if(abs(eigenmax)<abs(Ax[i][0]))
				eigenmax=Ax[i][0];
	return eigenmax;
	
}
//Check for Identity Matrix
bool Matrix::isIdentity()
{
	//Check if square matrix
	if(!isSquare())
		return false;
	for(int i=0;i<row;++i)
	{
		for(int j=0;j<column;++j)
		{
			//Check for diagonal elements are 1
			if(i==j && matrix[i][j]!=1)
				return false;
			//Check if non-diagonal elements are 0
			else
				if(i!=j && matrix[i][j]!=0)
					return false;
		}
	}
	return true;
} 

//Check for Square Matrix
bool Matrix::isSquare() const
{
	//Check nom of rows and no of columns are equal
	if(row==column)
		return true;
	return false;
}

//Check for Symmentric Matrix
bool Matrix::isSymmetric()
{
	//Check if square matrix
	if(!isSquare())
		return false;
		for(int i=0;i<row;++i)
			for(int j=0;j<column;++j)
					if(i!=j && matrix[i][j]!=matrix[j][i])
						return false;
	return true;
}

//Check for Null Matrix
bool Matrix::isNull()
{
	for(int i=0;i<row;++i)
		for(int j=0;j<column;++j)
			//Check all elements are 0
			if( matrix[i][j]!=0)
				return false;
	return true;
}

//Check for Diagonal Matrix
bool Matrix::isDiagonal()
{
	//Check if square matrix
	if(!isSquare())
		return false;
	for(int i=0;i<row;++i)
		for(int j=0;j<column;++j)
			//Check if non-diagonal elements are 0
			if(i!=j && matrix[i][j]!=0)
				return false;
	return true;
}

//Returning dominant value's position of each row
Matrix Matrix::findDiagonalyDominant() const

{
	//Check if square matrix
	if(!isSquare())
		return false;
	double sum;
	Matrix position(1,row);
	/*int *position=new int[row];
	if(position==NULL)
		exit(0);*/
	for(int i=0;i<row;++i)
	{
		double max=abs(matrix[i][0]);
		position.matrix[0][i]=0;
		for(int j=1;j<column;++j)
		{
			if(abs(matrix[i][j])>max)
			{
				max=abs(matrix[i][j]);
				position.matrix[0][i]=j;
			}
		}
	}
	for(int i=0;i<row;++i)
	{
		sum=0.0;
		for(int j=0;j<column;++j)
		{
			if(j!=position.matrix[0][i])
				sum+=abs(matrix[i][j]);
		}
		if(abs(matrix[i][int(position.matrix[0][i])])<sum)
			position.matrix[0][i]=-1;
	}
	return position;
}

//Check if diagonal element is greater than sum of other elements in the row
bool Matrix::isDiagonalyDominant() const
{	
	for(int i=0;i<row;++i)
	{
		cout<<"position"<<i<<matrix[0][i]<<"\n";
		int pos=matrix[0][i];	
		for(int j=i+1;j<row;++j)
			if(matrix[0][j]==-1 || pos==matrix[0][j]||pos==-1)
				return false;
	}
	return true;
}

//Display elements of matrix
void Matrix::display()const
{
	cout<<"Matrix\n Row="<<row<<"\tColumn="<<column<<"\n";
	for(int i=0;i<row;i++)
	{	for(int j=0;j<column;j++)
			cout<<"\t"<<matrix[i][j];
		cout<<"\n";
	}
}
Matrix Matrix ::getSubMatrix(unsigned rowIndex, unsigned columnIndex, unsigned size)
{
    Matrix temp(size, size);
    int icount, jcount, tempColumn = columnIndex;
    for (icount = 0; icount < size; icount++)
    {
        columnIndex = tempColumn;
        for (jcount = 0; jcount < size; jcount++)
        {
            temp.matrix[icount][jcount] = this->matrix[rowIndex][columnIndex++];
        }
        rowIndex++;
    }
    return temp;
}
double Matrix ::determinant()
{
    double det;
    if (this->row == 2 && this->column == 2)
        det = ((this->matrix[0][0] * this->matrix[1][1]) - (this->matrix[0][1] * this->matrix[1][0]));
    else
        det = -0;
    return det;
}
Matrix Matrix ::operator/(const Matrix &obj)
{
    Matrix temp(obj.row, obj.column);
    for (int i = 0; i < obj.row; i++)
    {
        for (int j = 0; j < obj.column; j++)
        {
            temp.matrix[i][j] = this->matrix[i][j] / obj.matrix[i][j];
        }
    }
    return temp;
}
/*ifstream&  operator>>(ifstream& in,const Matrix &obj)
{
    unsigned long tempVar;
    in>>obj.row;
    in>>obj.column;
    obj(obj.row,obj.column);
   // obj.matrix=new double*[obj.row];
    //for(unsigned int i=0;i<obj.row;i++)
     //   obj.matrix[i]=new double[obj.column];
   
    for(int i=0;i<obj.row;i++)
    {
        for(int j=0;j<obj.column;j++)
        {
            if(in>>tempVar)
                obj.matrix[i][j]=tempVar;
            else
                break;
        }
    }
    return in;
}*/
