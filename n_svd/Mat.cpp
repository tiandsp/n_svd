#include "Mat.h"
#include <memory>
#include <iostream>
using namespace std;
Mat::Mat(long r,long c)
{
	row=r;
	col=c;
	size=r*c;
	data=new double *[row];
	for(long i=0;i<row;i++)
	{
		data[i]=new double[col];
		memset(data[i],0,c*sizeof(double));
	}
}

Mat::Mat(const Mat &M)
{
	row=M.row;
	col=M.col;
	size=M.size;
	data=new double *[row];
	for(long i=0;i<row;i++)
	{
		data[i]=new double[col];
	}

	for (long i=0;i<row;i++)
	{
		for(long j=0;j<col;j++)
			data[i][j]=M.data[i][j];
	}
}


Mat::~Mat()
{
	for (long i=0;i<row;i++)
	{
		delete[] data[i];
	}
	delete[] data;
}

long Mat::getRow() const
{
	return row;
}

long Mat::getCol() const
{
	return col;
}

long Mat::getSize() const
{
	return size;
}

Mat &Mat::operator=(const Mat &M)
{
	for (long i=0;i<row;i++)
	{
		delete[] data[i];
	}
	delete[] data;

	row=M.row;
	col=M.col;
	size=M.size;
	data=new double *[row];
	for(long i=0;i<row;i++)
	{
		data[i]=new double[col];
	}

	for (long i=0;i<row;i++)
	{
		for(long j=0;j<col;j++)
			data[i][j]=M.data[i][j];
	}
	return *this;
}

bool Mat::operator==(const Mat &M) const
{
	if (col!=M.col)
		return false;
	if (row!=M.row)
		return false;
	if (size!=M.size)
		return false;
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			if (data[i][j]!=M.data[i][j])
			{
				return false;
			}
		}
	}
	return true;
}

bool Mat::operator!=(const Mat &M) const
{
	if (*this==M)
	{
		return false;
	}
	return true;
}

Mat Mat::operator~()
{
	Mat tmp(col,row);

	for (long j=0;j<col;j++)
	{
		for (long i=0;i<row;i++)
		{
			tmp.data[j][i]=data[i][j];
		}
	}
	return tmp;
}

Mat Mat::operator+(Mat &M)
{
	if (row!=M.row || col!=M.col)
	{
		cout<<"row or col not match"<<endl;
		return *this;
	}
	
	Mat tmp(row,col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j] +=M.data[i][j];
		}
	}
	return tmp;
}

Mat Mat::operator-(Mat &M)
{
	if (row!=M.row || col!=M.col)
	{
		cout<<"row or col not match"<<endl;
		return *this;
	}

	Mat tmp(row,col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j] -=M.data[i][j];
		}
	}
	return tmp;
}

Mat Mat::operator*(Mat &M)
{
	if (col!=M.row)
	{
		cout<<"row or col not match"<<endl;
		return *this;
	}

	Mat tmp(row,M.col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			for (long k=0;k<M.col;k++)
			{
				tmp.data[i][k] += (data[i][j] *M.data[j][k]);
			}
		}
	}
	return tmp;
}

Mat Mat::operator/(Mat &M)
{

	return *this;
}

Mat Mat::dotMultiplication(Mat &M)
{
	if (row!=M.row || col!=M.col)
	{
		cout<<"row or col not match"<<endl;
		return *this;
	}

	Mat tmp(row,col);
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j]=data[i][j]*M.data[i][j];
		}
	}
	return tmp;
}

Mat Mat::dotDivision(Mat &M)
{
	if (row!=M.row || col!=M.col)
	{
		cout<<"row or col not match"<<endl;
		return *this;
	}
	
	Mat tmp(row,col);
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j]=data[i][j]/M.data[i][j];
		}
	}
	return tmp;	
}

Mat Mat::inv()
{


	return *this;
}

void Mat::setElement(double n,long r,long c)
{
	if (r<0 || c<0 || r>row || c>col)
	{
		cout<<"overflow"<<endl;
		return;
	}
	data[r][c]=n;
}

double Mat::getElement(long r,long c)
{
	if (r<0 || c<0 || r>row || c>col)
	{
		cout<<"overflow"<<endl;
		return 0;
	}

	return data[r][c];
}

void Mat::print()
{
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
			cout<<data[i][j]<<"  ";
		cout<<endl;
	}
}