#include "Tensor.h"
#include <iostream>
#include <memory>
using namespace std;

Tensor::Tensor(long a,...)
{
	long *p=&a;
	n=0;
	while(*p!=0)
	{
		n++;
		p++;
	}

	p=&a;
	N=new long[n];
	data=new double *[n];
	long i=0;
	while(*p!=0)
	{
		N[i]=*p;
		data[i]=new double[*p];
		memset(data[i],0,(*p)*sizeof(double));
		i++;
		p++;
	}

}

Tensor::Tensor(const Tensor &T)
{
	n=T.n;
	N=new long[n];
	data=new double *[n];
	for (long i=0;i<n;i++)
	{
		N[i]=T.N[i];
		data[i]=new double[N[i]];
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=T.data[i][j];
		}
	}
}

Tensor::~Tensor()
{
	for (long i=0;i<n;i++)
	{
		delete[] data[i];
	}

	delete[] data;
	delete[] N;

}

Tensor &Tensor::operator=(const Tensor &T)
{
	for (long i=0;i<n;i++)
	{
		delete[] data[i];
	}
	delete[] data;
	delete[] N;

	n=T.n;
	N=new long[n];
	data=new double *[n];

	for (long i=0;i<n;i++)
	{
		N[i]=T.N[i];
		data[i]=new double[N[i]];
	}
	
	return *this;
}

bool Tensor::operator==(const Tensor &T) const
{
	if (n!=T.n)
	{
		return false;
	}

	for (long i=0;i<n;i++)
	{
		if (N[i]!=T.N[i])
		{
			return false;
		}
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			if (data[i][j]!=T.data[i][j])
			{
				return false;
			}
		}
	}

	return true;
}


