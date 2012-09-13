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


bool Tensor::operator!=(const Tensor &T) const
{
	if (*this==T)
	{
		return false;
	}
	return true;
}

Tensor Tensor::operator+(Tensor &T)
{
	if (n!=T.n)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<n;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=T.data[i][j]+data[i][j];
		}
	}
	return *this;
}

Tensor Tensor::operator+(double a)
{
	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=data[i][j]+a;
		}
	}
	return *this;
}

Tensor Tensor::operator-(Tensor &T)
{
	if (n!=T.n)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<n;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=data[i][j]-T.data[i][j];
		}
	}
	return *this;
}

Tensor Tensor::operator-(double a)
{
	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=data[i][j]-a;
		}
	}
	return *this;
}

Tensor Tensor::operator*(Tensor &T)
{
	if (n!=T.n)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<n;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=data[i][j]*T.data[i][j];
		}
	}
	return *this;

}

Tensor Tensor::operator*(double a)
{
	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=data[i][j]*a;
		}
	}
	return *this;
}

Tensor Tensor::operator/(Tensor &T)
{
	if (n!=T.n)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<n;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}
	Tensor tmp(*this);

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			if (T.data[i][j]==0)
			{
				cout<<"divison zero"<<endl;
				return tmp;
			}
			else
			{
				data[i][j]=data[i][j]/T.data[i][j];
			}
		}
	}
	return *this;
}

Tensor Tensor::operator/(double a)
{
	if (a==0)
	{
		cout<<"division zero"<<endl;
		return *this;
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{
			data[i][j]=data[i][j]/a;
		}
	}
	return *this;
}

N_SVD *Tensor::m_mode_svd()
{
	N_SVD *nsvd;

	return nsvd;
}

Mat Tensor::getSize()
{
	Mat tmp(1,n);
	for (long i=0;i<n;i++)
	{
		tmp.setElement(N[i],0,i);
	}
	return tmp;
}

long Tensor::getDim()
{
	return n;
}

void Tensor::setElement(double a,long b,...)
{
	




}

double Tensor::getElement(long b,...)
{

	return 0;
}