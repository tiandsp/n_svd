#include "Tensor.h"
#include <iostream>
#include <memory>
using namespace std;

Tensor::Tensor(long a,...)
{
	long *p=&a;
	n=1;
	dim=0;
	while(*p!=0)
	{
		dim++;
		p++;
	}

	p=&a;
	N=new long[dim];
	long i=0;
	while(*p!=0)
	{
		N[i]=*p;
		n*=*p;
		i++;
		p++;
	}

	data=new double[n];
	memset(data,0,n*sizeof(double));

}

Tensor::Tensor(const Tensor &T)
{
	dim=T.dim;
	n=T.n;
	N=new long[dim];
	data=new double[n];

	for (long i=0;i<dim;i++)
	{
		N[i]=T.N[i];
	}

	for (long i=0;i<n;i++)
	{
		data[i]=T.data[i];
	}
}

Tensor::~Tensor()
{
	delete[] data;
	delete[] N;
}

Tensor &Tensor::operator=(const Tensor &T)
{

	delete[] data;
	delete[] N;

	n=T.n;
	dim=T.dim;
	N=new long[dim];
	data=new double[n];

	for (long i=0;i<dim;i++)
	{
		N[i]=T.N[i];
	}

	for (long i=0;i<n;i++)
	{
		N[i]=T.N[i];
	}
	
	return *this;
}

bool Tensor::operator==(const Tensor &T) const
{
	if (n!=T.n)
	{
		return false;
	}

	if (dim!=T.dim)
	{
		return false;
	}

	for (long i=0;i<dim;i++)
	{
		if (N[i]!=T.N[i])
		{
			return false;
		}
	}

	for (long i=0;i<n;i++)
	{
		if (data[i]!=T.data[i])
		{
			return false;
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
	if (n!=T.n || dim!=T.dim)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<dim;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}

	for (long i=0;i<n;i++)
	{
		data[i]=T.data[i]+data[i];
	}
	return *this;
}

Tensor Tensor::operator+(double a)
{
	for (long i=0;i<n;i++)
	{
		data[i]=data[i]+a;
	}
	return *this;
}

Tensor Tensor::operator-(Tensor &T)
{
	if (n!=T.n || dim!=T.dim)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<dim;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}

	for (long i=0;i<n;i++)
	{
			data[i]=data[i]-T.data[i];
	}
	return *this;
}

Tensor Tensor::operator-(double a)
{
	for (long i=0;i<n;i++)
	{
		data[i]=data[i]-a;
	}
	return *this;
}

Tensor Tensor::operator*(Tensor &T)
{
	if (n!=T.n || dim!=T.dim)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<dim;i++)
	{
		if (N[i]!=T.N[i])
		{
			cout<<"dim must equal"<<endl;
			return *this;
		}
	}

	for (long i=0;i<n;i++)
	{
		data[i]=data[i]*T.data[i];
	}
	return *this;

}

Tensor Tensor::operator*(double a)
{
	for (long i=0;i<n;i++)
	{
		data[i]=data[i]*a;
	}
	return *this;
}

Tensor Tensor::operator/(Tensor &T)
{
	if (n!=T.n || dim!=T.dim)
	{
		cout<<"dim must equal"<<endl;
		return *this;
	}

	for (long i=0;i<dim;i++)
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
		if (T.data[i]==0)
		{
			cout<<"divison zero"<<endl;
			return tmp;
		}
		else
		{
			data[i]=data[i]/T.data[i];
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
		data[i]=data[i]/a;
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
	Mat tmp(1,dim);
	for (long i=0;i<dim;i++)
	{
		tmp.setElement(N[i],0,i);
	}
	return tmp;
}

long Tensor::getDim()
{
	return dim;
}

void Tensor::setElement(double a,long b,...)
{
	long num=0;
	long *NUM;
	long *p;
	p=&b;
	NUM=new long[n];
	while (*p!=0)
	{
		num++;
		if (num>n)
		{
			cout<<"dim must be equal";
			delete[] NUM;
			return;
		}

		NUM[num]=*p;
		p++;
	}
	if (num!=n)
	{
		cout<<"dim must be equal";
		delete[] NUM;
		return;
	}

	for (long i=0;i<n;i++)
	{
		for (long j=0;j<N[i];j++)
		{


		}

	}





}

double Tensor::getElement(long b,...)
{

	return 0;
}