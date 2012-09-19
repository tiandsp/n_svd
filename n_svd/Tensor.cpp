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
	yuji=new long[dim];
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
	int tmp;
	for (int i=0;i<dim;i++)
	{
		tmp=1;
		for (int j=i+1;j<dim;j++)
		{
			tmp*=N[j];
		}
		yuji[i]=tmp;
	}

}

Tensor::Tensor(const Tensor &T)
{
	dim=T.dim;
	n=T.n;
	N=new long[dim];
	yuji=new long[dim];
	data=new double[n];
	
	for (long i=0;i<dim;i++)
	{
		N[i]=T.N[i];
		yuji[i]=T.yuji[i];
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
	delete[] yuji;
}

Tensor &Tensor::operator=(const Tensor &T)
{

	delete[] data;
	delete[] N;
	delete[] yuji;

	n=T.n;
	dim=T.dim;
	N=new long[dim];
	yuji=new long[dim];
	data=new double[n];

	for (long i=0;i<dim;i++)
	{
		N[i]=T.N[i];
		yuji[i]=T.yuji[i];
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
	NUM=new long[dim];
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
	if (num!=dim)
	{
		cout<<"dim must be equal";
		delete[] NUM;
		return;
	}

	long index=0;
	for (long i=0;i<dim;i++)
	{
		index+=yuji[i]*NUM[i];
	}
	delete[] NUM;
	data[index]=a;

}

double Tensor::getElement(long b,...)
{
	long num=0;
	long *NUM;
	long *p=&b;

	while(*p!=0)
	{
		num++;
		p++;
	}

	if (dim!=num)
	{
		cout<<"dim wrong"<<endl;
		return 0;
	}

	NUM=new long[num];
	int i=0;
	p=&b;
	while(*p!=0)
	{
		NUM[i]=*p;
		p++;
		i++;
	}
	long index=0;
	for (int i=0;i<dim;i++)
	{
		index+=NUM[i]*yuji[i];
	}

	delete[] NUM;
	return data[index];
}

N_SVD *Tensor::m_mode_svd()
{
	N_SVD *nsvd;

	return nsvd;
}
