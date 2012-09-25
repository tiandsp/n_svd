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

Tensor::Tensor(Mat a)
{

	n=1;
	dim=a.getSize();

	N=new long[dim];
	yuji=new long[dim];
	long i=0;
	for (long i=0;i<dim;i++)
	{
		N[i]=a.getElement(0,i);
		n*=a.getElement(0,i);
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

long Tensor::getSize(long i)
{
	return N[i];
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

	if (num==1)
	{
		data[b]=a;
		delete[] NUM;
		return;
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

void Tensor::setElement(double a,Mat M)
{
	long num=0;
	num=M.getSize();
	if (num!=dim)
	{
		cout<<"dim must be equal";
		return;
	}

	for (long i=0;i<num;i++)
	{
		if (M.getElement(0,i)>N[i])
		{
			cout<<"setDlement wrong";
			return;
		}
	}

	long index=0;
	for (long i=0;i<dim;i++)
	{
		index+=yuji[i]*M.getElement(0,i);
	}
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

	if (num==1)
	{
		return data[b];
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

double Tensor::getElement(Mat M)
{
	long num=0;
	num=M.getSize();

	if (dim!=num)
	{
		cout<<"dim wrong"<<endl;
		return 0;
	}

	for (long i=0;i<num;i++)
	{
		if (M.getElement(0,i)>N[i])
		{
			cout<<"setDlement wrong";
			return 0;
		}
	}

	long index=0;
	for (int i=0;i<dim;i++)
	{
		index+=M.getElement(0,i)*yuji[i];
	}
	return data[index];
}


N_SVD *Tensor::m_mode_svd()
{
	N_SVD *nsvd;
	nsvd=new N_SVD;
	Mat d;
	Mat M;
	d=this->getSize();

	Tensor cdata(*this);
	long modes=d.getSize();
	nsvd->U=new Mat[modes];

	for (long i=0;i<modes;i++)
	{
		M=cdata.matricize(i);
		long r=M.getRow();
		long c=M.getCol();
		Mat tmp;
		tmp=M.mean(2);
		tmp=tmp.repmat(1,c);
		M=M-tmp;

		nsvd->U[i]=M.tfastsvd();
	
		cdata=cdata.mode_m_prod(M,i);

	}
	nsvd->core=new Tensor(cdata);

	return nsvd;
}

Tensor Tensor::shiftdim(long i)
{



	return *this;
}

Mat Tensor::reshape(long r,long c)
{
	Mat re;
	long k=0;

	for (long i=0;i<r;i++)
	{
		for (long j=0;j<c;j++)
		{
			re.setElement(data[k],r,c);
			k++;
		}
	}
	return re;
}

Tensor Tensor::reshape(Mat M)
{
	Tensor tmp(*this);




	return tmp;
}

Mat Tensor::matricize(long i)
{
	Mat M;
	Tensor cdata(*this);
	Mat s=cdata.getSize();
	long ss=cdata.getSize(i);
	long dim=s.prod();

	cdata=cdata.shiftdim(i);

	M=cdata.reshape(ss,dim/ss);
	return M;
}

Tensor Tensor::mode_m_prod(Mat M,long n)
{
	Mat dims(1,dim);
	long maxn=getDim();
	long r=M.getRow();
	long c=M.getCol();

	for (long i=0;i<dim;i++)
	{
		if (n==i)
			dims.setElement(r,0,i);
		else
			dims.setElement(N[i],0,i);
	}

	Mat t;
	t=matricize(n);
	t=t*M;
	Tensor *pt;
	Tensor re(*this);
	pt=new Tensor(*this);
	
	pt=M.tensorize(n,dims);
	re=*pt;
	delete pt;

	return re;
}
