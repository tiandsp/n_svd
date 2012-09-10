#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

typedef struct _DATA
{
	vector< double* > U;
	vector< vector<double> > core_tensor;
}Data;

typedef struct _SVD_Data
{
	vector< vector<double> > u;
	vector< vector<double> > s;
	vector< vector<double> > v;
}SVD_Data;

typedef struct _Size_Vec
{
	double r;	//二维vector的行的大小
	double c;	//二维vector的列的大小
}Size_Vec;

Size_Vec size(vector< vector<double> > M)	//返回一个二维vector行与列的大小
{
	Size_Vec s;
	vector< vector<double> >::iterator p;
	s.r=M.size();
	p=M.begin();
	s.c=p->size();

	return s;
}


double *vector2matrix(vector< vector<double> > M)  //将二维的vector转换为二维的数组
{
	vector< vector<double> >::iterator p_M;
	vector<double>::iterator pp_M;

	double *m;
	double r,c;

	r=M.size();
	p_M=M.begin();
	c=p_M->size();
	m=new double[long(r*c)];

	double i,j;
	for (i=0,p_M=M.begin();p_M!=M.end();i++,p_M++)
	{
		for (j=0,pp_M=p_M->begin();pp_M!=p_M->end();j++,pp_M++)
		{
			m[long(i*c+j)]=*pp_M;
		}
	}

	return m;

}

vector< vector<double> > matrix2vector(double *m,double r,double c)  //将二维的数组转换为二维的vector
{
	vector< vector<double> > M;
	vector<double> data;

	for (double i=0;i<r;i++)
	{
		for (double j=0;j<c;j++)
		{
			data.push_back(m[long(i*c+j)]);
		}
		M.push_back(data);
		data.clear();
	}

	return M;
}

double *inv(vector< vector<double> > s)
{
	double *re;

	return re;
}

vector< vector<double> > shiftdim(vector< vector<double> > tensor,double n)
{
	vector< vector<double> > re;
	vector< vector<double> >::iterator p_re;
	vector<double>::iterator pp_re;
	vector<double> data;


	return re;
}

vector< vector<double> > reshape(vector< vector<double> > tensor,double r, double c)
{
	vector< vector<double> > re;

	return re;

}

vector< vector<double> > matricize(vector< vector<double> > tensor,double n)
{

	double ss;		//张量tensor中第n维的大小
	double dim;
	vector<double> s;		//张量每一维的大小，是一个一维数组
	vector<double>::iterator p_s;	//用来迭代s这个数组
	vector< vector<double> >::iterator p_tensor;

	vector< vector<double> > shift_tensor; //对原tensor进行n的移位，以使分解中要保持的维数在第一维
	vector< vector<double> > re_tensor;	//返回的结果，这里是个二维矩阵
	vector<double> tmp;

	for (p_tensor=tensor.begin();p_tensor!=tensor.end();p_tensor++)
	{
		s.push_back(p_tensor->size());		//得到tensor每一维的大小
	}

	p_tensor=tensor.begin()+n;	//指向tensor展开中保持不变的那一维
	ss=p_tensor->size();		//所指向维数的大小

	dim=1;				//tensor中所有维数的大小的乘积
	for (p_s=s.begin();p_s!=s.end();p_s++)
	{
		dim *= (*p_s);
	}


	shift_tensor=shiftdim(tensor,n-1);
	re_tensor=reshape(shift_tensor,ss,dim/ss);

	/*
	for (p_tensor=tensor.begin();p_tensor!=tensor.end();p_tensor++)
	{
		for (p_s=p_tensor->begin();p_s!=p_tensor->end();p_s++)
		{
			tmp.push_back(*p_s);
		}
	}

	p_s=tmp.begin();
	s.clear();
	for (double i=0;i<ss;i++)
	{
		for (double j=0;j<dim/ss;j++)
		{
			s.push_back(*p_s);
			p_s++;
		}
		re_tensor.push_back(s);
		s.clear();
	}
	*/
	return re_tensor;

}

SVD_Data svd(vector< vector<double> > M,int n)
{

}

SVD_Data svd(vector< vector<double> > M)
{

}



double *mat_multiplier(double *a,double *b,double a_r,double a_c,double b_r,double b_c)
{
	double *c;
	c=new double[long(a_r*b_c)];
	memset(c,0,a_r*b_c*sizeof(double));

	for (double i=0;i<a_r;i++)
	{
		for (double j=0;j<a_c;j++)
		{
			for (double k=0;k<b_c;k++)
			{
				c[long(i*b_c+k)]=c[long(i*b_c+k)]+a[long(i*a_c+j)]*b[long(j*b_c+k)];
			}
		}
	}
	return c;
}



vector< vector<double> > vec_multiplier(double *a,double*b,double a_r,double a_c,double b_r,double b_c)
{
	vector< vector<double> > re;
	double *c;
	c=new double[long(a_r*b_c)];

	c=mat_multiplier(a,b,a_r,a_c,b_r,b_c);
	re=matrix2vector(c,a_r,b_c);
	delete[] c;

	return re;
}

vector< vector<double> > abs_sqrt_dot(vector< vector<double> > a,vector< vector<double> > b)
{
	vector< vector<double> > re;
	vector<double> data;
	vector< vector<double> >::iterator p_a;
	vector< vector<double> >::iterator p_b;
	vector<double>::iterator pp_a;
	vector<double>::iterator pp_b;

	for (p_a=a.begin(),p_b=b.begin();
		p_a!=a.end(),p_b!=b.end();p_a++,p_b++)
	{
		for (pp_a=p_a->begin(),pp_b=p_b->begin();
			pp_a!=p_a->end(),pp_b!=p_b->end();pp_a++,pp_b++)
		{
			data.push_back(abs(sqrt((*pp_a)*(*pp_b))));
		}
		re.push_back(data);
		data.clear();
	}
	return re;
}

vector< vector<double> > tfastsvd(vector< vector<double> > M)
{
	vector< vector<double> >::iterator p_M;
	vector<double>::iterator pp_M;
	vector< vector<double> > M2;
	SVD_Data re_svd;
	double r,c;
	double u;

	r=M.size();
	p_M=M.begin();
	c=p_M->size();

	if (r==1)
	{
		u=1;
		re_svd=svd(M,0);
	}
	else if (r<c)
	{
		double *m1;
		double *m2;
		m1=new double[long(r*c)];
		m2=new double[long(c*r)];

		m1=vector2matrix(M);
		for (double j=0;j<c;j++)
		{
			for (double i=0;i<r;i++)
			{
				m2[long(j*r+i)]=m1[long(i*c+j)];
			}
		}

		M2=vec_multiplier(m1,m2,r,c,c,r);
		delete[] m1;
		delete[] m2;

		re_svd=svd(M2);
	}
	else if (c==1)
	{
		vector< vector<double> > tmp_m;
		vector<double> data;
		data.push_back(1);
		re_svd.v.push_back(data);

		tmp_m=abs_sqrt_dot(M,M);
		re_svd.s=tmp_m;

		double *m;
		double *v;
		double *s;
		double *u;
		m=new double[long(r*c)];
		v=new double;
		s=new double[long(r*c)];
		u=new double[long()];
		m=vector2matrix(M);
		v=vector2matrix(re_svd.v);
		s=inv(re_svd.s);



		re_svd=svd(M,0);

		delete[] m;
		delete v;
		delete[] s;
		delete[] u;
	}
	else if (c<r)
	{

	}
	else
	{
		re_svd=svd(M,0);
	}
	return re_svd.u;

}

vector< vector<double> > mode_m_prod(vector< vector<double> > cdata,double *u,double i,int flag)
{


}


Data m_mode_svd(vector< vector<double> > dtensor,vector<double> modes)
{
	Data re;			//函数返回的结果，结构体
	vector<double> d;	//detnsor张量每一维的大小
	vector<double>::iterator p_modes;	//指向传入的要分解的维的迭代器
	vector<double>::iterator p_d;		//指向d,迭代每一维的detensor
	vector< vector<double> >::iterator pit_dtensor;		//对dtensor的列进行迭代
	vector< vector<double> > cdata;		//dtensor的一个副本，主要对他进行处理


	for (pit_dtensor=dtensor.begin();pit_dtensor!=dtensor.end();pit_dtensor++)	//统计dtensor每一维的大小，存入d中	
	{
		d.push_back(pit_dtensor->size());		
	}

	cdata=dtensor;
	vector< vector<double> > M;		//按dtensor中第i维展开后的结果，是一个二维数组
	Size_Vec size_M;		//M行与列分别的大小
	vector< vector<double> >::iterator p_M;		//对M的列进行迭代的迭代器
	vector<double>::iterator pp_M;				//对M的行进行迭代的迭代器
	vector< vector<double> > u;					//相应的按第i维展开后的特征，是一个二维数组
	double i;								//用来跟踪d的每个值

	for (i=0,p_modes=modes.begin();p_modes!=modes.end();i++,p_modes++)
	{
		M=matricize(cdata,*p_modes);		//对cdata张量按其第i维展开，比如cdata为（3,4,5,3,2）的张量，按第三维展开M就是（5,72）的矩阵

		for (p_M=M.begin();p_M!=M.begin();p_M++)
		{
			double tmp;
			double i;
			tmp=0;
			for (i=0,pp_M=p_M->begin();pp_M!=p_M->end();i++,pp_M++)		//求M中每一行的平均值
			{
				tmp=tmp*i;
				tmp=(tmp+(*pp_M))/(i+1);
			}

			for (pp_M=p_M->begin();pp_M!=p_M->end();pp_M++)		//这一行的所有的数再减去这个平均值
			{
				*pp_M=*pp_M-tmp;
			}
		}


		u=tfastsvd(M);	    //快速svd分解求出每一维dtensor的特征值
		Size_Vec s_u;
		double *trans_u;		//为压入到re.U做准备
		s_u=size(u);		//得到的特征矩阵的行列的大小
		trans_u=new double[long(s_u.r*s_u.c)];
		trans_u=vector2matrix(u);

		re.U.push_back(trans_u);	//存储按每一维展开的特征矩阵
		delete[] trans_u;

		double threld;
		threld=1;
		for (p_d=d.begin();p_d!=d.end();p_d++)
		{
			if (p_d!=p_modes)
			{
				threld=threld*(*p_d);
			}
		}

		if (threld>350000)
		{
			cdata=mode_m_prod(cdata,re.U.at(i),i,1);
		} 
		else
		{
			cdata=mode_m_prod(cdata,re.U.at(i),i,0);
		}

	}

	return re;
}


int main()
{
	//	Data re;




	return 0;
}