#include "n_svd.h"
#include <iostream>
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
	double r;	//��άvector���еĴ�С
	double c;	//��άvector���еĴ�С
}Size_Vec;

Size_Vec size(vector< vector<double> > M)	//����һ����άvector�����еĴ�С
{
	Size_Vec s;
	vector< vector<double> >::iterator p;
	s.r = M.size();
	p = M.begin();
	s.c = p -> size();

	return s;
}

double *vector2matrix(vector< vector<double> > M)  //����ά��vectorת��Ϊ��ά������
{
	vector< vector<double> >::iterator p_M;
	vector<double>::iterator pp_M;

	double *m;
	double r, c;

	r = M.size();
	p_M = M.begin();
	c = p_M -> size();
	m = new double[long(r * c)];

	double i, j;
	for (i = 0, p_M = M.begin(); p_M != M.end(); i++, p_M++)
	{
		for (j = 0, pp_M = p_M -> begin(); pp_M != p_M -> end(); j++, pp_M++)
		{
			m[long(i * c + j)] = *pp_M;
		}
	}

	return m;
}

double **vector2matrix2(vector< vector<double> > M)
{
	double **m;
	Size_Vec size_M;
	vector< vector<double> >::iterator p_M;
	vector<double>::iterator pp_M;	
	size_M = size(M);

	m = new double *[long(size_M.r)];
	for (double i = 0; i < size_M.r; i++)
	{
		m[long(i)] = new double[long(size_M.c)];
	}

	double i, j;
	for (i = 0, p_M = M.begin(); p_M != M.end(); i++, p_M++)
	{
		for (j = 0, pp_M = p_M -> begin(); pp_M != p_M -> end(); j++, pp_M++)
		{
			m[long(i)][long(j)] = *pp_M;
		}
	}

	return m;

}

vector< vector<double> > matrix2vector(double *m, double r, double c)  //����ά������ת��Ϊ��ά��vector
{
	vector< vector<double> > M;
	vector<double> data;

	for (double i = 0; i < r; i++)
	{
		for (double j = 0; j < c; j++)
		{
			data.push_back(m[long(i * c + j)]);
		}
		M.push_back(data);
		data.clear();
	}

	return M;
}

vector< vector<double> > matrix2vector2(double **m,double r,double c)
{
	vector< vector<double> > M;
	vector<double> data;
	for (double i = 0 ; i < r; i++)
	{
		for (double j = 0; j < c; j++)
		{
			data.push_back(m[long(i)][long(j)]);
		}
		M.push_back(data);
		data.clear();

	}
	return M;

}

void brinv(double *a, double n)
{ 
	double *is, *js, i, j, k, l, u, v;
	double d,p;
	is = new double[long(n)];
	js = new double[long(n)];
	for (k = 0; k < n; k++)
	{ 
		d = 0.0;
		for (i = k; i < n; i++)
			for (j = k; j < n; j++)
			{	
				l = i * n + j; 
				p = fabs(a[long(l)]);
				if (p > d) 
				{ 
					d = p;
					is[long(k)] = i; 
					js[long(k)] = j;
				}
			}
			if (d + 1.0 == 1.0)
			{
				delete[] is;
				delete[] js;
				printf("err**not inv\n");
				cout<<"error not inv"<<endl;
				exit(1);
			}
			if (is[long(k)] != k)
				for (j = 0; j < n; j++)
				{
					u = k * n + j;
					v = is[long(k)] * n + j;
					p = a[long(u)];
					a[long(u)] = a[long(v)];
					a[long(v)] = p;
				}
				if (js[long(k)] != k)
					for (i = 0; i < n; i++)
					{
						u = i * n + k;
						v = i * n + js[long(k)];
						p = a[long(u)];
						a[long(u)] = a[long(v)];
						a[long(v)] = p;
					}
					l = k * n + k;
					a[long(l)] = 1.0 / a[long(l)];
					for (j = 0; j < n; j++)
						if (j != k)
						{ 
							u = k * n + j; 
							a[long(u)] = a[long(u)] * a[long(l)];
						}
						for (i = 0; i < n; i++)
							if (i != k)
								for (j = 0; j < n; j++)
									if (j != k)
									{
										u = i * n + j;
										a[long(u)] = a[long(u)] - a[long(i * n + k)] * a[long(k * n + j)];
									}
									for (i = 0; i < n; i++)
										if (i != k)
										{
											u = i * n + k; 
											a[long(u)] = -a[long(u)] * a[long(l)];
										}
	}
	for (k = n-1; k >= 0; k--)
	{
		if (js[long(k)] != k)
			for (j = 0; j < n; j++)
			{
				u = k * n + j;
				v = js[long(k)] * n + j;
				p = a[long(u)]; 
				a[long(u)] = a[long(v)];
				a[long(v)] = p;
			}
			if (is[long(k)] != k)
				for (i = 0; i < n; i++)
				{
					u = i * n + k; 
					v = i * n + is[long(k)];
					p = a[long(u)]; 
					a[long(u)] = a[long(v)];
					a[long(v)] = p;
				}
	}

	delete[] is;
	delete[] js;

}

double *inv(vector< vector<double> > s)
{
	double *re;
	double *data;
	Size_Vec size_s;
	size_s = size(s);

	re = new double[long(size_s.r * size_s.c)];
	data = new double[long(size_s.r * size_s.c)];
	data = vector2matrix(s);

	brinv(data, size_s.r);

	return data;
}

void reverse_vec(vector< vector<double> > tensor, double p, double q)
{
	vector<double> tmp;
	while(p < q)
	{
		tmp = tensor.at(long(p));
		tensor.at(long(p)) = tensor.at(long(q));
		tensor.at(long(q)) = tmp;
		p++;
		q--;
	}
}

vector< vector<double> > shiftdim(vector< vector<double> > tensor, double n)		//ʹ�������Ҫ�ֽ����һά���ڵ�һά
{
	double len;
	len = tensor.size();
	
	reverse_vec(tensor, 0, len-n-1);			//һ��С�ı任���ɣ��������������֮��
	reverse_vec(tensor, len-n, len-1);
	reverse_vec(tensor, 0, len-1);

	return tensor;
}

vector< vector<double> > reshape(vector< vector<double> > tensor, double r, double c)	//����ά��tensor��Ϊ��ά��
{
	vector< vector<double> > re;
	vector< vector<double> >::iterator p_tensor;
	vector<double>::iterator pp_tensor;
	vector<double> tmp;
	vector<double> data;
	vector<double>::iterator p_tmp;

	for (p_tensor = tensor.begin(); p_tensor != tensor.end(); p_tensor++)
	{
		for (pp_tensor = p_tensor->begin(); pp_tensor != p_tensor->end(); pp_tensor++)
		{
			tmp.push_back(*pp_tensor);
		}
	}

	p_tmp = tmp.begin();
	for (double i = 0; i < r; i++)
	{
		for (double j = 0; j < c; j++)
		{
			data.push_back(*p_tmp);
			p_tmp++;
		}
		re.push_back(data);
		data.clear();
	}
	return re;

}

vector< vector<double> > matricize(vector< vector<double> > tensor, double n)
{

	double ss;		//����tensor�е�nά�Ĵ�С
	double dim;
	vector<double> s;		//����ÿһά�Ĵ�С����һ��һά����
	vector<double>::iterator p_s;	//��������s�������
	vector< vector<double> >::iterator p_tensor;

	vector< vector<double> > shift_tensor; //��ԭtensor����n����λ����ʹ�ֽ���Ҫ���ֵ�ά���ڵ�һά
	vector< vector<double> > re_tensor;	//���صĽ���������Ǹ���ά����
	vector<double> tmp;

	for (p_tensor = tensor.begin(); p_tensor != tensor.end(); p_tensor++)
	{
		s.push_back(p_tensor -> size());		//�õ�tensorÿһά�Ĵ�С
	}

	p_tensor = tensor.begin() + n;	//ָ��tensorչ���б��ֲ������һά
	ss = p_tensor -> size();		//��ָ��ά���Ĵ�С

	dim = 1;				//tensor������ά���Ĵ�С�ĳ˻�
	for (p_s = s.begin(); p_s != s.end(); p_s++)
	{
		dim *= (*p_s);
	}


	shift_tensor = shiftdim(tensor, n);
	re_tensor = reshape(shift_tensor, ss, dim / ss);

	return re_tensor;

}

SVD_Data svd(vector< vector<double> > M, int n)
{
	SVD_Data re;



	
	return re;
}

SVD_Data svd(vector< vector<double> > M)
{
	SVD_Data re;	
	Size_Vec size_M;
	double **a;
	double **p;
	double *d;
	double **q;

	size_M = size(M);
	if (size_M.r >= size_M.c)
	{
		a=new double *[long(size_M.r)];
		for (double i = 0; i < size_M.r; i++)
		{
			a[long(i)]=new double[long(size_M.c)];
		}

		p=new double *[long(size_M.r)];
		for (double i=0;i< size_M.r;i++)
		{
			p[long(i)]=new double[long(size_M.r)];
		}

		d=new double[long(size_M.c)];

		q=new double *[long(size_M.c)];
		for (double i=0;i<size_M.c;i++)
		{
			q[long(i)]=new double[long(size_M.c)];
		}

		a=vector2matrix2(M);

//		svd(size_M.r, size_M.c, a, p, d, q);

		re.s=matrix2vector2(p,size_M.r,size_M.r);
		re.u=matrix2vector(d,1,0);
		re.v=matrix2vector2(q,size_M.c,size_M.c);

		for (double i=0;i<size_M.r;i++)
		{
			delete[] a[long(i)];
		}
		delete[] a;

		for (double i=0;i<size_M.r;i++)
		{
			delete[] p[long(i)];
		}
		delete[] p;

		delete[] d;

		for (double i=0;i<size_M.c;i++)
		{
			delete[] q[long(i)];
		}
		delete[] q;
	}
	else if (size_M.r < size_M.c)
	{


	}



	return re;

}



double *mat_multiplier(double *a, double *b, double a_r, double a_c, double b_r, double b_c)
{
	double *c;
	c=new double[long(a_r * b_c)];
	memset(c, 0, a_r * b_c * sizeof(double));

	for (double i = 0; i < a_r; i++)
	{
		for (double j = 0; j < a_c; j++)
		{
			for (double k = 0; k < b_c; k++)
			{
				c[long(i * b_c + k)] = c[long(i * b_c + k)] + a[long(i * a_c + j)] * b[long(j * b_c + k)];
			}
		}
	}
	return c;
}



vector< vector<double> > vec_multiplier(double *a, double*b, double a_r, double a_c, double b_r, double b_c)
{
	vector< vector<double> > re;
	double *c;
	c = new double[long(a_r * b_c)];

	c = mat_multiplier(a, b, a_r, a_c, b_r, b_c);
	re = matrix2vector(c, a_r, b_c);
	delete[] c;

	return re;
}

vector< vector<double> > abs_sqrt_dot(vector< vector<double> > a, vector< vector<double> > b)
{
	vector< vector<double> > re;
	vector<double> data;
	vector< vector<double> >::iterator p_a;
	vector< vector<double> >::iterator p_b;
	vector<double>::iterator pp_a;
	vector<double>::iterator pp_b;

	for (p_a = a.begin(), p_b = b.begin();
		p_a != a.end(), p_b != b.end(); p_a++, p_b++)
	{
		for (pp_a = p_a -> begin(), pp_b = p_b -> begin();
			pp_a != p_a -> end(), pp_b != p_b -> end(); pp_a++,pp_b++)
		{
			data.push_back(abs(sqrt((*pp_a) * (*pp_b))));
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

	SVD_Data re_svd;
	Size_Vec size_M;
	double u;

	size_M = size(M);

	if (size_M.r == 1)
	{
		u = 1;
		re_svd = svd(M, 0);
	}
	else if (size_M.r < size_M.c)
	{
		vector< vector<double> > M2;
		double *m1;
		double *m2;
		m1 = new double[long(size_M.r * size_M.c)];
		m2 = new double[long(size_M.c * size_M.r)];

		m1 = vector2matrix(M);
		for (double j = 0; j < size_M.c; j++)
		{
			for (double i = 0; i < size_M.r; i++)
			{
				m2[long(j * size_M.r + i)] = m1[long(i * size_M.c + j)];
			}
		}

		M2 = vec_multiplier(m1, m2, size_M.r, size_M.c, size_M.c, size_M.r);
		delete[] m1;
		delete[] m2;

		re_svd=svd(M2);
	}
	else if (size_M.c == 1)
	{
		vector< vector<double> > tmp_m;
		vector<double> data;
		data.push_back(1);
		re_svd.v.push_back(data);

		tmp_m = abs_sqrt_dot(M, M);
		re_svd.s = tmp_m;

		double *m;
		double *v;
		double *s;
		double *u;
		m = new double[long(size_M.r * size_M.c)];
		v = new double;
		s = new double[long(size_M.r * size_M.c)];
		u = new double[long()];
		m = vector2matrix(M);
		v = vector2matrix(re_svd.v);
		s = inv(re_svd.s);
		re_svd = svd(M, 0);

		delete[] m;
		delete v;
		delete[] s;
		delete[] u;
	}
	else if (size_M.c < size_M.r)
	{

	}
	else
	{
		re_svd = svd(M, 0);
	}
	return re_svd.u;

}

vector< vector<double> > mode_m_prod(vector< vector<double> > cdata, double *u, double i, int flag)
{
	vector< vector<double> > re;

	return re;
}


Data m_mode_svd(vector< vector<double> > dtensor, vector<double> modes)
{
	Data re;			//�������صĽ�����ṹ��
	vector<double> d;	//detnsor����ÿһά�Ĵ�С
	vector<double>::iterator p_modes;	//ָ�����Ҫ�ֽ��ά�ĵ�����
	vector<double>::iterator p_d;		//ָ��d,����ÿһά��detensor
	vector< vector<double> >::iterator pit_dtensor;		//��dtensor���н��е���
	vector< vector<double> > cdata;		//dtensor��һ����������Ҫ�������д���


	for (pit_dtensor = dtensor.begin(); pit_dtensor != dtensor.end(); pit_dtensor++)	//ͳ��dtensorÿһά�Ĵ�С������d��	
	{
		d.push_back(pit_dtensor -> size());		
	}

	cdata = dtensor;
	vector< vector<double> > M;		//��dtensor�е�iάչ����Ľ������һ����ά����
	Size_Vec size_M;		//M�����зֱ�Ĵ�С
	vector< vector<double> >::iterator p_M;		//��M���н��е����ĵ�����
	vector<double>::iterator pp_M;				//��M���н��е����ĵ�����
	vector< vector<double> > u;					//��Ӧ�İ���iάչ�������������һ����ά����
	double i;								//��������d��ÿ��ֵ

	for (i = 0, p_modes = modes.begin(); p_modes != modes.end(); i++, p_modes++)
	{
		M = matricize(cdata, *p_modes);		//��cdata���������iάչ��������cdataΪ��3,4,5,3,2����������������άչ��M���ǣ�5,72���ľ���

		for (p_M = M.begin(); p_M != M.begin(); p_M++)
		{
			double tmp;
			double i;
			tmp = 0;
			for (i = 0, pp_M = p_M -> begin(); pp_M != p_M -> end(); i++, pp_M++)		//��M��ÿһ�е�ƽ��ֵ
			{
				tmp = tmp * i;
				tmp = (tmp + (*pp_M)) / (i + 1);
			}

			for (pp_M = p_M -> begin(); pp_M != p_M -> end(); pp_M++)		//��һ�е����е����ټ�ȥ���ƽ��ֵ
			{
				*pp_M = *pp_M - tmp;
			}
		}


		u = tfastsvd(M);	    //����svd�ֽ����ÿһάdtensor������ֵ
		Size_Vec s_u;
		double *trans_u;		//Ϊѹ�뵽re.U��׼��
		s_u = size(u);		//�õ���������������еĴ�С
		trans_u = new double[long(s_u.r * s_u.c)];
		trans_u = vector2matrix(u);

		re.U.push_back(trans_u);	//�洢��ÿһάչ������������
		delete[] trans_u;

		double threld;
		threld = 1;
		for (p_d = d.begin(); p_d != d.end(); p_d++)
		{
			if (p_d != p_modes)
			{
				threld = threld * (*p_d);
			}
		}

		if (threld > 350000)
		{
			cdata = mode_m_prod(cdata, re.U.at(i), i, 1);
		} 
		else
		{
			cdata = mode_m_prod(cdata, re.U.at(i), i, 0);
		}

	}

	return re;
}
