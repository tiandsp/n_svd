#include "Mat.h"
#include <memory>
#include <iostream>
#include <cmath>
#include <stdlib.h>

//svd function
#define SIGN(u, v)     ( (v)>=0.0 ? fabs(u) : -fabs(u) )
#define MAX(x, y)     ( (x) >= (y) ? (x) : (y) )  

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

	Svd.u=new double *[row];
	for (long i=0;i<row;i++)
	{
		Svd.u[i]=new double[row];
		memset(Svd.u[i],0,row*sizeof(double));
	}

	Svd.s=new double[col];
	memset(Svd.s,0,col*sizeof(double));

	Svd.v=new double *[col];
	for (long i=0;i<col;i++)
	{
		Svd.v[i]=new double[col];
		memset(Svd.v[i],0,col*sizeof(double));
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

	Svd.u=new double *[row];
	for (long i=0;i<row;i++)
	{
		Svd.u[i]=new double[row];
	}
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<row;j++)
		{
			Svd.u[i][j]=M.Svd.u[i][j];
		}
	}

	Svd.s=new double[col];
	for (long i=0;i<col;i++)
	{
		Svd.s[i]=M.Svd.s[i];
	}
	
	Svd.v=new double *[col];
	for (long i=0;i<col;i++)
	{
		Svd.v[i]=new double[col];
	}
	for (long i=0;i<col;i++)
	{
		for (long j=0;j<col;j++)
		{
			Svd.v[i][j]=M.Svd.v[i][j];
		}
	}
}

Mat::~Mat()
{
	for (long i=0;i<row;i++)
	{
		delete[] data[i];
		delete[] Svd.u[i];
	}
	delete[] data;
	delete[] Svd.u;

	delete[] Svd.s;

	for (long i=0;i<col;i++)
	{
		delete[] Svd.v[i];
	}
	delete Svd.v;

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
		delete[] Svd.u[i];
	}
	delete[] data;
	delete[] Svd.u;

	delete[] Svd.s;

	for (long i=0;i<col;i++)
	{
		delete[] Svd.v[i];
	}
	delete Svd.v;

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

	Svd.u=new double *[row];
	for (long i=0;i<row;i++)
	{
		Svd.u[i]=new double[row];
	}
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<row;j++)
		{
			Svd.u[i][j]=M.Svd.u[i][j];
		}
	}

	Svd.s=new double[col];
	for (long i=0;i<col;i++)
	{
		Svd.s[i]=M.Svd.s[i];
	}

	Svd.v=new double *[col];
	for (long i=0;i<col;i++)
	{
		Svd.v[i]=new double[col];
	}
	for (long i=0;i<col;i++)
	{
		for (long j=0;j<col;j++)
		{
			Svd.v[i][j]=M.Svd.v[i][j];
		}
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

	for (long i=0;i<row;i++)
	{
		delete[] Svd.u[i];
	}
	delete[] Svd.u;

	delete[] Svd.s;

	for (long i=0;i<col;i++)
	{
		delete[] Svd.v[i];
	}
	delete Svd.v;

	for (long j=0;j<col;j++)
	{
		for (long i=0;i<row;i++)
		{
			tmp.data[j][i]=data[i][j];
		}
	}

	Svd.u=new double *[col];
	for (long i=0;i<col;i++)
	{
		Svd.u[i]=new double[col];
		memset(Svd.u[i],0,col*sizeof(double));
	}

	Svd.s=new double[row];
	memset(Svd.s,0,row*sizeof(double));

	Svd.v=new double *[row];
	for (long i=0;i<row;i++)
	{
		Svd.v[i]=new double[row];
		memset(Svd.v[i],0,row*sizeof(double));
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
			tmp.data[i][j] =data[i][j]+ M.data[i][j];
		}
	}
	return tmp;
}

Mat Mat::operator+(double a)
{
	Mat tmp(row,col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j] =data[i][j]+ a;
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
			tmp.data[i][j] =data[i][j]-M.data[i][j];
		}
	}
	return tmp;
}

Mat Mat::operator-(double a)
{
	Mat tmp(row,col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j] =data[i][j]-a;
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

Mat Mat::operator*(double a)
{
	Mat tmp(row,col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j] =data[i][j]*a;
		}
	}
	return tmp;
}

Mat Mat::operator/(Mat &M)
{
	if (col!=M.row)
	{
		cout<<"row or col not match"<<endl;
		return *this;
	}

	Mat tmp(row,M.col);
	M.inv();

	tmp=*this * M;
	return tmp;
}

Mat Mat::operator/(double a)
{
	if (a==0)
	{
		cout<<"dive zero wrong"<<endl;
		system("pause");
		return *this;
	}

	Mat tmp(row,col);

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			tmp.data[i][j] =data[i][j]/a;
		}
	}
	return tmp;
}


void Mat::dotMultiplication(Mat &M)
{
	if (row!=M.row || col!=M.col)
	{
		cout<<"row or col not match"<<endl;
		return ;
	}

	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			data[i][j]=data[i][j]*M.data[i][j];
		}
	}
}

void Mat::dotDivision(Mat &M)
{
	if (row!=M.row || col!=M.col)
	{
		cout<<"row or col not match"<<endl;
		return ;
	}
	
	for (long i=0;i<row;i++)
	{
		for (long j=0;j<col;j++)
		{
			data[i][j]=data[i][j]/M.data[i][j];
		}
	}
}

void Mat::inv()
{
	if (col!=row)
	{
		cout<<"col must equal row"<<endl;
		system("pause");
		return ;
	}
	
	long n=row;
	double *is, *js, i, j, k, l, u, v;
	double d,p;
	double *a;
	a=new double[n*n];

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			a[long(i*n+j)]=data[long(i)][long(j)];
		}
	}


	is = new double[n];
	js = new double[n];
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
				cout<<"error not inv"<<endl;
				system("pause");
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

	for (i=0;i<n;i++)
	{
		for (j=0;j<n;j++)
		{
			data[long(i)][long(j)]=a[long(i*n+j)];
		}
	}
	delete[] is;
	delete[] js;
	delete[] a;
}

void Mat::svd()
{

	if (row>=col)
	{
		svd(row,col,data,Svd.u,Svd.s,Svd.v);
	}
	else
	{
		Mat tmp(*this);
		tmp=~tmp;
		svd(tmp.row,tmp.col,tmp.data,tmp.Svd.u,tmp.Svd.s,tmp.Svd.v);

		for (long i=0;i<row;i++)
		{
			for (long j=0;j<row;j++)
			{
				Svd.u[i][j]=-tmp.Svd.v[i][j];
			}
		}

		for (long i=0;i<col;i++)
		{
			Svd.s[i]=tmp.Svd.s[i];
		}

		for (long i=0;i<col;i++)
		{
			for (long j=0;j<col;j++)
			{
				Svd.v[i][j]=-tmp.Svd.u[i][j];
			}
		}
	}

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
	cout<<endl;
}


static double   radius(double u, double v)
{
        double          w;
        u = fabs(u);
        v = fabs(v);
        if (u > v) {
                w = v / u;
                return (u * sqrt(1. + w * w));
    } else {
        if (v) {
                 w = u / v;
             return (v * sqrt(1. + w * w));
                } else
                        return 0.0;
        }
}

/*
 Given matrix a[m][n], m>=n, using svd decomposition a = p d q' to get
 p[m][n], diag d[n] and q[n][n].
*/
void Mat::svd(int m, int n, double **a, double **p, double *d, double **q)
{
        int             flag, i, its, j, jj, k, l, nm, nm1 = n - 1, mm1 = m - 1;
        double          c, f, h, s, x, y, z;
        double          anorm = 0, g = 0, scale = 0;
        //double         *r = tvector_alloc(0, n, double);
        double          *r = (double*)malloc(sizeof(double)*n);

        for (i = 0; i < m; i++)
                for (j = 0; j < n; j++)
                        p[i][j] = a[i][j];
        //for (i = m; i < n; i++)
        //                p[i][j] = 0;

        /* Householder reduction to bidigonal form */
        for (i = 0; i < n; i++)
        {
                l = i + 1;
                r[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m)
                {
                        for (k = i; k < m; k++)
                                scale += fabs(p[k][i]);
                        if (scale)
                        {
                                for (k = i; k < m; k++)
                                {
                                        p[k][i] /= scale;
                                        s += p[k][i] * p[k][i];
                                }
                                f = p[i][i];
                                g = -SIGN(sqrt(s), f);
                                h = f * g - s;
                                p[i][i] = f - g;
                                if (i != nm1)
                                {
                                        for (j = l; j < n; j++)
                                        {
                                                for (s = 0.0, k = i; k < m; k++)
                                                        s += p[k][i] * p[k][j];
                                                f = s / h;
                                                for (k = i; k < m; k++)
                                                        p[k][j] += f * p[k][i];
                                        }
                                }
                                for (k = i; k < m; k++)
                                        p[k][i] *= scale;
                        }
                }
                d[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m && i != nm1)
                {
                        for (k = l; k < n; k++)
                                scale += fabs(p[i][k]);
                        if (scale)
                        {
                                for (k = l; k < n; k++)
                                {
                                        p[i][k] /= scale;
                                        s += p[i][k] * p[i][k];
                                }
                                f = p[i][l];
                                g = -SIGN(sqrt(s), f);
                                h = f * g - s;
                                p[i][l] = f - g;
                                for (k = l; k < n; k++)
                                        r[k] = p[i][k] / h;
                                if (i != mm1)
                                {
                                        for (j = l; j < m; j++)
                                        {
                                                for (s = 0.0, k = l; k < n; k++)
                                                        s += p[j][k] * p[i][k];
                                                for (k = l; k < n; k++)
                                                        p[j][k] += s * r[k];
                                        }
                                }
                                for (k = l; k < n; k++)
                                        p[i][k] *= scale;
                        }
                }
                anorm = MAX(anorm, fabs(d[i]) + fabs(r[i]));
        }

        /* Accumulation of right-hand transformations */
        for (i = n - 1; i >= 0; i--)
        {
                if (i < nm1)
                {
                        if (g)
                        {
                                for (j = l; j < n; j++)
                                        q[j][i] = (p[i][j] / p[i][l]) / g;
                                for (j = l; j < n; j++)
                                {
                                        for (s = 0.0, k = l; k < n; k++)
                                                s += p[i][k] * q[k][j];
                                        for (k = l; k < n; k++)
                                                q[k][j] += s * q[k][i];
                                }
                        }
                        for (j = l; j < n; j++)
                                q[i][j] = q[j][i] = 0.0;
                }
                q[i][i] = 1.0;
                g = r[i];
                l = i;
        }
        /* Accumulation of left-hand transformations */
        for (i = n - 1; i >= 0; i--)
        {
                l = i + 1;
                g = d[i];
                if (i < nm1)
                        for (j = l; j < n; j++)
                                p[i][j] = 0.0;
                if (g)
                {
                        g = 1.0 / g;
                        if (i != nm1)
                        {
                                for (j = l; j < n; j++)
                                {
                                        for (s = 0.0, k = l; k < m; k++)
                                                s += p[k][i] * p[k][j];
                                        f = (s / p[i][i]) * g;
                                        for (k = i; k < m; k++)
                                                p[k][j] += f * p[k][i];
                                }
                        }
                        for (j = i; j < m; j++)
                                p[j][i] *= g;
                } else
                        for (j = i; j < m; j++)
                                p[j][i] = 0.0;
                ++p[i][i];
        }
        /* diagonalization of the bidigonal form */
        for (k = n - 1; k >= 0; k--)
        {                       /* loop over singlar values */
                for (its = 0; its < 30; its++)
                {               /* loop over allowed iterations */
                        flag = 1;
                        for (l = k; l >= 0; l--)
                        {       /* test for splitting */
                                nm = l - 1;     /* note that r[l] is always
                                                 * zero */
                                if (fabs(r[l]) + anorm == anorm)
                                {
                                        flag = 0;
                                        break;
                                }
                                if (fabs(d[nm]) + anorm == anorm)
                                        break;
                        }
                        if (flag)
                        {
                                c = 0.0;        /* cancellation of r[l], if
                                                 * l>1 */
                                s = 1.0;
                                for (i = l; i <= k; i++)
                                {
                                        f = s * r[i];
                                        if (fabs(f) + anorm != anorm)
                                        {
                                                g = d[i];
                                                h = radius(f, g);
                                                d[i] = h;
                                                h = 1.0 / h;
                                                c = g * h;
                                                s = (-f * h);
                                                for (j = 0; j < m; j++)
                                                {
                                                        y = p[j][nm];
                                                        z = p[j][i];
                                                        p[j][nm] = y * c + z * s;
                                                        p[j][i] = z * c - y * s;
                                                }
                                        }
                                }
                        }
                        z = d[k];
                        if (l == k)
                        {       /* convergence */
                                if (z < 0.0)
                                {
                                        d[k] = -z;
                                        for (j = 0; j < n; j++)
                                                q[j][k] = (-q[j][k]);
                                }
                                break;
                        }
                        if (its == 30)
                        {
                                //error("svd: No convergence in 30 svd iterations", non_fatal);
                                return;
                        }
                        x = d[l];       /* shift from bottom 2-by-2 minor */
                        nm = k - 1;
                        y = d[nm];
                        g = r[nm];
                        h = r[k];
                        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                        g = radius(f, 1.0);
                        /* next QR transformation */
                        f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
                        c = s = 1.0;
                        for (j = l; j <= nm; j++)
                        {
                                i = j + 1;
                                g = r[i];
                                y = d[i];
                                h = s * g;
                                g = c * g;
                                z = radius(f, h);
                                r[j] = z;
                                c = f / z;
                                s = h / z;
                                f = x * c + g * s;
                                g = g * c - x * s;
                                h = y * s;
                                y = y * c;
                                for (jj = 0; jj < n; jj++)
                                {
                                        x = q[jj][j];
                                        z = q[jj][i];
                                        q[jj][j] = x * c + z * s;
                                        q[jj][i] = z * c - x * s;
                                }
                                z = radius(f, h);
                                d[j] = z;       /* rotation can be arbitrary
                                                 * id z=0 */
                                if (z)
                                {
                                        z = 1.0 / z;
                                        c = f * z;
                                        s = h * z;
                                }
                                f = (c * g) + (s * y);
                                x = (c * y) - (s * g);
                                for (jj = 0; jj < m; jj++)
                                {
                                        y = p[jj][j];
                                        z = p[jj][i];
                                        p[jj][j] = y * c + z * s;
                                        p[jj][i] = z * c - y * s;
                                }
                        }
                        r[l] = 0.0;
                        r[k] = f;
                        d[k] = x;
                }
        }
        free(r);

                // dhli add: the original code does not sort the eigen value
                // should do that and change the eigen vector accordingly

}



