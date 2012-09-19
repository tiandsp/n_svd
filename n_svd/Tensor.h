#ifndef _TENSOR_H
#define _TENSOR_H
#include "Mat.h"
#include "N_SVD.h"

class Mat;
class N_SVD;
class Tensor
{
public:
	Tensor(long a,...);
	Tensor(const Tensor &T);
	~Tensor();

	Tensor &operator=(const Tensor &);
	bool operator==(const Tensor &) const;
	bool operator!=(const Tensor &) const;

	Tensor operator+(Tensor &);
	Tensor operator+(double);

	Tensor operator-(Tensor &);
	Tensor operator-(double);

	Tensor operator*(Tensor &);	//点乘
	Tensor operator*(double);

	Tensor operator/(Tensor &); //点除
	Tensor operator/(double);

	N_SVD *m_mode_svd();


	Mat getSize();
	long getDim();

	void setElement(double,long,...);
	double getElement(long,...);

private:
	Mat tfastsvd();

private:
	long *N;	//每个维数的值
	long dim;	//数据的维数
	double *data;	//每个数据的值
	long n;	//数据的个数
	long *yuji;
};


#endif