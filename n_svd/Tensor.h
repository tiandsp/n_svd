#ifndef _TENSOR_H
#define _TENSOR_H
#include "Mat.h"
#include "N_SVD.h"
#include <iostream>

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

	Tensor operator*(Tensor &);	//���
	Tensor operator*(double);

	Tensor operator/(Tensor &); //���
	Tensor operator/(double);

	N_SVD *m_mode_svd();


	Mat getSize();
	long getSize(long i);
	long getDim();

	void setElement(double,long,...);
	double getElement(long,...);

private:
	Mat matricize(Tensor,long);
	Tensor mode_m_prod(Mat,long);
	Tensor shiftdim(long);
	Mat reshape(long,long);

private:
	long *N;	//ÿ��ά����ֵ
	long dim;	//���ݵ�ά��
	double *data;	//ÿ�����ݵ�ֵ
	long n;	//���ݵĸ���
	long *yuji;
};


#endif