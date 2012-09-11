#ifndef _TENSOR_H
#define  _TENSOR_H


class Tensor
{
public:
	Tensor();
	Tensor(long a,...);
	Tensor(double **data,long a,...);

	~Tensor();



private:
	long *N;
	long n;
	double element;
	double **data;

};


#endif