#ifndef _TENSOR_H
#define _TENSOR_H
#include "Mat.h"

class Mat;
class Tensor
{
public:
	Tensor();
	Tensor(long a,...);
	Tensor(const Tensor &T);
	~Tensor();

	Tensor &operator=(const Tensor &);
	bool operator==(const Tensor &) const;
	bool operator!=(const Tensor &) const;




	Mat getSize();
	long getDim();


private:
	long *N;
	long n;
	double element;
	double **data;

};


#endif