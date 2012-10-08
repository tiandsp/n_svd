#ifndef _BASEDATA_H
#define _BASEDATA_H

#include "Mat.h"
#include "Tensor.h"
class Mat;
class Tensor;

class SVD
{
public:
	SVD(){};
	~SVD(){};

	Mat *u;
	Mat *s;
	Mat *v;

};

class N_SVD
{
public:
	N_SVD(){};
	~N_SVD(){};


	Mat **U;
	Tensor *core;

};



#endif





