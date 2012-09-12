#ifndef _BASEDATA_H
#define _BASEDATA_H

#include "Mat.h"

class Mat;

class SVD
{
public:
	SVD(){};
	~SVD(){};

	Mat *u;
	Mat *s;
	Mat *v;

};





#endif





