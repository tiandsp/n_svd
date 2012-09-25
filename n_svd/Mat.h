#ifndef _MAT_H
#define _MAT_H
#include "BaseData.h"
#include "Tensor.h"
#include <iostream>
using namespace std;

class Tensor;
class SVD;
class Mat
{
public:
	Mat(long r=0,long c=0);
	Mat(const Mat &M);
	~Mat();

	long getRow() const;
	long getCol() const;
	long getSize() const;

	Mat &operator=(const Mat &);
	bool operator==(const Mat &) const;
	bool operator!=(const Mat &) const;
	friend ostream& operator<<(ostream&,const Mat&);

	Mat operator~();		//����ת��
	Mat operator+(Mat &);
	Mat operator+(double);

	Mat operator-(Mat &);
	Mat operator-(double);

	Mat operator*(Mat &);	//������
	Mat operator*(double);

	Mat operator/(Mat &);	//������
	Mat operator/(double);

	Mat dotMultiplication(Mat &);	//���
	Mat dotDivision(Mat &);		//���
	Mat inv();	//��������
	void svd(SVD*);
	Mat tfastsvd();

	Mat dot();
	Mat sqrtM();
	Mat mean(int n);
	Mat reshape(long,long);
	Tensor *reshape(Mat);
	long prod();
	long prod(long,long);

	void setElement(double,long,long);
	double getElement(long,long);
	void print();
	Mat repmat(long,long);
	Tensor *tensorize(long,Mat);

	Mat wshift(long);


private:
	void svd(int m, int n, double **a, double **p, double *d, double **q);
	void reverse(Mat&,long,long);

private:
	double **data;
	long row;
	long col;
	long size;

	double **u;
	double *s;
	double **v;

};

#endif