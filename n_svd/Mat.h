#ifndef _MAT_H
#define _MAT_H
#include "BaseData.h"

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
	Mat dot();
	Mat sqrtM();
	Mat mean(int n);


	void setElement(double,long,long);
	double getElement(long,long);
	void print();

private:
	void svd(int m, int n, double **a, double **p, double *d, double **q);

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