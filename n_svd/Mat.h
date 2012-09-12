#ifndef _MAT_H
#define _MAT_H
#include "BaseData.h"

/*
typedef struct _SVD
{
	//	double **u;
	//	double *s;
	//	double **v;
	Mat u;
	Mat s;
	Mat v;

}SVD;
*/
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

	Mat operator~();		//¾ØÕó×ªÖÃ
	Mat operator+(Mat &);
	Mat operator+(double);

	Mat operator-(Mat &);
	Mat operator-(double);

	Mat operator*(Mat &);	//¾ØÕó²æ³Ë
	Mat operator*(double);

	Mat operator/(Mat &);	//¾ØÕó²æ³ý
	Mat operator/(double);

	void dotMultiplication(Mat &);	//µã»ý
	void dotDivision(Mat &);		//µã³ý
	void inv();	//¾ØÕóÇóÄæ
	SVD *svd();
//	Mat dot();
//	Mat sqrt();
//  Mat abs();

	void setElement(double,long,long);
	double getElement(long,long);
	void print();

public:
//	SVD Svd;

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