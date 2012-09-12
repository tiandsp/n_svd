#ifndef _MAT_H
#define _MAT_H

typedef struct _SVD
{
	double **u;
	double *s;
	double **v;

}SVD;

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

	Mat operator/(Mat &);	//¾ØÕó²æ³ı
	Mat operator/(double);

	void dotMultiplication(Mat &);
	void dotDivision(Mat &);
	void inv();	//¾ØÕóÇóÄæ
	void svd();

	void setElement(double,long,long);
	double getElement(long,long);
	void print();

public:
	SVD Svd;

private:
	void svd(int m, int n, double **a, double **p, double *d, double **q);

private:
	double **data;
	long row;
	long col;
	long size;

};

#endif