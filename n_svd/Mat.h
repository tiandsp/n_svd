#ifndef _MAT_H
#define _MAT_H

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
	Mat operator-(Mat &);
	Mat operator*(Mat &);	//������
	Mat operator/(Mat &);	//������

	Mat dotMultiplication(Mat &);
	Mat dotDivision(Mat &);
	Mat inv();	//��������

	void setElement(double,long,long);
	double getElement(long,long);
	void print();

private:
	double **data;
	long row;
	long col;
	long size;

};

#endif