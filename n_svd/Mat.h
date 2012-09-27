#ifndef _MAT_H
#define _MAT_H
#include "BaseData.h"
#include "Tensor.h"
#include <iostream>
#include <string>
using namespace std;

class Tensor;
class SVD;
class Mat
{
public:
	Mat(long r=0,long c=0);
	Mat(const Mat &M);
	~Mat();

	long getRow() const;	//���ؾ�������
	long getCol() const;	//���ؾ�������
	long getSize() const;	//���ؾ��������ܸ���

	Mat &operator=(const Mat &);		
	bool operator==(const Mat &) const;
	bool operator!=(const Mat &) const;
	friend ostream& operator<<(ostream&,const Mat&);

	Mat operator~();		//����ת��
	Mat operator+(Mat &);	//���������
	Mat operator+(double);	//�����һ������

	Mat operator-(Mat &);	//���������
	Mat operator-(double);	//�����һ������

	Mat operator*(Mat &);	//�����һ��˷�
	Mat operator*(double);	//�������һ������

	Mat operator/(Mat &);	//���������������һ���������
	Mat operator/(double);	//�������һ������

	void eye(long);		//����n*n�ĵ�λ����
	Mat dotMultiplication(Mat &);	//���������Ԫ�طֱ��Ӧ��
	Mat dotDivision(Mat &);			//���������Ԫ�طֱ��Ӧ��
	Mat inv();			//��������
	void svd(SVD*);		
	Mat tfastsvd();

	Mat dot();
	Mat sqrtM();		//����ÿ��Ԫ�طֱ����
	Mat mean(int n);	//n=1��2��Ϊ1ʱ��ÿһ�е�ƽ��ֵ������Ϊһ��ԭ�еľ���Ϊ2ʱ��ÿһ�е�ƽ��ֵ������Ϊһ��ԭ�еľ���
	Mat reshape(long,long);		//�ı������к��У��������㡰ԭ��*ԭ��=����*���С�
	Tensor *reshape(Mat);		//������ı�Ϊ��MatΪ����������������MatΪһ��n�У�����Ϊnά��ÿһά�����ݵĸ���ΪMat��Ӧ����ֵ
	long prod();		//���ؾ���������Ԫ�صĳ˻�
	long prod(long,long);	//���س��ˣ�long,long�����λ��֮��ľ�������Ԫ�صĳ˻�

	double norm();		//���ؾ�������Ԫ�ص�ƽ���͵ĸ�
	double det();		//���ؾ��������ʽ

	void setElement(double,long,long);
	double getElement(long,long);
	void MatCopy(double*,long,long);	//��һ�������п������ݵ�Mat��
	void ReadFileData(string);			//���ļ��ж�ȡ���ݵ�Mat��

	void print();				//��������Ԫ��
	Mat repmat(long,long);		//�Ѿ���Ԫ��������չ
	Tensor *tensorize(long,Mat);

	Mat wshift(long);		//��Matѭ������nλ����ʱֻ֧��һ�еľ���


private:
	void create(long,long);		//�����ĳ�ʼ��
	void clear();				//���Mat��ռ���ڴ�
	void svd(int m, int n, double **a, double **p, double *d, double **q);
	void reverse(Mat&,long,long);		//��wshift()ʹ�ã����ṩ�ⲿ����

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