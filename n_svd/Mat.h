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

	long getRow() const;	//返回矩阵行数
	long getCol() const;	//返回矩阵列数
	long getSize() const;	//返回矩阵数据总个数

	Mat &operator=(const Mat &);		
	bool operator==(const Mat &) const;
	bool operator!=(const Mat &) const;
	friend ostream& operator<<(ostream&,const Mat&);

	Mat operator~();		//矩阵转置
	Mat operator+(Mat &);	//两矩阵相加
	Mat operator+(double);	//矩阵加一个常数

	Mat operator-(Mat &);	//两矩阵相减
	Mat operator-(double);	//矩阵减一个常数

	Mat operator*(Mat &);	//矩阵乘一般乘法
	Mat operator*(double);	//矩阵乘以一个常数

	Mat operator/(Mat &);	//矩阵叉除，矩阵乘另一个矩阵的逆
	Mat operator/(double);	//矩阵除以一个常数

	void eye(long);		//产生n*n的单位矩阵
	Mat dotMultiplication(Mat &);	//点积，矩阵元素分别对应乘
	Mat dotDivision(Mat &);			//点除，矩阵元素分别对应除
	Mat inv();			//矩阵求逆
	void svd(SVD*);		
	Mat tfastsvd();

	Mat dot();
	Mat sqrtM();		//矩阵每个元素分别求根
	Mat mean(int n);	//n=1或2，为1时求每一列的平均值，返回为一行原列的矩阵；为2时求每一行的平均值，返回为一列原行的矩阵
	Mat reshape(long,long);		//改变矩阵的行和列，必须满足“原行*原列=新行*新列”
	Tensor *reshape(Mat);		//将矩阵改变为以Mat为索引的张量，这里Mat为一行n列，张量为n维，每一维的数据的个数为Mat对应的数值
	long prod();		//返回矩阵中所有元素的乘积
	long prod(long,long);	//返回除了（long,long）这个位置之外的矩阵所有元素的乘积

	double norm();		//返回矩阵所有元素的平方和的根
	double det();		//返回矩阵的行列式

	void setElement(double,long,long);
	double getElement(long,long);
	void MatCopy(double*,long,long);	//从一个数组中拷贝数据到Mat中
	void ReadFileData(string);			//从文件中读取数据到Mat中

	void print();				//输出矩阵的元素
	Mat repmat(long,long);		//把矩阵元素周期扩展
	Tensor *tensorize(long,Mat);

	Mat wshift(long);		//把Mat循环左移n位，暂时只支持一行的矩阵


private:
	void create(long,long);		//真正的初始化
	void clear();				//清空Mat所占的内存
	void svd(int m, int n, double **a, double **p, double *d, double **q);
	void reverse(Mat&,long,long);		//供wshift()使用，不提供外部调用

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