#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{


	int r=2;
	int c=2;
//	SVD *re;

//	re=new SVD;
//	re->u=new Mat(r,r);
//	re->s=new Mat(1,c);
//	re->v=new Mat(c,c);
	Mat m(r,c);
	m.setElement(1,0,0);
	m.setElement(4,0,1);
//	m.setElement(5,0,2);
	m.setElement(3,1,0);
	m.setElement(2,1,1);

	m.print();
//	m=~m;
//	m.svd(re);

//	delete re->u;
//	delete re->s;
//	delete re->v;

//	delete re;

	//m=m.dot();
	m=m.reshape(4,1);
	m.print();

	system("pause");
	return 0;
}