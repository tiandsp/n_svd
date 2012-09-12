#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{
	int r=3;
	int c=2;
	SVD *re;


	Mat m(r,c);
	m.setElement(1,0,0);
	m.setElement(4,0,1);
	m.setElement(3,1,0);
	m.setElement(2,1,1);
	m.setElement(5,2,0);
	m.setElement(4,2,1);

	m=~m;
	m.print();
	re=m.svd();

	re->u->print();
	re->s->print();
	re->v->print();




/*
	for (int i=0;i<r;i++)
	{
		for (int j=0;j<r;j++)
		{
			cout<<m.Svd.u[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;

	for (int i=0;i<c;i++)
	{
		cout<<m.Svd.s[i]<<"  ";
	}

	cout<<endl<<endl;

	for (int i=0;i<c;i++)
	{
		for (int j=0;j<c;j++)
		{
			cout<<m.Svd.v[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;
*/


	system("pause");
	return 0;
}