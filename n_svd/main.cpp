#include "Mat.h"
#include <iostream>
using namespace std;
int main()
{
	Mat m(3,2);
	m.setElement(1,0,0);
	m.setElement(4,0,1);
	m.setElement(2,1,0);
	m.setElement(5,1,1);
	m.setElement(3,2,0);
	m.setElement(6,2,1);

	m=~m;
	m.print();
	m.svd();

	for (int i=0;i<3;i++)
	{
		for (int j=0;j<3;j++)
		{
			cout<<m.Svd.u[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;

	for (int i=0;i<2;i++)
	{
		cout<<m.Svd.s[i]<<"  ";
	}

	cout<<endl<<endl;

	for (int i=0;i<2;i++)
	{
		for (int j=0;j<2;j++)
		{
			cout<<m.Svd.u[i][j]<<"  ";
		}
		cout<<endl;
	}
	cout<<endl;


	system("pause");
	return 0;
}