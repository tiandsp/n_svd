#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{


	int r=3;
	int c=3;

	Mat m(r,c);
	
	m.setElement(3,0,0);
	m.setElement(2,0,1);
	m.setElement(-2,0,2);
	
	m.setElement(1,1,0);
	m.setElement(1,1,1);
	m.setElement(-3,1,2);

	m.setElement(-2,2,0);
	m.setElement(-3,2,1);
	m.setElement(11,2,2);


	cout<<m<<endl;
//	cout<<m.det()<<endl;
	Mat s;

//s=m.wshift(4);
	s.eye(5);
	s=m;
	s.print();

	system("pause");
	return 0;
}