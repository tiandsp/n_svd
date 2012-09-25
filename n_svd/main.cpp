#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{


	int r=1;
	int c=5;

	Mat m(r,c);
	m.setElement(1,0,0);
	m.setElement(4,0,1);
	m.setElement(5,0,2);
	m.setElement(2,0,3);
	m.setElement(6,0,4);

//	m.print();
	cout<<m<<endl;
	Mat s;

//	s=m.repmat(3,3);
	s=m.wshift(4);

	s.print();

	system("pause");
	return 0;
}