#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{


	int r=3;
	int c=1;

	Mat m(r,c);
	m.setElement(1,0,0);
	m.setElement(4,1,0);
	m.setElement(3,2,0);

//	m.print();
	cout<<m;
	Mat s;

	s=m.repmat(1,2);


	s.print();

	system("pause");
	return 0;
}