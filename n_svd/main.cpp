#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{


	int r=1;
	int c=2;

	Mat m(r,c);
	m.setElement(1,0,0);
	m.setElement(4,0,1);

//	m.print();
	cout<<m<<endl;
	Mat s;

	s=m.repmat(3,3);


	s.print();

	system("pause");
	return 0;
}