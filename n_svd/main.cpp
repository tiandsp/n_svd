#include "Mat.h"
#include "BaseData.h"
#include <iostream>
using namespace std;
int main()
{
	double in[] = {
		 12, 51,   4 ,
		  6, 167, 68 ,
		 4,  24, 41 ,
	};

	int r=3;
	int c=3;

	Mat m(r,c);
	Mat s;
	Mat n;
	m.MatCopy(in,r,c);
	s=m.sqrtM();
//	s.eye(3);
//	n=m.dotMultiplication(s);

	cout<<m<<endl;
	cout<<s<<endl;
//	cout<<n<<endl;
//	cout<<m.det()<<endl;


//s=m.wshift(4);
//	s.eye(5);
//	s=m;
//	s.print();

//	string str="e://data.txt";
//	s.ReadFileData(str);
//	s.print();
	

	system("pause");
	return 0;
}