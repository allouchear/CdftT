#include<iostream>
#include<vector>
#include <analytic/MathFunction.h>

using namespace std;

int main()
{
	Factorial Table(100);
	cout<<"test pre-initialisation"<<endl;
	Binomial Bin(100, 50, Table);
	cout<<"Test post-initialisation"<<endl;

	cout<<Bin.binomial(18,7, Table)<<endl;
	cout<<Table.factorial(18)/Table.factorial(7)/Table.factorial(11)<<endl;
	cout<<endl;
	cout<<Bin.binomial(20,18, Table)<<endl;
	cout<<Table.factorial(20)/Table.factorial(18)/Table.factorial(2);

	return 0;
}
