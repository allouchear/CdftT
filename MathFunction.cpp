#include<iostream>
#include"MathFunction.h"

using namespace std;

Factorial::Factorial(int n)
{
	_tab.resize(n);
	_tab[0]=1;
	for(int i=1; i<=n; i++)
		_tab[i]=_tab[i-1]*i;
}

double Factorial::factorial(int n)
{
	if(n>_tab.size())
	{
		for(int i=_tab.size(); i<=n; i++)
			_tab.push_back(_tab[i-1]*i);
	}
	return _tab[n];
}

double Factorial::double_factorial(int n)
{
	return factorial(factorial(n));
}