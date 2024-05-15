#include<iostream>
#include<vector>
#include <analytic/MathFunction.h>

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

Binomial::Binomial(int i, int j, Factorial& F)
{
	_tab.resize(i, vector<double>(j));
	for(i=0; i<_tab.size(); i++)
		for(j=0; j<=i; j++)
			_tab[i][j] = F.factorial(i)/F.factorial(j)/F.factorial(i-j);
}

double Binomial::binomial(int i, int j, Factorial& F)
{
/*													A debbug
	if(i>_tab.size() || j>_tab[0].size())
	{
		cout<<"Test 1"<<endl;
		_tab=Binomial(i, j ,F).tab();
		cout<<"Test 2"<<endl;
	}
*/
	return _tab[i][j];
}

double power(double e, double n)
{
	double p=1.0;
	int k;
	if(fabs(e)<1e-10)
	{
		if(n==0)
			return 1.0;
		else
			return 0.0;
	}
	for(k=1; k<n; k++)
		p*=e;
	return p;
}

double f(int i, int l, int m, double A, double B, const Binomial& Bi)
{
	int j, jmin, jmax;
	double sum=0.0;

	jmin=0;
	if(jmin<i-m)
		jmin=i-m;
	jmax=i;
	if(jmax>l)
		jmax=l;
	for(j=jmin; j<=jmax; j++)
		sum+=Bi.binomial(l,j)*Bi.binomial(m, i-j)*power(-A, l-j)*power(-B, m-i+j);
}

double Theta(int i,int r,int l1,int l2, double A, double B, double g, const Factorial& Fa)
{
	return f(i,l1,l2,A,B)*Fa.factorial(i)/Fa.factorial(r)/Fa.factorial(i-2*r)/pow(g,i-r);
}

int m1p(int i)
{
	if(i%2==0) return 1;
	else return -1;
}

double A(int i,int r, int u,int l1,int l2, double A, double B, double C,double g, const Factorial& Fa)
{
	return m1p(i+u)*f(i,l1,l2,A,B)*Fa.factorial(i)*power(C,i-2*(r+u))/
		(Fa.factorial(r)*Fa.factorial(u)*Fa.factorial(i-2*r-2*u)*power(4*g,r+u));
}

double B(int i,int ip,int r, int rp, int u, double PQ, double d, double T1, double T2, const Factorial& Fa)
{
	int ii=i+ip-2*r-2*rp;
	return m1p(ip+u)*T1*T2*Fa.factorial(ii)/Fa.factorial(u)/Fa.factorial(ii-2*u)*pow(PQ,ii-2*u)/(pow(4.0,i+ip-r-rp)*pow(d,ii-u));
}
