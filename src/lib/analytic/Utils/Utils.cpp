#include<iostream>
#include<vector>
#include <analytic/Utils/Utils.h>

using namespace std;

Factorial::Factorial()
{
	_tab = vector<double>();
}

Factorial::Factorial(int n)
{
	_tab = vector<double>(n,1);
	for(int i=2; i<n; i++)
		for(int k=0; k<=i/2-1; k++)
			_tab[i]*= i-2*k;
}

double Factorial::factorial(int n)
{
	if(n==0)
		return 1;

	return double_factorial(n)*double_factorial(n-1);
}

double Factorial::double_factorial(int n)
{
	if(n<0)
		return 1;
	
	if(size_t(n)>=_tab.size())
	{
	double r;
		for(size_t i=_tab.size(); i<size_t(n); i++)
		{	
			r=1;
			for(size_t k=0; k<=i/2-1; k++)
				r*= i-2*k;	
			_tab.push_back(r);
		}
	}

	return _tab[n];
}

Binomial::Binomial()
{
	_fact=Factorial();
	_tab=vector<vector<double>>();
}

Binomial::Binomial(int i, Factorial& F) : _fact(F)
{
	vector<double>V(0);
	_tab=vector< vector<double> >(i,V);
	for(size_t k=0; k<_tab.size(); k++)
	{
		_tab[k].resize(k+1);
		for(size_t l=0; l<=k; l++)
			_tab[k][l] = _fact.factorial(k)/_fact.factorial(l)/_fact.factorial(k-l);
	}
}

double Binomial::binomial(int i, int j)
{
/*													A debug
	if(i>_tab.size() || j>_tab[0].size())
	{
		cout<<"Test 1"<<endl;
		_tab=Binomial(i, j ,F).tab();
		cout<<"Test 2"<<endl;
	}
*/
	return _tab[i][j];
}

double power(double e, int n)
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
	for(k=1; k<=n; k++)
		p*=e;

	return p;
}

double f(int i, int l, int m, double A, double B, Binomial& Bi)
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

	return sum;
}

double Theta(int i,int r,int l1,int l2, double A, double B, double g, Binomial& Bi)
{
	return f(i,l1,l2,A,B,Bi)*Bi.fact().factorial(i)/Bi.fact().factorial(r)/Bi.fact().factorial(i-2*r)/pow(g,i-r);
}

int m1p(int i)
{
	if(i%2==0) return 1;
	else return -1;
}

double A(int i,int r, int u,int l1,int l2, double A, double B, double C,double g, Binomial& Bi)
{
	return m1p(i+u)*f(i,l1,l2,A,B,Bi)*Bi.fact().factorial(i)*power(C,i-2*(r+u))/
		(Bi.fact().factorial(r)*Bi.fact().factorial(u)*Bi.fact().factorial(i-2*r-2*u)*power(4*g,r+u));
}

double B(int i,int ip,int r, int rp, int u, double PQ, double d, double T1, double T2, Factorial& Fa)
{
	int ii=i+ip-2*r-2*rp;
	return m1p(ip+u)*T1*T2*Fa.factorial(ii)/Fa.factorial(u)/Fa.factorial(ii-2*u)*pow(PQ,ii-2*u)/(pow(4.0,i+ip-r-rp)*pow(d,ii-u));
}

double myGamma(int n, Factorial& Fa)
{
	return Fa.double_factorial(2*n-1)*sqrt(M_PI)/power(2,n);
}

double F(int n,double t, Factorial& Fa)
{
	double et=exp(-t);
	double twot=2*t;
	double T=0.0;
	double x=1.0;
	int i=0;
    double DD=1.0;
	double TMAX = 50.0;
	int MAXFACT = 200;
	double acc = 1e-16;

	if(fabs(t)<=acc) 
		return 1/(double)(2*n+1);

	if(t>=TMAX)
		return myGamma(n, Fa)/power(t,n)/2/sqrt(t);


	while(fabs(x/T)>acc && (n+i)<MAXFACT)
	{	
        x=Fa.double_factorial(2*n-1)/Fa.double_factorial(2*(n+i+1)-1)*DD;
		T += x;
		i++;
        DD *= twot;
	}
	if(n+i>=MAXFACT)
	{
		cout<<"Divergence in F, Ionic integrals"<<endl;
		exit(1);
	}
	T *=et;

	return T;
}

vector<double> getFTable(int mMax, double t, Factorial& Fa)
{
	double tCritic = 30.0;
	vector<double> Fmt(mMax+1);
	int m;
	if(t>tCritic)
	{
		Fmt[0] = sqrt(M_PI/t) * 0.5;
		for(m=1; m<=mMax; m++)
			Fmt[m] = Fmt[m-1] * (m-0.5) / t;
		return Fmt;
	}
	Fmt[mMax] = F(mMax,t, Fa);
	double expt = exp(-t);
	double twot = 2*t;
	for(m = mMax-1; m>=0; m--)
		Fmt[m] = (twot * Fmt[m+1] + expt) / (m*2+1);
	return Fmt;
}