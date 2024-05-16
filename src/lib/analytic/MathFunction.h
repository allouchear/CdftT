#ifndef CDFTT_MATHFUNCTION_H_INCLUDED
#define CDFTT_MATHFUNCTION_H_INCLUDED

#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

	//! A factorial class.
	/*! This class will be used in the GTF class to for the normalisation. */

class Factorial
{
	private:
		vector<double> _tab;	
	public:
		Factorial();

			//! A real constructor.
			/*! This constructor is used to create a table from 0 to n factorial. */

		Factorial(int);

			//! A default desctructor.
			/*! We don't use it. */

		~Factorial(){}

			//! A normal member taking one argument and returning a double value.
			/*! \return The n factorial value. */

		double factorial(int);

			//! A normal member taking one argument and returning a double value.
			/*! \return The n double factorial value. */

		double double_factorial(int);
};

class Binomial
{
	private:
		vector<vector<double>> _tab;
		Factorial _fact;
	public:
		Binomial(int, int, Factorial&);
		Binomial();
		~Binomial(){}
		double binomial(int, int);
		vector<vector<double>>& tab()
		{
			return _tab;
		}
		Factorial& fact()
		{
			return _fact;
		}
};

double power(double, double);
double f(int, int, int, double, double, Binomial&);
double Theta(int, int, int, int, double, double, double, Binomial&);
int m1p(int);
double A(int, int, int, int, int, double, double, double, double, Binomial&);
double B(int,int,int , int , int , double , double , double , double , Binomial&);
double myGamma(int, Factorial&);
double F(int,double, Factorial&);
vector<double> getFTable(int, double);

#endif