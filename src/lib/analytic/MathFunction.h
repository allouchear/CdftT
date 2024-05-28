#ifndef CDFTT_MATHFUNCTION_H_INCLUDED
#define CDFTT_MATHFUNCTION_H_INCLUDED

#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

	//! A factorial class.
	/*! This class will be used in the GTF class for the normalisation. */

class Factorial
{
	private:
		vector<double> _tab;	
	public:

			//! A default constructor.
			/*! This constructor create a table without a size (0). */

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

	//! A binomial class.
	/*! This class will be used in the GTF class for some calculus. */

class Binomial
{
	private:
		vector<vector<double>> _tab;
		Factorial _fact;
	public:

		//! A real constructor.
		/*! This constructor is used to create a table from 0 to n binomial. */

		Binomial(int, Factorial&);

		//! A default constructor.
		/*! This constructor create a table without a size (0). */

		Binomial();

		//! A default desctructor.
		/*! We don't use it. */

		~Binomial(){}

		//! A normal member taking two arguments and returning a double value.
		/*! \return The i,j binomial value. */

		double binomial(int, int);

		//! A normal member taking no arguments and returning a vector<vector<double>> value.
		/*! \return The binomial table. */

		vector<vector<double>>& tab()
		{
			return _tab;
		}

		//! A normal member taking no arguments and returning a Factorial value.
		/*! \return The factorial table. */

		Factorial& fact()
		{
			return _fact;
		}
};

//Ne pas oublier de commenter
double power(double, int);
double f(int, int, int, double, double, Binomial&);
double Theta(int, int, int, int, double, double, double, Binomial&);
int m1p(int);
double A(int, int, int, int, int, double, double, double, double, Binomial&);
double B(int,int,int , int , int , double , double , double , double , Factorial&);
double myGamma(int, Factorial&);
double F(int,double, Factorial&);
vector<double> getFTable(int, double, Factorial&);

#endif