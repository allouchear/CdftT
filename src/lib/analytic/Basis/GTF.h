#ifndef CDFTT_GTF_H_INCLUDED
#define CDFTT_GTF_H_INCLUDED

#include<iostream>
#include<iomanip>
#include<vector>
#include<analytic/Utils/MathFunction.h>

using namespace std;

	//! A GTF class.
	/*! This class will be used in the CGTF class. */

class GTF
{
	private:
		double _exposant;
		double _coefficient;
		vector<double> _coord;
		vector<int> _l;
		Binomial _bino;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one GTF on 0 or "None" value. */

		GTF();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for one GTF. */

		GTF(const double&, const double&, const vector<double>&, const vector<int>&, Binomial&);

			//! A default desctructor.
			/*! We don't use it. */

		~GTF(){}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The exposant value. */

		double exposant()
		{
			return _exposant;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The coefficient value. */

		double coefficient()
		{
			return _coefficient;
		}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The coordinates of the atom's GTF. */

		vector<double> coord()
		{
			return _coord;
		}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The quantic number of the atom's GTF. */

		vector<int> l()
		{
			return _l;
		}

			//! A normal member taking no arguments and returning a binomial value.
			/*! \return The binomial object using for calculus. */

		Binomial& bino()
		{
			return _bino;
		}

			//! A normal member taking one argument and returning a double value.
			/*! \return norme value. */

		double normeGTF();

			//! A normal member taking one argument and returning a double value.
			/*! \return norme value. */

		double normeGTF(GTF& q);

			//! A normal member taking one argument.
			/*! This member normalise the radial's composant of the GTF. */

		void normalise_radialGTF();

			//! A normal member taking one argument.
			/*! This member normalise the GTF. */

		void normaliseGTF();

			//! A normal member taking one argument and returning a double value.
			/*! \return The GTF overlap. */

		double overlapGTF(GTF&);

			//! A normal member taking two arguments and returning a double value.
			/*! \return The GTF overlap. */

		double overlap3GTF(GTF&, GTF&);

			//! A normal member taking three arguments and returning a double value.
			/*! \return The GTF overlap. */

		double overlap4GTF(GTF&, GTF&, GTF&);

			//! A normal member taking one argument and returning a double value.
			/*! \return The result of an integral with two GTF. */

		double GTFstarGTF(GTF&);
		
			//! A normal member taking two arguments and returning a GTF value.
			/*! \return The result of an integral with three GTF. */
			
		double GTFstarGTFstarGTF(GTF&, GTF&);
		
			//! A normal member taking three arguments and returning a GTF value.
			/*! \return The result of an integral with four GTF. */
			
		double GTFstarGTFstarGTFstarGTF(GTF&, GTF&, GTF&);
		
			//! A normal member taking four arguments and returning a double value.
			/*! \return The result of an integral with three ???. */
			
		double GTFxyzGTF(GTF&, int, int, int);
		
			//! A normal member taking one argument and returning a double value.
			/*! \return The kinetic ???. */
			
		double kineticGTF(GTF&);
		
			//! A normal member taking three arguments and returning a double value.
			/*! \return The ionic potential ???. */
			
		double ionicPotentialGTF(GTF&, vector<double>, double);
		
			//! A normal member taking three arguments and returning a double value.
			/*! \return The eri ???. */
			
		double ERIGTF(GTF&, GTF&, GTF&);
			
		double func(double x, double y, double z);
		
			//! An operator member taking one argument and returning a void value.
			
		void operator*=(double);
		
			//! An operator member taking one argument and returning a void value.
			
		void operator/=(double);

			//! A normal membre taking five arguments and returning a void value.
			/*! Insert all the data in the GTF. */
		void push_back(const double&, const double&, const vector<double>&, const vector<int>&, Binomial&);

};

	//! An operator member taking two arguments and returning a bool value.
	/*! \return The bool value of an equality between two CGTF. */
bool operator==(GTF, GTF);

	//! An operator member taking two arguments and returning an ostream value.
	/*! Print all the data of one GTF */
ostream& operator<<(ostream&, GTF&);
double operator*(vector<GTF>, vector<double>);

#endif
