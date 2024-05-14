#ifndef CDFTT_GTF_H_INCLUDED
#define CDFTT_GTF_H_INCLUDED

#include<iostream>
#include<cmath>
#include"MathFunction.h"

using namespace std;

	//! A GTF class.
	/*! This class will be used in the CGTF class. */

class GTF
{
	private:
		double _exposant;
		double _coefficient;
		vector<double> _coord;
		vector<double> _l;
	public:

			//! A constructor.
			/*! This constructor is used to add all of the parameters for one GTF. */

		GTF(const double&, const double&, const vector<double>&, const vector<double>&);

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

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The quantic number of the atom's GTF. */

		vector<double> l()
		{
			return _l;
		}

			//! A normal member taking one argument and returning a double value.
			/*! \return norme value. */

		double normeGTF(const Factorial&);

			//! A normal member taking one argument.
			/*! This member normalise the radial's composant of the GTF. */

		void normalise_radialGTF(const Factorial&);

			//! A normal member taking one argument.
			/*! This member normalise the GTF. */

		void normaliseGTF(const Factorial&);

			//! A normal member taking two arguments and returning a GTF value.
			/*! \return The GTF overlap. */

		GTF overlapGTF(const GTF&, const GTF&);

			//! An operator member taking two arguments and returning a vector<double>.
			/*! \return The sum of two vector. */

		operator+(const vector<double>&, const vector<double>&)

			//! An operator member taking two arguments and returning a GTF value.
			/*! \return The product of two GTF. */

		operator*(const GTF&, const GTF&);
};

#endif