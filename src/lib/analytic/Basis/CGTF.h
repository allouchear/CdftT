#ifndef CDFTT_CGTF_H_INCLUDED
#define CDFTT_CGTF_H_INCLUDED

#include<iostream>
#include<analytic/Basis/GTF.h>

using namespace std;

	//! A CGTF class.
	/*! This class will be used in the LCAO class. */

class CGTF
{
	private:
		int _num_center;
		int _numberOfFunctions;
		vector<double> _coefficients;
		vector<GTF> _gtf;
		int _L;
		int _M;
		Binomial _bino;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one CGTF on 0 or "None" value. */

		CGTF();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for one CGTF. */

		CGTF(vector<GTF>);

			//! A default desctructor.
			/*! We don't use it. */

		~CGTF(){}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of GTF. */

		void setCoef(double);

		vector<double> coefficients() {return _coefficients;}

		int numberOfFunctions()
		{
			return _numberOfFunctions;
		}

			//! A normal member taking no arguments and returning a vector<GTF> value.
			/*! \return The table of GTF which compose de CGTF value. */

		vector<GTF> gtf()
		{
			return _gtf;
		}

			//! A normal member taking no arguments and returning a Binomial value.
			/*! \return The binomial value. */

		Binomial bino()
		{
			return _bino;
		}

			//! A normal member taking three arguments and returning a double value.
			/*! \return The ERI value   ???? */

		double ERICGTF(CGTF&, CGTF&, CGTF&);

			//! A normal member taking no arguments and returning a void value.
			/*! Normalise the CGTF */

		void normaliseCGTF();

			//! A normal member taking one argument and returning a double value.
			/*! \return The overlap value between two CGTF. */

		double overlapCGTF(CGTF&);
		
			//! A normal member taking two arguments and returning a double value.
			/*! \return The overlap value between three CGTF. */

		double overlap3CGTF(CGTF&, CGTF&);

			//! A normal member taking three arguments and returning a double value.
			/*! \return The overlap value between four CGTF. */

		double overlap4CGTF(CGTF&, CGTF&, CGTF&);

			//! A normal member taking four arguments and returning a double value.
			/*! \return ??? */

		double CGTFxyzCGTF(CGTF&, int, int, int);

			//! A normal member taking one argument and returning a double value.
			/*! \return The kinetic value ??? */

		double kineticCGTF(CGTF&);

			//! A normal member taking two arguments and returning a double value.
			/*! \return The ionic potential value ??? */

		double ionicPotentialCGTF(CGTF&, vector<double>, double);
		
			//! A normal member taking one argument and returning a double value.
			/*! \return ???. */

		double CGTFstarCGTF(CGTF&);

		//bool CGTFEqCGTF(CGTF&);

			//! A normal membre taking one argument and returning a void value.
			/*! Insert all the data in the CGTF. */
		void push_back(GTF&);

		void setNumCenter(int);

		double func(double x, double y, double z) const;
};

	//! An operator member taking two arguments and returning a bool value.
	/*! \return The bool value of an equality between two CGTF. */
bool operator==(CGTF, CGTF);

	//! An operator member taking two arguments and returning an ostream value.
	/*! Print all the data of one CGTF */
ostream& operator<<(ostream&, CGTF&);

double operator*(const vector<CGTF>&, const vector<double>&);

#endif
