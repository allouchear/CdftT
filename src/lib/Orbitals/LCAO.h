#ifndef CDFTT_LCAO_H_INCLUDED
#define CDFTT_LCAO_H_INCLUDED

#include<iostream>
#include <Basis/CGTF.h>
#include <Utils/FCHK.h>

using namespace std;

	//! A LCAO class.
	/*! This class will be used in the Orbitals class. */

class LCAO
{
	private:
		vector<CGTF> _cgtf;
		int _numberOfFunctions;
		Binomial _bino;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one LCAO on 0 or "None" value. */

		LCAO();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for one LCAO. */

		LCAO(vector<CGTF>);

		LCAO(FCHK&, Binomial&);

			//! A default desctructor.
			/*! We don't use it. */

		~LCAO(){}

			//! A normal member taking no arguments and returning a vector<CGTF> value.
			/*! \return The table of CGTF which compose the LCAO. */

		vector<CGTF> cgtf() const &
		{
			return _cgtf;
		}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of CGTF in one LCAO. */

		int numberOfFunctions() const &
		{
			return _numberOfFunctions;
		}

			//! A normal member taking three arguments and returning a double value.
			/*! \return The ERI value ??? */

		double ERILCAO(LCAO&, LCAO&, LCAO&);

			//! A normal member taking no arguments and returning a void value.
			/*! Normalise the LCAO. */

		void normaliseLCAO();

			//! A normal member taking one argument and returning a double value.
			/*! \return The overlap value between two LCAO. */

		double overlapLCAO();

			//! A normal member taking two arguments and returning a double value.
			/*! \return The overlap value between three LCAO. */

		double overlap3LCAO(const vector<double>& c1, const vector<double>& c2, const vector<double>& c3);

			//! A normal member taking three arguments and returning a double value.
			/*! \return The overlap value between four LCAO. */

		double overlap4LCAO(const vector<double>& c1, const vector<double>& c2, const vector<double>& c3, const vector<double>& c4);

			//! A normal member taking one argument and returning a double value.
			/*! \return The kinetic value ??? */

		double kineticLCAO(LCAO&);

			//! A normal member taking two arguments and returning a double value.
			/*! \return The ionic potential value ??? */

		double ionicPotentialLCAO(LCAO&, vector<double>, double);

			//! A normal member taking one argument and returning a double value.
			/*! \return ???. */

		double LCAOstarLCAO(LCAO&);

			//! A normal member taking four arguments and returning a double value.
			/*! \return ???. */

		double LCAOxyzLCAO(LCAO&, int, int, int);

		//bool LCAOEqLCAO(LCAO&);

			//! A normal membre taking two arguments and returning a void value.
			/*! Insert all the data in the LCAO. */
		void push_back(const vector<CGTF>&);

		double func(const vector<double>& c, double x, double y, double z) const;

			//! An operator member taking two arguments and returning an ostream value.
			/*! Print all the data of one LCAO. */
		friend ostream& operator<<(ostream&, const LCAO&);
};

//! An operator member taking two arguments and returning a bool value.
/*! \return The bool value of an equality between two LCAO. */

bool operator==(LCAO, LCAO);

double operator*(const LCAO&, const vector<vector<double>>&);

#endif