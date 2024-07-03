#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include<iostream>
#include<analytic/Utils/WFX.h>
#include<analytic/Utils/FCHK.h>
#include<analytic/Utils/MOLDENGAB.h>
#include<analytic/Basis/CGTF.h>
#include<common/PeriodicTable.h>
#include<common/Descriptors.h>

using namespace std;


	//! An Orbitals class.
	/*! This class will be used to calculate descriptors. */
class Orbitals
{
	private:
		vector<CGTF> _vcgtf;
		vector<vector<vector<double>>> _coefficients;
		int _numberOfAo;
		int _numberOfMo;
		int _number_of_alpha_electrons;
		int _number_of_beta_electrons;
		int _number_of_atoms;
		vector<int> _primitive_centers;
		vector<int> _atomic_numbers;
		vector<string> _symbol;
		vector<vector<double>> _orbital_energy;
		vector<vector<double>> _all_f;
		vector<int> _numOrb;
		vector<vector<double>> _occupation_number;
		bool _alpha_and_beta;
		Binomial _bino;
		Descriptors _descriptors;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one LCAO on 0 or "None" value. */
		Orbitals();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for one LCAO. */
		Orbitals(WFX&, Binomial&, const PeriodicTable&);

		Orbitals(FCHK&, Binomial&, const PeriodicTable&);

		Orbitals(MOLDENGAB&, Binomial&, const PeriodicTable&);

		Orbitals(LOG&, Binomial&, const PeriodicTable&);

			//! A default desctructor.
			/*! We don't use it. */
		~Orbitals() {}

			//! A normal member taking no arguments and returning a vector<LCAO> value.
			/*! \return The table of LCAO which compose the Orbitals. */
		vector<CGTF>& vcgtf() {return _vcgtf;}

		vector<vector<vector<double>>>& coefficients() {return _coefficients;}

		vector<vector<double>>& Energy()  {return _orbital_energy;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of LCAO in the Orbitals. */
		int NumberOfAo() {return _numberOfAo;}

		int NumberOfMo() {return _numberOfMo;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of atoms. */
		int NumberOfAtoms() const {return _number_of_atoms;} 

			//! A normal member taking one argument and returning an int value.
			/*! \return The primitive center. */
		int PrimitiveCenter(int i) const {return _primitive_centers[i];}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The table of primitive centers. */
		vector<int> PrimitiveCenters() const {return _primitive_centers;}

			//! A normal member taking no arguments and returning a vector<string> value.
			/*! \return The table of symbol of atoms. */
		vector<string> symbol() const {return _symbol;}

		void NormaliseAllBasis();

		double ERIorbitals(Orbitals& q, Orbitals& r, Orbitals& s);

			//! A normal member taking two arguments and returning a double value.
			/*! \return The overlap between two LCAO which compose the Orbitals. */
		double Overlap(int i, int j, int alpha=0);

			//! A normal member taking two arguments and returning a void value.
			/*! Print the value of an overlap. */
		void PrintOverlap(int, int, int alpha=0);

		double Overlap3Orbitals(int i, int j, int k, int alpha=0);
		double Overlap4Orbitals(int i, int j, int k, int l, int alpha=0);
		double kinetic();
		double ionicPotential(vector<double> C, double Z);
		double OrbstarOrb();
		double OrbxyzOrb(int ix, int iy, int iz);
		void normaliseBasis();
		double func(double x, double y, double z) const;


			//! A normal member taking no arguments and returning an int value.
			/*! \return The HOMO Molecular Orbital number. */
		void HOMO();

			//! A normal member taking no arguments and returning an int value.
			/*! \return The LUMO Molecular Orbital number. */
		void LUMO();

			//! A normal member taking no arguments and returning a double value.
			/*! \return The HOMO's energy. */
		double eHOMO(int alpha=0) const {return _orbital_energy[alpha][_numOrb[0]];}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The LUMO's energy. */
		double eLUMO(int alpha=0) const {return _orbital_energy[alpha][_numOrb[1]];}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The matrix of overlaps. */
		vector<vector<double>> get_S();

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of f value. */
		vector<double> get_f(int orb, int alpha=0);

			//! A normal member taking one argument and returning a void value.
			/*! Actualise _all_f value. */
		void get_f(int alpha=0);

			//! A normal member taking no arguments and returning a void value.
			/*! Actualise _all_f for HOMO and LUMO Orbitals. */
		void HOMO_LUMO();

			//! A normal member taking two arguments and returning a void value.
			/*! Actualise _all_f for i and j Orbitals. */
		void HOMO_LUMO(int i, int j);

		void PrintDescriptors();

		void PrintDescriptors(int i, int j);

		vector<vector<double>>& OccupationNumber() {return _occupation_number;}

		Descriptors Descripteurs() {return _descriptors;}

		bool AlphaAndBeta() {return _alpha_and_beta;}

				//! An operator member taking two arguments and returning an ostream value.
				/*! Print all the data of two Orbitals */
		friend ostream& operator<<(ostream&, Orbitals&);
};

double operator*(const Orbitals& a, const vector<double>& coord);

#endif