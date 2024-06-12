#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include<iostream>
#include<analytic/Utils/WFX.h>
#include<analytic/Orbitals/LCAO.h>
#include<common/PeriodicTable.h>
#include<common/Descriptors.h>

using namespace std;


	//! An Orbitals class.
	/*! This class will be used to calculate descriptors. */
class Orbitals
{
	private:
		vector<LCAO> _lcao;
		int _numberOfFunctions;
		int _number_of_alpha_electrons;
		int _number_of_beta_electrons;
		int _number_of_atoms;
		vector<int> _primitive_centers;
		vector<string> _symbol;
		vector<double> _orbital_energy;
		vector<vector<double>> _all_f;
		vector<int> _numOrb;
		vector<double> _occupation_number;
		Descriptors _descriptors;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one LCAO on 0 or "None" value. */
		Orbitals();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for one LCAO. */
		Orbitals(WFX&, Binomial&, const PeriodicTable&);

			//! A default desctructor.
			/*! We don't use it. */
		~Orbitals() {};

			//! A normal member taking no arguments and returning a vector<LCAO> value.
			/*! \return The table of LCAO which compose the Orbitals. */
		vector<LCAO> vlcao() {return _lcao;}

			//! A normal member taking one argument and returning a LCAO value.
			/*! \return The LCAO i which compose the Orbitals. */
		LCAO lcao(int i) {return _lcao[i];}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of LCAO in the Orbitals. */
		int NumberOfFunctions() {return _numberOfFunctions;}

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

			//! A normal member taking two arguments and returning a double value.
			/*! \return The overlap between two LCAO which compose the Orbitals. */
		double Overlap(int i, int j) {return _lcao[i].overlapLCAO(_lcao[j]);}

			//! A normal member taking two arguments and returning a void value.
			/*! Print the value of an overlap. */
		void PrintOverlap(int, int);

			//! A normal member taking no arguments and returning an int value.
			/*! \return The HOMO Molecular Orbital number. */
		void HOMO();

			//! A normal member taking no arguments and returning an int value.
			/*! \return The LUMO Molecular Orbital number. */
		void LUMO();

			//! A normal member taking no arguments and returning a double value.
			/*! \return The HOMO's energy. */
		double eHOMO() const {return _orbital_energy[_numOrb[0]];}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The LUMO's energy. */
		double eLUMO() const {return _orbital_energy[_numOrb[1]];}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The matrix of overlaps. */
		vector<vector<double>> get_S() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of f value. */
		vector<double> get_f(int alpha) const;

			//! A normal member taking one argument and returning a void value.
			/*! Actualise _all_f value. */
		void get_f();

			//! A normal member taking no arguments and returning a void value.
			/*! Actualise _all_f for HOMO and LUMO Orbitals. */
		void HOMO_LUMO();

			//! A normal member taking two arguments and returning a void value.
			/*! Actualise _all_f for i and j Orbitals. */
		void HOMO_LUMO(int i, int j);

		void PrintDescriptors();

		void PrintDescriptors(int i, int j);

		vector<double> OccupationNumber() {return _occupation_number;}

		Descriptors Descripteurs() {return _descriptors;}

				//! An operator member taking two arguments and returning an ostream value.
				/*! Print all the data of two Orbitals */
		friend ostream& operator<<(ostream&, const Orbitals&);
};

#endif