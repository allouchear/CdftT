#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include<iostream>
#include <Utils/WFX.h>
#include <Utils/FCHK.h>
#include <Utils/MOLDENGAB.h>
#include <Basis/CGTF.h>
#include <Common/PeriodicTable.h>
#include <Common/Descriptors.h>

using namespace std;


	//! An Orbitals class.
	/*! This class will be used to calculate descriptors. */
class Orbitals
{
	private:
		vector<CGTF> _vcgtf;
		vector<CGTF> _vcgtf_non_normalise;
		vector<vector<vector<double>>> _coefficients;
		int _numberOfAo;
		int _numberOfMo;
		int _number_of_alpha_electrons;
		int _number_of_beta_electrons;
		int _number_of_atoms;
		int _number_of_gtf;
		vector<int> _primitive_centers;
		vector<int> _atomic_numbers;
		vector<double> _coordinates;
		vector<string> _symbol;
		vector<vector<double>> _orbital_energy;
		vector<vector<double>> _all_f;
		vector<int> _numOrb;
		vector<vector<double>> _occupation_number;
		double _energy;
		bool _alpha_and_beta;
		Binomial _bino;
		Descriptors _descriptors;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for Orbitals on 0 or "None" value. */
		Orbitals();

			//! Constructor.
			/*! This constructor is used to add all of the parameters for Orbitals from a .wfx file. */
		Orbitals(WFX&, Binomial&, const PeriodicTable&);

			//! Constructor.
			/*! This constructor is used to add all of the parameters for Orbitals from a .fchk file. */
		Orbitals(FCHK&, Binomial&, const PeriodicTable&);

			//! Constructor.
			/*! This constructor is used to add all of the parameters for Orbitals from a .molden or a .gab file. */
		Orbitals(MOLDENGAB&, Binomial&, const PeriodicTable&);

			//! Constructor.
			/*! This constructor is used to add all of the parameters for Orbitals from a .log file. */
		Orbitals(LOG&, Binomial&, const PeriodicTable&);

			//! A default desctructor.
			/*! We don't use it. */
		~Orbitals() {}

			//! A normal member taking no arguments and returning a vector<CGTF> value.
			/*! \return The table of CGTF which compose the Orbitals. */
		vector<CGTF>& vcgtf() {return _vcgtf;}

			//! A normal member taking no arguments and returning a vector<vector<vector<double>>> value.
			/*! \return The table of all coefficients (alpha/beta, number of molecular orbitals, and number of CGTF) which compose the Orbitals. */
		vector<vector<vector<double>>>& coefficients() {return _coefficients;}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The table of energy (alpha/beta, number of molecular orbitals) of each molecular orbital. */
		vector<vector<double>>& Energy()  {return _orbital_energy;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of atomic orbitals in the Orbitals. */
		int NumberOfAo() {return _numberOfAo;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of molecular orbitals in the Orbitals. */
		int NumberOfMo() {return _numberOfMo;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of atoms. */
		int NumberOfAtoms() const {return _number_of_atoms;} 

			//! A normal member taking one argument and returning an int value.
			/*! \return The primitive center i. */
		int PrimitiveCenter(int i) const {return _primitive_centers[i];}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The table of primitive centers. */
		vector<int> PrimitiveCenters() const {return _primitive_centers;}

			//! A normal member taking no arguments and returning a vector<string> value.
			/*! \return The table of symbol of atoms. */
		vector<string> symbol() const {return _symbol;}

			//! A normal member taking no arguments and returning a void value.
			/*! Normalise all the basis. */
		void NormaliseAllBasis();

			//! A normal member taking no arguments and returning a void value.
			/*! Unnormalise all the basis. */
		void DenormaliseAllBasis();

			//! A normal member taking three arguments and returning a double value.
			/*! Return the ERI value between three Orbitals ??? */
		double ERIorbitals(Orbitals& q, Orbitals& r, Orbitals& s);

			//! A normal member taking three arguments and returning a double value.
			/*! \return The overlap between two Orbitals i and j. */
		double Overlap(int i, int j, int alpha=0);

			//! A normal member taking three arguments and returning a void value.
			/*! Calculate and print the value of an overlap between two Orbitals i and j. */
		void PrintOverlap(int, int, int alpha=0);

			//! A normal member taking four arguments and returning a double value.
			/*! \return The overlap between three Orbitals i, j, and k. */
		double Overlap3Orbitals(int i, int j, int k, int alpha=0);

			//! A normal member taking five arguments and returning a double value.
			/*! \return The overlap between four Orbitals i, j, k, and l. */
		double Overlap4Orbitals(int i, int j, int k, int l, int alpha=0);

			//! A normal member taking no arguments and returning a double value.
			/*! \return The total kinetic value of Orbitals. */
		double kinetic();

			//! A normal member taking no arguments and returning a double value.
			/*! \return The total ionic potential value of Orbitals. */
		double ionicPotential(vector<double> C, double Z);

			//! A normal member taking no arguments and returning a double value.
			/*! \return The integral value of Orbitals ??? */
		double OrbstarOrb();

			//! A normal member taking three arguments and returning a double value.
			/*! \return The integral value of Orbitals ??? */
		double OrbxyzOrb(int ix, int iy, int iz);

			//! A normal member taking three arguments and returning a double value.
			/*! \return The value of Orbitals at the coordinates x,y,z. */
		double func(double x, double y, double z) const;

			//! A normal member taking no arguments and returning an int value.
			/*! \return The HOMO Molecular Orbital number. */
		void HOMO();

			//! A normal member taking no arguments and returning an int value.
			/*! \return The LUMO Molecular Orbital number. */
		void LUMO();

			//! A normal member taking one argument and returning a double value.
			/*! \return The HOMO's energy. */
		double eHOMO(int alpha=0) const {return _orbital_energy[alpha][_numOrb[0]];}

			//! A normal member taking one argument and returning a double value.
			/*! \return The LUMO's energy. */
		double eLUMO(int alpha=0) const {return _orbital_energy[alpha][_numOrb[1]];}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The matrix of overlaps. */
		vector<vector<double>> get_S();

			//! A normal member taking two arguments and returning a vector<double> value.
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

			//! A normal member taking no arguments and returning a void value.
			/*! Print all the descriptors with the FMO method. */
		void PrintDescriptors();

			//! A normal member taking two arguments and returning a void value.
			/*! Select HOMO (orbital i) and LUMO (orbital j) and print all the descriptors with the FMO method. */
		void PrintDescriptors(int i, int j);

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The table of occupation (alpha/beta, number of molecular orbital) for each molecular orbital. */
		vector<vector<double>>& OccupationNumber() {return _occupation_number;}

			//! A normal member taking no arguments and returning a Descriptors value.
			/*! \return The Descriptors object. */
		Descriptors Descripteurs() {return _descriptors;}

			//! A normal member taking no arguments and returning a bool value.
			/*! \return Boolean of alpha and beta (true if alpha and beta electrons have the same coefficients in the read file). */
		bool AlphaAndBeta() {return _alpha_and_beta;}

			//! An operator member taking two arguments and returning an ostream value.
			/*! Print all the data of two Orbitals */
		friend ostream& operator<<(ostream&, Orbitals&);

			//! A normal member taking one argument and returning a void value.
			/*! Save the data for the format choose in the string with the name (example : name.format) 
			 *  Call the appropriate member for the format choose. */
		void Save(string&);

			//! A normal member taking one argument and returning a void value.
			/*! Save the data in a .wfx file */
		void Save_wfx(string&);

			//! A normal member taking one argument and returning a void value.
			/*! Save the data in a .molden file */
		void Save_molden(string&);

			//! A normal member taking one argument and returning a void value.
			/*! Save the data in a .gab file */
		void Save_gab(string&);

			//! A normal member taking no arguments and returning a void value.
			/*! Sort the Orbitals because the format wfx is completely different from other and can swap of place some molecular orbitals. */
		void Sorting();
};

	//! An operator member taking two arguments and returning a double value.
	/*! \return The double value of an Orbitals at the coordinates x,y,z*/
double operator*(const Orbitals& a, const vector<double>& coord);

#endif