#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include<iostream>
#include <Utils/WFX.h>
#include <Utils/FCHK.h>
#include <Utils/MOLDENGAB.h>
#include <Basis/CGTF.h>
#include <Common/PeriodicTable.h>
#include <Common/Descriptors.h>
#include <Common/Structure.h>

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
		Structure _struct;
		vector<int> _atomic_numbers;
		vector<string> _symbol;
		vector<vector<double>> _orbital_energy;
		vector<vector<double>> _all_f;
		vector<int> _numOrb;
		vector<vector<double>> _occupation_number;
		bool _alpha_and_beta;
		Binomial _bino;
		Descriptors _descriptors;
		vector<CGTF> _vcgtf_non_normalise;
		int _number_of_gtf;
		double _energy;
		vector<double> _coordinates;
		bool _mixte;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for Orbitals on 0 or "None" value. */
		Orbitals();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for Orbitals with the data in .wfx file. */
		Orbitals(WFX&, Binomial&, const PeriodicTable&);

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for Orbitals with the data in .fchk file. */
		Orbitals(FCHK&, Binomial&, const PeriodicTable&);

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for Orbitals with the data in .molden or .gab file. */
		Orbitals(MOLDENGAB&, Binomial&, const PeriodicTable&);

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for Orbitals with the data in .log file. */
		Orbitals(LOG&, Binomial&, const PeriodicTable&);

			//! A default desctructor.
			/*! We don't use it. */
		~Orbitals() {}

			//! A normal member taking no arguments and returning a vector<CGTF> value.
			/*! \return The table of CGTF which compose the Orbitals. */
		vector<CGTF>& vcgtf() {return _vcgtf;}

			//! A normal member taking no arguments and returning a vector<vector<vector<double>>> value.
			/*! [spinType][nMO][mCGTF]
			 * \return The table of coefficients of all spin-orbitals which compose the Orbitals. */
		vector<vector<vector<double>>>& coefficients() {return _coefficients;}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The table of molecular orbital energy which compose the Orbitals. */
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
			/*! \return The primitive center. */
		int PrimitiveCenter(int i) const {return _primitive_centers[i];}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The table of primitive centers. */
		vector<int> PrimitiveCenters() const {return _primitive_centers;}

			//! A normal member taking no arguments and returning a vector<string> value.
			/*! \return The table of symbol of atoms. */
		vector<string> symbol() const {return _symbol;}

			//! A normal member taking no arguments and returning a boolean value.
			/*! \return If alpha and beta are separed in the file or not (true for yes, false for no). */
		bool alphaAndBeta() const {return _alpha_and_beta;}

			//! A normal member taking no arguments and returning a void value.
			/*! Normalise all the CGTF which compose the Orbitals. */
		void NormaliseAllBasis();

			//! A normal member taking three arguments and returning a double value.
			/*! \return The value of ERI ???. */
		double ERIorbitals(Orbitals& q, Orbitals& r, Orbitals& s);

			//! A normal member taking three arguments and returning a double value.
			/*! \return The overlap between two orbitals i and j. */
		double Overlap(int i, int j, int alpha=0);

			//! A normal member taking three arguments and returning a void value.
			/*! Print the value of an overlap between two orbitals i and j. */
		void PrintOverlap(int, int, int alpha=0);

			//! A normal member taking four arguments and returning a double value.
			/*! \return The overlap between three orbitals i, j, and k. */
		double Overlap3Orbitals(int i, int j, int k, int alpha=0);

			//! A normal member taking five arguments and returning a double value.
			/*! \return The overlap between four orbitals i, j, k, and l. */
		double Overlap4Orbitals(int i, int j, int k, int l, int alpha=0);

			//! A normal member taking no arguments and returning a double value.
			/*! \return The value of the kinetic energy. */
		double kinetic();

			//! A normal member taking two arguments and returning a double value.
			/*! \return The value of the ionic potential. */
		double ionicPotential(vector<double> C, double Z);

			//! A normal member taking no arguments and returning a double value.
			/*! \return The value of the integral of Orbitals * Orbitals. */
		double OrbstarOrb();

			//! A normal member taking three arguments and returning a double value.
			/*! \return The value of the integral of Orbitals * Orbitals(ix, iy, iz) * Orbitals. */
		double OrbxyzOrb(int ix, int iy, int iz);

			//! A normal member taking no arguments and returning a void value.
			/*! Normalise the basis.
			 *  Can be delete. It was developpe for LCAO. */
		void normaliseBasis();

			//! A normal member taking three arguments and returning a double value.
			/*! \return The value of the Orbitals at the point (x,y,z). */
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

			//! A normal member taking no arguments and returning a void value.
			/*! Print all the descriptors with the FMO method. */
		void PrintDescriptors();

			//! A normal member taking no arguments and returning a void value.
			/*! Change the HOMO on the orbital i and the LUMO on the orbitals j.
			 *  Print all the descriptors with the FMO method. */
		void PrintDescriptors(int i, int j);

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! OccNum[spinType][nMO]
			 *  \return The table of occupation number. */
		vector<vector<double>>& OccupationNumber() {return _occupation_number;}

			//! A normal member taking no arguments and returning a Descriptors value.
			/*! \return the Descriptors object. */
		Descriptors Descripteurs() {return _descriptors;}

			//! A normal member taking no arguments and returning a boolean value.
			/*! \return If alpha and beta are separed in the file or not (true for yes, false for no). */
		bool AlphaAndBeta() {return _alpha_and_beta;}

			//! An operator member taking two arguments and returning an ostream value.
			/*! Print all the data of the Orbitals */
		friend ostream& operator<<(ostream&, Orbitals&);


			//! get function
			/*! returns the structure attribute*/
		Structure get_struct();
			//! Make a grid of electronic density
			/*! creates a grid of electronic density. Values are calculated with Orbitals::density(x, y, z)*/
		Grid makeGrid(const Domain& d);

			//! Electronic density
			/*! Calculates and returns the electronic density from molecular orbitals */
		double density(double x, double y, double z);
			
			//! Make a grid of Molecular orbitals
			/*! Creates a grid of molecular orbitals. Values are calculated with Orbitals::phis()*/
		Grid makeOrbGrid(const Domain& d, const vector<int>& nums, const vector<int>& typesSpin);

			//! Electronic density
			/*! Calculates and returns the electronic density from molecular orbitals */
		vector<double> phis(double x, double y, double z, const vector<int>& nums, const vector<int>& typesSpin);
			
			//! Electron localisation function
			/*! Calculates and returns the ELF*/
		double ELF(const double& x, const double& y, const double& z, double epsilon=2.87e-5);

			//! Make ELF grid
			/*! Make an ELF grid using Orbitals::ELF()*/
		Grid makeELFgrid(const Domain& d,const double& epsilon=2.87e-5);


			//! A normal member taking no arguments and returning a void value.
			/*! Sort the CGTF. 
			 *  This method need to be delete later. 
			 *  This was developped to solve a problem when you want save a .wfx file to .molden or .gab. */
		void Sorting();

			//! A normal member taking no arguments and returning a void value.
			/*! Denormalise all the CGTF which compose the Orbitals. */
		void DenormaliseAllBasis();

			//! A normal member taking one argument and returning a void value.
			/*! Save all the data in Orbitals in a file (tag = name.format). */
		void Save(string& tag);

			//! A normal member taking one argument and returning a void value.
			/*! Save all the data in Orbitals in a file .wfx. */
		void Save_wfx(string& tag);

			//! A normal member taking one argument and returning a void value.
			/*! Save all the data in Orbitals in a file .molden. */
		void Save_molden(string& tag);

			//! A normal member taking one argument and returning a void value.
			/*! Save all the data in Orbitals in a file .gab. */
		void Save_gab(string& tag);
};

	//! An operator member taking two arguments and returning a double value.
	/*! coord = (x,y,z). 
	 *  \return The value of Orbitals * Orbitals at the point (x,y,z). */
double operator*(const Orbitals& a, const vector<double>& coord);

#endif
