#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED
 
/*
------------------------------------------------------------------------------------------------------------------------------
 Energies (hardness, mu, w, Xi, DEmin, wk-, wk+, hardnessk-, hardnessk+, hardnessk) are given in eV
 Softnesses (S, sk-, sk+) are given in eV^-1
------------------------------------------------------------------------------------------------------------------------------
   mu-         = eHOMO
   mu+         = eLUMO
   mu          = Chemical potential = (mu+ + mu-)/2
   hardness    = Chemical hardness = (mu+  -  mu-)
   Xi          = Electronegativity = -mu
   w           = Electrophilicity index = mu^2/(2 hardness) 
   w-          = propensity to donate electron = mu-^2/(2 hardness) 
   w+          = propensity to accept electron = mu+^2/(2 hardness) 
   S           = Global softness = 1/hardness
   Qmax        = Maximal electronic charge accepted by an electrophile = -mu/hardness
   DEmin       = Energy decrease if the electrophile take Qmax = -mu^2/(2 hardness) 
   fk-         = Local Fukui electrophilic attack
   fk+         = Local Fukui nucleophilic attack
   sk-         = Local softness electrophilic attack = S fk-
   sk+         = Local softness nucleophilic attack = S fk+
   wk-         = Local philicity index of electrophilic attack = w fk-
   wk+         = Local philicity index of nucleophilic attack = w fk+
   hardnessk-  = Local hardness = mu+ fk+ - mu- fk- - (mu+- mu-)*(fk+-fk-)
   hardnessk+  = Local hardness = mu+ fk+ - mu- fk- + (mu+- mu-)*(fk+-fk-)
   hardnessk   = Local hardness = mu+ fk+ - mu- fk-
   Deltafk     = Dual descripor = (fk+ - fk-) : 
                 >0 => site favored for a nucleophilic attack
                 <0 => site favored for an electrophilic attack
------------------------------------------------------------------------------------------------------------------------------
 References:
  - Revisiting the definition of local hardness and hardness kernel 
    C. A. Polanco-Ramrez et al
    Phys. Chem. Chem. Phys., 2017, 19, 12355-12364
    DOI: 10.1039/c7cp00691h
  - Applications of the Conceptual Density Functional Theory 
    Indices to Organic Chemistry Reactivity
    Luis R. Domingo, Mar Ríos-Gutiérrez and Patricia Pérez 
    Molecules 2016, 21, 748; doi:10.3390/molecules21060748
  - Electrodonating and Electroaccepting Powers
    José L. Gazquez, André Cedillo, and Alberto Vela
    J. Phys. Chem. A 2007, 111, 1966-1970, DOI: 10.1021/jp065459f
  - Introducing “UCA-FUKUI” software: reactivity-index calculations
    Jesús Sánchez-Márquez et al.
    J Mol Model (2014) 20:2492, DOI 10.1007/s00894-014-2492-1
  - Dual descriptor and molecular electrostatic potential: 
    complementary tools for the study of the coordination 
    chemistry of ambiphilic ligands
    F.  Guégan et al.
    Phys.Chem.Chem.Phys., 2014, 16 , 15558-15569, 
    DOI: 10.1039/c4cp01613k
  - New Dual Descriptor for Chemical Reactivity
    Ch. Morell et al.
    J. Phys. Chem. A 2005, 109, 205-212, DOI: 10.1021/jp046577a
------------------------------------------------------------------------------------------------------------------------------
*/

#include<iostream>
#include<analytic/WFX.h>
#include<analytic/LCAO.h>

using namespace std;

	//! A Density class.
	/*! This class will be used to calculate the density. */
class Density
{
	private:
		vector<GTF> _gtf;
		vector<double> _occupation_number;
		vector<vector<double>> _orbital_coefficients;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one LCAO on 0 or "None" value. */
		Density();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters to calculate the density. */
		Density(WFX&, Binomial&);

			//! A default desctructor.
			/*! We don't use it. */
		~Density() {}

			//! A normal member taking no arguments and returning a vector<GTF> value.
			/*! \return The table of GTF which be use to calculate the density. */
		vector<GTF> gtf() {return _gtf;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of occupation number which be use to calculate the Orbitals. */
		vector<double> OccupationNumber() {return _occupation_number;}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The table of orbital coefficients which be use to calculate the Orbitals. */
		vector<vector<double>> OrbitalCoefficients() {return _orbital_coefficients;}
};

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
		Density _density;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters for one LCAO on 0 or "None" value. */
		Orbitals();

			//! A real constructor.
			/*! This constructor is used to add all of the parameters for one LCAO. */
		Orbitals(WFX&, Binomial&);

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

			//! A normal member taking no arguments and returning a double value.
			/*! \return The mu- value. */
		double muminus() const {return eHOMO();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The mu+ value. */
		double muplus() const {return eLUMO();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The mu value. */
		double mu() const {return (eHOMO()+eLUMO())/2;}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The electronegativity value. */
		double electronegativity() const {return -mu();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The hardness value. */
		double hardness() const {return eLUMO()-eHOMO();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The softness value. */
		double softness() const {return 1/hardness();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The omega value. */
		double w() const {return mu()*mu()/(2*hardness());}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The omega- value. */
		double wminus() const {return muminus()*muminus()/(2*hardness());}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The omega+ value. */
		double wplus() const {return muplus()*muplus()/(2*hardness());}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The Qmax value. */
		double Qmax() const {return -mu()/hardness();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The Delta E min value. */
		double DEmin() const {return -w();}

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return The matrix of overlaps. */
		vector<vector<double>> get_S() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of f value. */
		vector<double> get_f(int alpha) const;

			//! A normal member taking one argument and returning a void value.
			/*! Actualise _all_f value. */
		void get_f();

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of fk- value. */
		vector<double> fkminus() const {return _all_f[0];}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of fk+ value. */
		vector<double> fkplus() const {return _all_f[1];}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of fk0 value. */
		vector<double> fk0() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of Delta fk value. */
		vector<double> Deltafk() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of sk- value. */
		vector<double> skminus() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of sk+ value. */
		vector<double> skplus() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of sk-/sk+ value. */
		vector<double> skfrac() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of omegak- value. */
		vector<double> wkminus() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of omegak- value. */
		vector<double> wkplus() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of hardnessk- value. */
		vector<double> hardnesskminus() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of hardnessk+ value. */
		vector<double> hardnesskplus() const;

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of hardnessk value. */
		vector<double> hardnessk() const;

			//! A normal member taking no arguments and returning a void value.
			/*! Actualise _all_f for HOMO and LUMO Orbitals. */
		void HOMO_LUMO();

			//! A normal member taking two arguments and returning a void value.
			/*! Actualise _all_f for i and j Orbitals. */
		void HOMO_LUMO(int i, int j);

				//! An operator member taking two arguments and returning an ostream value.
				/*! Print all the data of two Orbitals */
		friend ostream& operator<<(ostream&, const Orbitals&);
};

#endif