#ifndef CDFTT_DESCRIPTORS_H_INCLUDED
#define CDFTT_DESCRIPTORS_H_INCLUDED

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

#include <vector>
#include <common/Structure.h>

using namespace std;

class Descriptors
{	
	private:

		Structure _str;

		vector<double> _Qp;
		vector<double> _Q0;
		vector<double> _Qm;
	
		double _mu;
		double _mup;
		double _mum;
	
		vector<double> _fk0;
		vector<double> _fkm;
		vector<double> _fkp;

		double _xi;
		double _hardness;
		double _w;
		double _wp;
		double _wm;
		double _S;
		double _Qmax;
		double _DEmin;
		vector<double> _Deltafk;
		vector<double> _wkm;
		vector<double> _wkp;
		vector<double> _Skm;
		vector<double> _Skp;
		vector<double> _Skfrac;
		vector<double> _hardnessk;
		vector<double> _hardnesskm;
		vector<double> _hardnesskp;
	
	public:

		Descriptors();
		Descriptors(WFX&, const PeriodicTable&);
		~Descriptors() {}

		void compute_all();
		
			//! Mu, Mum, Mup
			/*! Sets the values of mu given the ionisation potential and the elctron affinity*/
		void set_all_mu(const double& I,const double& A);
		
			//! fukui
			/*! Calculates and sets the values of the fukui functions*/
		void compute_fk();

		void set_mu_fk_data(vector<vector<double>>, double, double);
		void set_mu_fk_data(vector<vector<double>>);
	
		/********************************************************************************************/
		friend ostream& operator<<(ostream& flux, const Descriptors&);
};

#endif