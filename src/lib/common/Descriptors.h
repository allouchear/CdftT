#ifndef CDFTT_DESCRIPTORS_H_INCLUDED
#define CDFTT_DESCRIPTORS_H_INCLUDED
/*
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
#include <numeric/Grid.h>
#include <numeric/GridCP.h>
#include <vector>
#include <common/Structure.h>

using namespace std;


class Descriptors
{	
	private:

		Structure _str;

		/******************************************GLOBAL DESCRIPTORS**************************************************/
		double _mu;
		double _mup;
		double _mum;
		
		double _xi;
		double _hardness;
		double _w;
		double _wp;
		double _wm;
		double _S;
		double _Qmax;
		double _DEmin;
		
		/******************************************LOCAL DESCRIPTORS**************************************************/
		vector<double> _fk0;
		vector<double> _fkm;
		vector<double> _fkp;
		
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
		void reset();
		Descriptors();
		Descriptors(GridCP& g0, vector<double> Q0, vector<double> Qm, vector<double> Qp);
		~Descriptors() {}

		void compute_all();
	

	
			//! Mu, Mum, Mup
			/*! Sets the values of mu given the ionisation potential and the elctron affinity*/
		void set_all_mu(const double& I,const double& A);
	
			//! fukui
			/*! Calculates and sets the values of the fukui functions*/
		void compute_fk_From_Charge( vector<double> Q0, vector<double> Qm, vector<double> Qp);
	
		/********************************************************************************************/
		friend ostream& operator<<(ostream& flux, const Descriptors&);
};

#endif
