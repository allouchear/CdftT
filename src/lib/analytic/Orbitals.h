#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include<iostream>
#include<analytic/WFX.h>
#include<analytic/LCAO.h>

using namespace std;

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
		double _HeV;
	public:
		Orbitals();
		Orbitals(WFX&, Binomial&);
		~Orbitals() {};
		vector<LCAO> vlcao() {return _lcao;}
		LCAO lcao(int i) {return _lcao[i];}
		int NumberOfFunctions() {return _numberOfFunctions;}
		int NumberOfAtoms() const {return _number_of_atoms;} 
		int PrimitiveCenter(int i) const {return _primitive_centers[i];}
		vector<int> PrimitiveCenters() const {return _primitive_centers;}
		vector<string> symbol() const {return _symbol;}
		double Overlap(int i, int j) {return _lcao[i].overlapLCAO(_lcao[j]);}
		void PrintOverlap(int, int);
		int HOMO() const;
		int LUMO() const;
		double eHOMO() const {return _orbital_energy[HOMO()];}
		double eLUMO() const {return _orbital_energy[LUMO()];}
		double muminus() const {return eHOMO();}
		double muplus() const {return eLUMO();}
		double mu() const {return (eHOMO()+eLUMO())/2;}
		double electronegativity() const {return -mu();}
		double hardness() const {return eLUMO()-eHOMO();}
		double softness() const {return 1/hardness();}
		double w() const {return mu()*mu()/(2*hardness());}
		double wminus() const {return muminus()*muminus()/(2*hardness());}
		double wplus() const {return muplus()*muplus()/(2*hardness());}
		double Qmax() const {return -mu()/hardness();}
		double DEmin() const {return -w();}
		vector<vector<double>> get_S() const;
		vector<double> get_f(int alpha) const;
		void get_f(vector<int>);
		vector<double> fkminus() const {return _all_f[0];}
		vector<double> fkplus() const {return _all_f[1];}
		vector<double> fk0() const;
		vector<double> Deltafk() const;
		vector<double> skminus() const;
		vector<double> skplus() const;
		vector<double> skfrac() const;
		vector<double> wkminus() const;
		vector<double> wkplus() const;
		vector<double> hardnesskminus() const;
		vector<double> hardnesskplus() const;
		vector<double> hardnessk() const;

		friend ostream& operator<<(ostream&, const Orbitals&);
};

#endif