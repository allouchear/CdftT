#ifndef CDFTT_WFX_H_INCLUDED
#define CDFTT_WFX_H_INCLUDED

#include<iostream>
#include<string>
#include<fstream>

using namespace std;

class NCEG
{
	private:
		string _symbol;
		vector<double> _gradients;
	public:
		NCEG();
		~NCEG(){}
		string symbol() {return _symbol;}
		vector<double> gradients() {return _gradients;}
		void add(NCEG);
};

class AEDF
{
	private:
		int _Number_of_EDF_Primitives;
		vector<int> _EDF_Primitives_Centers;
		vector<int> _EDF_Primitives_Types;
		vector<double> _EDF_Primitives_Exponents;
		vector<double> _EDF_Primitives_Coefficients;
	public:
		AEDF();
		~AEDF(){}
		int Number_of_EDF_Primitives() {return _Number_of_EDF_Primitives;}
		vector<int> EDF_Primitives_Centers() {return _EDF_Primitives_Centers;}
		vector<int> EDF_Primitives_Types() {return _EDF_Primitives_Types;}
		vector<double> EDF_Primitives_Exponents() {return _EDF_Primitives_Exponents;}
		vector<double> EDF_Primitives_Coefficients() {return _EDF_Primitives_Coefficients;}
		void add(AEDF);
};

class WFX
{
	private:
		string _Title;														
		string _Keywords;													
		int _Number_of_Nuclei;												
		int _Number_of_Primitives;	
		int _Number_of_Occupied_Molecular_Orbital;	
		int _Number_of_Perturbations;									
		vector<string> _Nuclear_Names;	
		vector<int> _Atomic_Number;
		vector<double> _Nuclear_Charges;	
		vector<vector<double>> _Nuclear_Cartesian_Coordinates;
		double _Net_Charge;													
		int _Number_of_Electrons;											
		int _Number_of_Alpha_Electrons;										
		int _Number_of_Beta_Electrons;										
		int _Electronic_Spin_Multiplicity;
		string _Model;							
		vector<int> _Primitive_Centers;										
		vector<int> _Primitive_Types;											
		vector<double> _Primitive_Exponents;									
		vector<double> _Molecular_Orbital_Occupation_Numbers;				
		vector<double> _Molecular_Orbital_Energies;
		vector<string> _Molecular_Orbital_Spin_Types;					
		vector<vector<double>> _Molecular_Orbital_Primitive_Coefficients;
		vector<int> _MO_Number;
		double _Energy; // Energy = T + Vne + Vee + Vnn
		double _Virial_Ratio; // (-V/T)
		vector<NCEG> _Nuclear_Cartesian_Energy_Gradients;
		double _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei; // ,W
		double _Full_Virial_Ratio; // ,-(V-W)/T
		int _Number_of_Core_Electrons;
		vector<AEDF> _Additionnal_Electron_Density_Function; //(EDF)
	public:
		WFX();
		~WFX(){};
		string Title() {return _Title;}
		string Keywords() {return _Keywords;}
		int Number_of_Nuclei() {return _Number_of_Nuclei;}
		int Number_of_Primitives() {return _Number_of_Primitives;}	
		int Number_of_Occupied_Molecular_Orbital() {return _Number_of_Occupied_Molecular_Orbital;}	
		int Number_of_Perturbations() {return _Number_of_Perturbations;}									
		vector<string> Nuclear_Names() {return _Nuclear_Names;}	
		vector<int> Atomic_Number() {return _Atomic_Number;}
		vector<double> Nuclear_Charges() {return _Nuclear_Charges;}	
		vector<vector<double>> Nuclear_Cartesian_Coordinates() {return _Nuclear_Cartesian_Coordinates;}
		double Net_Charge() {return _Net_Charge;}													
		int Number_of_Electrons() {return _Number_of_Electrons;}											
		int Number_of_Alpha_Electrons() {return _Number_of_Alpha_Electrons;}							
		int Number_of_Beta_Electrons() {return _Number_of_Beta_Electrons;}										
		int Electronic_Spin_Multiplicity() {return _Electronic_Spin_Multiplicity;}
		string Model() {return _Model;}							
		vector<int> Primitive_Centers() {return _Primitive_Centers;}										
		vector<int> Primitive_Types() {return _Primitive_Types;}											
		vector<double> Primitive_Exponents() {return _Primitive_Exponents;}							
		vector<double> Molecular_Orbital_Occupation_Numbers() {return _Molecular_Orbital_Occupation_Numbers;}			
		vector<double> Molecular_Orbital_Energies() {return _Molecular_Orbital_Energies;}
		vector<string> Molecular_Orbital_Spin_Types() {return _Molecular_Orbital_Spin_Types;}					
		vector<vector<double>> Molecular_Orbital_Primitive_Coefficients() {return _Molecular_Orbital_Primitive_Coefficients;}
		vector<int> MO_Number() {return _MO_Number;}
		double Energy() {return _Energy;} // Energy = T + Vne + Vee + Vnn
		double Virial_Ratio() {return _Virial_Ratio;} // (-V/T)
		vector<NCEG> Nuclear_Cartesian_Energy_Gradients() {return _Nuclear_Cartesian_Energy_Gradients;}
		double Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei() {return _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei;} // ,W
		double Full_Virial_Ratio() {return _Full_Virial_Ratio;} // ,-(V-W)/T
		int Number_of_Core_Electrons() {return _Number_of_Core_Electrons;}
		vector<AEDF> Additionnal_Electron_Density_Function() {return _Additionnal_Electron_Density_Function;} //(EDF)
};












/*
string get_one_block_from_wfx_file(istream&, string, int);
int get_one_block_int_from_wfx_file(istream&, string,  int);
double get_one_block_real_from_wfx_file(istream&, string, int);
bool get_one_int_from_wfx_file(istream&, string, int);
double get_one_orbital_from_wfx_file(istream&, int, int);
string readFile(string);
string mystrcasestr(const string, const string);
*/

#endif