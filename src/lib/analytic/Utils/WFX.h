#ifndef CDFTT_WFX_H_INCLUDED
#define CDFTT_WFX_H_INCLUDED

#include<iostream>
#include<vector>
#include<string>
#include<fstream>
#include<analytic/Utils/Utils.h>

using namespace std;

	//! A MOPC class.
	/*! This class will be used in the WFX class. */

class MOPC
{
	private:
		int _MO_Number;
		vector<double> _Coefficients;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters on 0 or "None" value. */

		MOPC();

			//! A default desctructor.
			/*! We don't use it. */

		~MOPC(){}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The MO number. */

		int MO_Number() {return _MO_Number;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The MO coefficients. */

		vector<double> Coefficients() {return _Coefficients;}

			//! A normal member taking two arguments and returning a void value.
			/*! Push back one MO number and their coefficients */

		void push_back(int, vector<double>);

};

	//! A NCEG class.
	/*! This class will be used in the WFX class. */

class NCEG
{
	private:
		string _symbol;
		vector<double> _gradient;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters on 0 or "None" value. */

		NCEG();

			//! A default desctructor.
			/*! We don't use it. */

		~NCEG(){}

			//! A normal member taking no arguments and returning a string value.
			/*! \return The symbol of the atom. */

		string symbol() {return _symbol;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The gradient of one atom. */

		vector<double> gradient() {return _gradient;}

			//! A normal member taking two arguments and returning a void value.
			/*! Push back the symbol of an atom and their gradient */

		void push_back(string, vector<double>);
};

	//! A AEDF class.
	/*! This class will be used in the WFX class. */

class AEDF
{
	private:
		int _Number_of_EDF_Primitives;
		vector<int> _EDF_Primitives_Centers;
		vector<int> _EDF_Primitives_Types;
		vector<double> _EDF_Primitives_Exponents;
		vector<double> _EDF_Primitives_Coefficients;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters on 0 or "None" value. */

		AEDF();

			//! A default desctructor.
			/*! We don't use it. */

		~AEDF(){}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of EDF primitives. */

		int Number_of_EDF_Primitives() {return _Number_of_EDF_Primitives;}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The EDF primitives centers. */

		vector<int> EDF_Primitives_Centers() {return _EDF_Primitives_Centers;}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The EDF primitives types. */

		vector<int> EDF_Primitives_Types() {return _EDF_Primitives_Types;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The EDF primitives exponents. */

		vector<double> EDF_Primitives_Exponents() {return _EDF_Primitives_Exponents;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The EDF primitives coefficients. */

		vector<double> EDF_Primitives_Coefficients() {return _EDF_Primitives_Coefficients;}

			//! A normal member taking five arguments and returning a void value.
			/*! Push back the number, the centers, the types, the expoenents and the coefficients of EDF primitives. */

		void push_back(int, vector<int>, vector<int>, vector<double>, vector<double>);
};

	//! A WFX class.
	/*! This class will be used to read and write in the wfx format. */

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
		vector<double> _Nuclear_Cartesian_Coordinates;
		double _Net_Charge;													
		int _Number_of_Electrons;											
		int _Number_of_Alpha_Electrons;										
		int _Number_of_Beta_Electrons;										
		int _Electronic_Spin_Multiplicity;
		string _Model;							
		vector<int> _Primitive_Centers;										
		vector<int> _Primitive_Types;
		vector<vector<int>> _Lxyz;											
		vector<double> _Primitive_Exponents;									
		vector<vector<double>> _Molecular_Orbital_Occupation_Numbers;				
		vector<vector<double>> _Molecular_Orbital_Energies;
		vector<string> _Molecular_Orbital_Spin_Types;					
		vector<vector<MOPC>> _Molecular_Orbital_Primitive_Coefficients;
		double _Energy; // Energy = T + Vne + Vee + Vnn
		double _Virial_Ratio; // (-V/T)
		vector<NCEG> _Nuclear_Cartesian_Energy_Gradients;
		double _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei; // ,W
		double _Full_Virial_Ratio; // ,-(V-W)/T
		int _Number_of_Core_Electrons;
		AEDF _Additionnal_Electron_Density_Function; //(EDF)

		bool _alpha_and_beta;
	public:

			//! A default constructor.
			/*! This constructor is used to set all of the parameters on 0 or "None" value. */

		WFX();

			//! A constructor.
			/*! This constructor is used to set all of the parameters with the data in the file. */

		WFX(ifstream&);

			//! A default desctructor.
			/*! We don't use it. */

		~WFX(){};

			//! A normal member taking no arguments and returning a string value.
			/*! \return The title. */

		string Title() {return _Title;}

			//! A normal member taking no arguments and returning a string value.
			/*! \return The keywords. */

		string Keywords() {return _Keywords;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of nuclei. */

		int Number_of_Nuclei() {return _Number_of_Nuclei;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of primitives. */

		int Number_of_Primitives() {return _Number_of_Primitives;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of occupied molecular orbital. */

		int Number_of_Occupied_Molecular_Orbital() {return _Number_of_Occupied_Molecular_Orbital;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of perturbations. */

		int Number_of_Perturbations() {return _Number_of_Perturbations;}	

			//! A normal member taking no arguments and returning a vector<string> value.
			/*! \return The table of nuclear names. */

		vector<string> Nuclear_Names() {return _Nuclear_Names;}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The table of atomic number. */

		vector<int> Atomic_Number() {return _Atomic_Number;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of nuclear charges. */

		vector<double> Nuclear_Charges() {return _Nuclear_Charges;}	

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of nuclear cartesian coordinates. */

		vector<double> Nuclear_Cartesian_Coordinates() {return _Nuclear_Cartesian_Coordinates;}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The net charge. */

		double Net_Charge() {return _Net_Charge;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of electrons. */

		int Number_of_Electrons() {return _Number_of_Electrons;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of alpha electrons. */

		int Number_of_Alpha_Electrons() {return _Number_of_Alpha_Electrons;}							
		
			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of beta electrons. */

		int Number_of_Beta_Electrons() {return _Number_of_Beta_Electrons;}

			//! A normal member taking no arguments and returning an int value.
			/*! \return The electronic spin multiplicity. */

		int Electronic_Spin_Multiplicity() {return _Electronic_Spin_Multiplicity;}

			//! A normal member taking no arguments and returning a string value.
			/*! \return The model. */

		string Model() {return _Model;}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The table of primitives centers. */

		vector<int> Primitive_Centers() {return _Primitive_Centers;}

			//! A normal member taking no arguments and returning a vector<int> value.
			/*! \return The table of primitives types. */

		vector<int> Primitive_Types() {return _Primitive_Types;}

			//! A normal member taking no arguments and returning a vector<vector<int>> value.
			/*! \return The table of Lx, Ly, and Lz values for each primitives. */

		vector<vector<int>> Lxyz() {return _Lxyz;}	

			//! A normal member taking one argument and returning a vector<int> value.
			/*! \return The table of Lx, Ly, and Lz values for one primitive. */

		vector<int> Lxyz(int i) {return _Lxyz[i];}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of primitive exponents. */

		vector<double> Primitive_Exponents() {return _Primitive_Exponents;}			

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of molecular orbital occupation numbers. */

		vector<vector<double>> Molecular_Orbital_Occupation_Numbers() {return _Molecular_Orbital_Occupation_Numbers;}

			//! A normal member taking no arguments and returning a vector<double> value.
			/*! \return The table of molecular orbital energies. */

		vector<vector<double>> Molecular_Orbital_Energies() {return _Molecular_Orbital_Energies;}

			//! A normal member taking no arguments and returning a vector<string> value.
			/*! \return The table of molecular orbital spin types. */

		vector<string> Molecular_Orbital_Spin_Types() {return _Molecular_Orbital_Spin_Types;}

			//! A normal member taking no arguments and returning a vector<MOPC> value.
			/*! \return The table of molecular orbital primitive coefficients. */

		vector<vector<MOPC>> Molecular_Orbital_Primitive_Coefficients() {return _Molecular_Orbital_Primitive_Coefficients;}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The energy. */

		double Energy() {return _Energy;} // Energy = T + Vne + Vee + Vnn

			//! A normal member taking no arguments and returning a double value.
			/*! \return The viriral ratio. */

		double Virial_Ratio() {return _Virial_Ratio;} // (-V/T)

			//! A normal member taking no arguments and returning a vector<NCEG> value.
			/*! \return The table of nuclear cartesian energy gradients. */

		vector<NCEG> Nuclear_Cartesian_Energy_Gradients() {return _Nuclear_Cartesian_Energy_Gradients;}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The nuclear virial of energy gradient based forces on nuclei. */

		double Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei() {return _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei;} // ,W

			//! A normal member taking no arguments and returning a double value.
			/*! \return The full virial ratio. */

		double Full_Virial_Ratio() {return _Full_Virial_Ratio;} // ,-(V-W)/T

			//! A normal member taking no arguments and returning an int value.
			/*! \return The number of core electrons. */

		int Number_of_Core_Electrons() {return _Number_of_Core_Electrons;}

			//! A normal member taking no arguments and returning an AEDF value.
			/*! \return The additionnal electron density function. */

		AEDF Additionnal_Electron_Density_Function() {return _Additionnal_Electron_Density_Function;} //(EDF)

		bool AlphaAndBeta() {return _alpha_and_beta;}


			//! A normal member taking four arguments and returning a vector<int> value.
			/*! \return The one block of int read. */

		vector<int> read_one_block_int(ifstream&, string, bool, int);

			//! A normal member taking three arguments and returning a vector<double> value.
			/*! \return The one block of real read. */

		vector<double> read_one_block_real(ifstream&, string, bool);

			//! A normal member taking four arguments and returning a vector<double> value.
			/*! \return The one block of real read. */

		vector<double> read_one_block_real(ifstream&, string, bool, int);

			//! A normal member taking three arguments and returning a vector<string> value.
			/*! \return The one block of string read. */

		vector<string> read_one_block_string(ifstream&, string, bool);

			//! A normal member taking three arguments and returning an int value.
			/*! \return One int read. */

		int read_int(ifstream&, string, bool);

			//! A normal member taking three arguments and returning a double value.
			/*! \return One real read. */

		double read_real(ifstream&, string, bool);

			//! A normal member taking three arguments and returning a string value.
			/*! \return One string read. */

		string read_string(ifstream&, string, bool);

			//! A normal member taking three arguments and returning a vector<MOPC> value.
			/*! \return The MOPC block read. */

		vector<MOPC> read_MOPC_block(ifstream&, string, bool);

			//! A normal member taking three arguments and returning a vector<NCEG> value.
			/*! \return The NCEG block read. */

		vector<NCEG> read_NCEG_block(ifstream&, string, bool);

			//! A normal member taking three arguments and returning an AEDF value.
			/*! \return The AEDF block read. */

		AEDF read_AEDF_block(ifstream&, string, bool);

			//! A normal member taking three arguments and returning a void value.
			/*! Read the file and set parameters on it. */

		void read_file_wfx(ifstream&);


			//! A normal member taking five arguments and returning a void value.
			/*! Write one block of int */

		void write_one_block_int(ofstream&, vector<int>, string, bool, int);

			//! A normal member taking four arguments and returning a void value.
			/*! Write one block of real */

		void write_one_block_real(ofstream&, vector<double>, string, bool);

			//! A normal member taking five arguments and returning a void value.
			/*! Write one block of real */

		void write_one_block_real(ofstream&, vector<double>, string, bool, int);

			//! A normal member taking four arguments and returning a void value.
			/*! Write one block of string */

		void write_one_block_string(ofstream&, vector<string>, string, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write one int */

		void write_one_matrix_real(ofstream&, vector<vector<double>>, string, bool);

		void write_int(ofstream&, int, string, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write one real */

		void write_real(ofstream&, double, string, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write one string */

		void write_string(ofstream&, string ,string, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write MOPC block */

		void write_MOPC_block(ofstream&, vector<vector<MOPC>>, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write NCEG block */

		void write_NCEG_block(ofstream&, vector<NCEG>, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write AEDF block */

		void write_AEDF_block(ofstream&, AEDF, bool);

			//! A normal member taking four arguments and returning a void value.
			/*! Write all of the parameters in a wfx file */

		void write_file_wfx(ofstream&);
};

	//! A function taking three arguments and returning a long int value.
	/*! \return The position of a block and the number of elements in. */

long int LocaliseBlock(ifstream&, int&, string);

	//! A function taking four arguments and returning a long int value.
	/*! \return The position of a MO block, their number, and the number of elements in. */

long int LocaliseMO(ifstream&, int&, int&, string);

	//! A function taking one argument and returning a vector<int> value.
	/*! \return The Lx, Ly, and Lz values of one primitive type*/

vector<int> setLxyz(int);


#endif