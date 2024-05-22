#ifndef _CDFTT_DOMAIN_H_INCLUDED
#define _CDFTT_DOMAIN_H_INCLUDED

using namespace std;
#include <fstream>
#include <vector>

//! Domain class
/*! This class is used to define the domain in 3D space containing a structure, molecule,...*/
class Domain
{	
	int _Nval;
	int _N1;
	int _N2;
	int _N3;
	double _O[3];
	vector<vector<double>> _T;

	public:
		
		//! Default Constructor
		/*!
			Sets all attributes to 0 except _Nval which is 1 by default.
		*/	
	Domain();
	
		//! Constructor
		/*!
			Calls read_From_Cube to construct new Domain
		*/
	Domain(ifstream& nameFile);
		//! read .cube file
		/*!
			Reads from a .cube file and initializes the number of atoms, the geometry, etc...
		*/
	void read_From_Cube(ifstream& nameFile);
	
		//! get() function
		/*!
			Returns the number of of values per point Nval
		*/
	int Nval() const;
	
		//! get() function
		/*!
			Returns the number of points in the X1 direction
		*/
	int N1() const;
	
		//! get() function
		/*!
			Returns the number of points in the X2 direction
		*/
	int N2() const;
	
		//! get() function
		/*!
			Returns the number of points in the X3 direction
		*/
	int N3() const;
	
		//! get() function
		/*!
			Returns the coordinates of the origin
		*/
	double* O();
	
		//! get() function
		/*!
			Returns the "translation matrix??"
		*/
	vector<vector<double>> T() const;
	
		//!Operator ==
		/*! 
			Overload of operator ==
		*/
	bool operator==(const Domain& D) const;
	
		//!Operator ==
		/*! 
			Overload of operator ==
		*/
	bool operator!=(const Domain& D) const; 
};

#endif //_CDFTT_DOMAIN_H_INCLUDED
