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
	double _dx;
	double _dy;
	double _dz;
	double _dv;

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
	
		//! Constructor
		/*!
			Calls set_N..() to set attributes and sets T and O to 0*/
	Domain(int i, int n, int m, int l, double* O);
	
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
		
		//! set() function
		/*! sets the number of values per point to N*/
	void set_Nval(int N);
	
		//! get() function
		/*!
			Returns the number of points in the X1 direction
		*/
	int N1() const;
	
		//! set() function
		/*! sets the number of points in the x direction to N*/
	void set_N1(int N);
	
		//! get() function
		/*!
			Returns the number of points in the X2 direction
		*/
	int N2() const;
	
		//! set() function
		/*! sets the number of points in the y direction to N*/
	void set_N2(int N);
	
		//! get() function
		/*!
			Returns the number of points in the X3 direction
		*/
	int N3() const;
	
		//! set() function
		/*! sets the number of points in the z direction to N*/
	void set_N3(int N);
	
		//! get() function
		/*!
			Returns the coordinates of the origin
		*/
	double* O();
	
		//! get() function
		/*!
			Returns the "translation matrix"
		*/
	vector<vector<double>> T() const;
	
		//! get() function
		/*!
			Returns Tij
		*/
	double Tij(int i , int j) const;
		
		//! set() function
		/*! sets the i j component of the translation vector to v*/	
	void set_T(double v, int i, int j);
	
		//! get() function
		/*!
			Returns the element distance in the x direction
		*/
	double dx() const;
	
		//! get() function
		/*!
			Returns the element distance in the y direction
		*/
	double dy() const;
	
		//! get() function
		/*!
			Returns the element distance in the z direction
		*/
	double dz() const;
		
		//! get() function
		/*!
			Returns the element volume
		*/
	double dv() const;
	
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
