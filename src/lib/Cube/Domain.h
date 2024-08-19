#ifndef _CDFTT_DOMAIN_H_INCLUDED
#define _CDFTT_DOMAIN_H_INCLUDED

using namespace std;
#include <fstream>
#include <vector>
#include <Common/Structure.h>
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
	vector<vector<double>> _inv_T;
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
	Domain(int Nval, int N1, int N2,int N3, double xmax, double ymax, double zmax, vector<vector<double>> T);
	
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
		
		//! x
		/*! returns the x value in space of the point i, j ,k */
	double x(int i, int j, int k) const;
	
		//! y
		/*! returns the y value in space of the point i, j ,k */
	double y(int i, int j, int k) const;
	
		//! z
		/*! returns the z value in space of the point i, j ,k */
	double z(int i, int j, int k) const;
	
		//! i
		/*! returns the i value on grid of a point in space x, y, z. Rounds up*/
	int i(double x, double y, double z) const;
	
		//! j
		/*! returns the j value on grid of a point in space x, y, z. Rounds up*/
	int j(double x, double y, double z) const;

		//! k
		/*! returns the k value on grid of a point in space x, y, z. Rounds up*/
	int k(double x, double y, double z) const;
	
		//! Inverse T
		/*! Inverts the translation matrix. used to init _invT*/
	void inverse_T();
		
		//! Measure molecule
		/*! Returns the maximum length betweens the atoms of S. best to use Grid::sizeUpMol()*/
	double sizeUpMol(const Structure& S, double scale);
		
		//! Sets all
		/*! Sets all the attributes of Domain*/
	void set_all(int Nval, int N1, int N2,int N3, double xmax, double ymax, double zmax, vector<vector<double>> T);
};

#endif //_CDFTT_DOMAIN_H_INCLUDED
