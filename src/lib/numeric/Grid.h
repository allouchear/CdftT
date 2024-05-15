#ifndef _CDFTT_GRID_H_INCLUDED
#define _CDFTT_GRID_H_INCLUDED

using namespace std;
#include <numeric/Domain.h>
#include <common/Structure.h>


class Grid
{
	Domain _dom;
	Structure _str;
	vector<vector<vector<vector<double>>>> _V;
	
	public:
		//! Default Constructor
		/*! 
			Sets all attributes to 0. _Nval=1 by default so V will be 1 value of 0 /point 
		*/
	Grid();
	
		//! Constructor
		/*! 
			Calls Grid::read_From_Cube(..) to initialize values of grid 
		*/
	Grid(ifstream& nameFile, const PeriodicTable& Table);
	
		/*! 
			Initializes Grid using Domain::Domain(ifstream& nameFile) and read_From_Cube(ifstream& nameFile, int Natoms, const PeriodicTable& Table) and assigns the data from .cube to V
		*/
	void read_From_Cube(ifstream& nameFile, const PeriodicTable& Table);
	
		//! get() function
		/*!
			Returns the data contained in a vector<..<vector
		*/
	vector<vector<vector<vector<double>>>> V() const;
	
		//! get() function
		/*!
			Returns the domain
		*/
	Domain dom() const;
	
		//! get() function
		/*!
			Returns the molecular/atomic structure on the grid
		*/
	Structure str() const;
	
};

#endif //_CDFTT_GRID_H_INCLUDED
