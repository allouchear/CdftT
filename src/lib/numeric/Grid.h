#ifndef _CDFTT_GRID_H_INCLUDED
#define _CDFTT_GRID_H_INCLUDED

using namespace std;
#include <numeric/Domain.h>
#include <common/Structure.h>
#include <functional>


class Grid
{
	Domain _dom;
	Structure _str;
	vector<vector<vector<vector<double>>>> _V;
	
	void next_Density(int i, int j, int k, double& rhocenter, vector<vector<int>>& trajectory);
	
	void addSurroundingDensity(int i,int j,int k, vector<vector<int>>& equals, double& rhocenter);
	
	//! sets boundary values
	/*! sets the boundary values to closest value of the interior. Used by gradient and laplacian*/
	void reset_Boundary(int nBound);
	
	public:
		//! reset
		/*! resets all values _V to zero and resizes to _dom size*/
	void reset();
	
		//! Default Constructor
		/*! 
			Sets all attributes to 0. _Nval=1 by default so V will be 1 value of 0 /point 
		*/
	Grid();
	
		//! constructor from domain d
		/*! builds a grid with domain d*/
	Grid(const Domain& d);
	
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
			Returns the data contained in a vector<..<vector. We recommend avoiding this method as much as possible as it is very slow.
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
	
		//! Set function
		/*! Sets _dom*/
	void set_dom(const Domain& d);
	
		//! Set function
		/*! Sets _str*/
	void set_str(const Structure& S);
	
		//! Set function
		/*! Sets _V*/
	void set_V(const vector<vector<vector<vector<double>>>>& U);
	
		//! operator +
		/*!
			Operator + overload using Structure::operator+(g)
		*/
	Grid operator+(const Grid& g);
	
		//! addition
		/*! adds grid g to this using Structure::add(..)*/
	Grid add(const Grid& g);
	
	
		//! operator *
		/*!
			returns the product _V and g._V point by point. Structure result is that of the LHS.
		*/
	Grid operator*(const Grid& g);
	
		//! operator -
		/*!
			returns the difference of _V and g._V point by point. Struccture result is that of the LHS.
		*/
	Grid operator-(const Grid& g);
		
		//! Sum
		/*! sums the values of _V over the domain*/
	double sum();
	
		//!Constructor
		/*! Creates a Grid with the default structure, the domain of g, and _V=f(x,y,z). fpar contains information on the type of function used*/
	Grid coulomb_Grid(double q, vector<double> R);
	
		//! Integration
		/*! Interates _V over the domain sum()*dV*/
	double integrate_Over_Dom();
	
		//! Laplacian coefficients
		/*! Gets the coefficients of the laplcian operator in finite difference*/
	void coefs_Laplacian(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc) const;

		//! Laplacian grid
		/*! returns a grid g where g._V are the laplacian of _V.*/
		//* Note: the outer layers of thickness nBound of the grid will be zeroes as required by finite diff*/
	Grid laplacian(int nBound) const;
	
		//! Laplacian coefficients
		/*! Gets the coefficients of the laplcian operator in finite difference*/
	//void coefs_Gradient(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc) const; ALLOUCHE NON
	void coefs_Gradient(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz) const;

		//! Gradient grid
		/*! returns a grid g where g._V are the gradient of _V.*/
		/*! Note: the outer layers of thickness nBound of the grid will be zeroes as required by finite diff*/
		/*! Note 2: number of values per point will change to 4*/
	Grid gradient(int nBound) const;
	
		//! fine grid
		/*! returns the same structure with a finer grid. Intermediate values calculated by cubic interpolation*/
	Grid finer_Grid();
	
		//!coarse grid
		/*! returns the same structure with a coarser grid. values are meaned across the grid.*/
	Grid coarser_Grid();
	
		//!save
		/*!save grid onto .cube file*/
	void save(ofstream& name);

	double value(int i, int j, int k) const;
	double value(int i, int j, int k, int l) const;
	
		
	void next(int i, int j, int k, double& Current, vector<vector<int>>& traj);
	
	
	vector<double> atom_attract_diff(const vector<vector<int>>& attract);
	
	
	void addSurroundingEqualPoints(int i,int j,int k, vector<vector<int>>& equals, double& current);
	
		//include reference Tang
	
	Grid aim_On_Grid(int nBound);
	Grid aim_On_Grid_Density();
	
	double value(double x, double y, double z) const;
};


#endif //_CDFTT_GRID_H_INCLUDED
