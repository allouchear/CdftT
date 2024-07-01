#ifndef _CDFTT_GRIDCP_H_INCLUDED
#define _CDFTT_GRIDCP_H_INCLUDED

using namespace std;
#include <numeric/Domain.h>
#include <common/Structure.h>
#include <functional>
#include <vector>
#include <numeric/Grid.h>

	//! Critical point
	/*! structure containing data on attractors in the AIM framework. */
struct CriticalPoint
{
	int index[3];
	int rank;
	int signature;
	double lambda[3];
	double integral;
	double nuclearCharge;
	int numVolume;
	int numCenter;
	double volume;
	double value;
};

	//! Critical point grid
	/*! Grid class used to find the domains associated to atoms in molecules using multiple methods specified in W. Tang et al J. Phys. Condens. Matter 21 (2009) 084204, DOI 10.1088/0953-8984/21/8/084204*/
class GridCP
{
	private:
		Domain _domain;
		Structure _str;
		vector<vector<vector<vector<double>>>> _V;
		vector< vector <  vector<int> > > _volumeNumberOfPoints;
		vector< vector <  vector<int> > > _known;
		vector< CriticalPoint > _criticalPoints;
		double _integral;
		double _nuclearCharge;
		
			//! New critical point
			/*! Creates a new critical point*/
		CriticalPoint  newCriticalPoint(int i, int j, int k, int numV);
		
			//! reset known
			/*! resets the vector containing known points to 0s.*/
		void resetKnown();
		
			//! Grid initialisation
			/* Initialises gridCP using a Grid. Copies domain, structure, and data into a GridCP and returns it. if near grid method is chosen, the method calculates the gradient of the grid and assigns those values to _V. */
		void initGridCP(const Grid& grid, bool ongrid);
		
			//! Set neighbouring points
			/*! Sets neighbouring points' _known values to kn. Takes the coordinates of the initial point as parameter*/
		int setArroundTo(int current[], int kn);
		
			//! Flag maximum
			/*! returns True if the point has the maximum value between it's neighbours*/
		bool isMax(int current[]);
			
			//! Chooses next point
			/*! Finds the neighbouring point ON the grid which maximizes the gradient density projection and equals next[] to it*/
		bool nextPointOnGrid(int current[], int next[]);
		
			//! Chooses next point
			/*! Same as On Grid but calculates the correction vector required t find the true assent trajectory for near grid method*/
		bool nextPoint(double deltaR[], int current[], int next[]);
		
			//! Add surrounding equal gradient points
			/*! Adds the surrounding points of equal gradient projection to the list of visited points*/
		bool addSurroundingEqualPoints(int current[], vector<vector<int>>& listOfVisitedPoints);
		
			//! Compute assent trajectory 
			/*! Uses previous methods to find the trajectory of points on grid or not leading to a maximum. Returns the list of visited points*/
		vector<vector<int>> assentTrajectory(int current[],bool ongrid, bool refine);
		
			//!Assign grid ints
			/*! Assigns the points of the grid to bader volumes with onGrid or near grid*/
		void assignPointsByGradient(bool ongrid);
		
			//! Calculate NumCenters
			/*! Finds the closest atom of the structure to each critical point and assigns them to it*/
		void computeNumCenters();
			
			//! reset function
			/*! Resets all attributes of GridCP to 0*/
		void reset();
		
			//! Domain Flag
			/*! Checks the values of _domain.N1, .N2, .N3 and .Nval and returns false if any are smaller than 2( 1 for Nval)*/
		bool okDomain();
		
			//! Init 3D vector
			/*! Returns a vector of size N1xN2xN3 of zeroes*/
		vector< vector < vector<int>  > > get3DIntVector();
		
			//! Edge of grid flag
			/*! Checks neighbouring points and compares the number of points in their volumes to current and returns true if they are different values*/
		bool isVolumeEdge(int current[]);
		
			//! Compute Volumes
			/*! Computes the volumes of each critical point*/
		void computeVolumes();
		
			//! Remove attractor
			/*! Removes attractors with density smaller than TOL*/
		void removeBasins0();
		
			//!Remove attractor
			/*!Removes attractors whose domain is smaller than Npoints/1000*/
		void removeNonSignificantBasins();
		
			//! Refine edge points
			/*! refines points adjacent to bader surface*/
		int refineEdgeNearGrad();

			//!Assign grid points
			/*! Assigns the points of the grid to VDD volumes.  Ref : Fonseca Guerra et al. https://doi.org/10.1002/jcc.10351 */
		void buildVDD();


		bool addSurroundingSignPoints(int current[], vector<vector<int>>& listOfVisitedPoints, double cutoff);

	public:
		
			//!Default constructor
			/*! calls rest() to set all attributes to 0*/
		GridCP();
			
			//! Build Basins
			/*! calls methods to build basins using Grid and method=0,1,2,3 ( on grid, near grid, near+refine, VDD resp.)*/
		void buildBasins(const Grid&, int method);

			//! Compute Integrals
			/*! Integrates over attractor indices*/
		void computeIntegrals(const Grid& grid);
		
			//! Compute AIM
			/*! Integrates over attractor indices and calculates AIM charges for each atom*/
		vector<double> computeAIMCharges(const Grid&);
		
			//! print function
			/*! prints coordinates, density and nuclear charge of attractors*/
		void printCriticalPoints();

		Structure str() const;

			//!Assign grid points
			/*! Assigns the points based on the sign of function.  Refs :  V. Tognetti : https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.23840 &  T. le Bahers : https://pubs.acs.org/doi/abs/10.1021/ct200308m */
		void build2BasinSign(const Grid& grid);
			//!Assign grid points
			/*! Assigns the points based on the sign of function.  Refs :  V. Tognetti : https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.23840 &  T. le Bahers : https://pubs.acs.org/doi/abs/10.1021/ct200308m */
		void buildBasinsBySign(const Grid& grid, double cutoff);
};

#endif //_CDFTT_GRIDCP_H_INCLUDED
