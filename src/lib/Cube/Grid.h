#ifndef _CDFTT_GRID_H_INCLUDED
#define _CDFTT_GRID_H_INCLUDED

#include <fstream>
#include <vector>
#include <Common/Structure.h>
#include <Cube/Domain.h>

using std::ifstream;
using std::ofstream;
using std::vector;

/**
 * @brief Grid class.
 *
 * This class represents a 3D grid for storing scalar or vector fields, typically used for molecular densities and related quantities.
 */
class Grid
{
    private:
        /** @brief Domain object describing the grid geometry and axes. */
        Domain _domain;

        /** @brief Structure object describing the molecule/structure on the grid. */
        Structure _structure;

        /** @brief Values for each grid point. The first three dimensions are spatial dimensions (x,y,z). The fourth dimensions stores the values. */
        vector<vector<vector<vector<double>>>> _values;

        /**
         * @brief Advances the density ascent trajectory to the next point.
         *
         * @param i Grid index along the first (x) axis.
         * @param j Grid index along the second (y) axis.
         * @param k Grid index along the third (z) axis.
         * @param rhocenter Reference to current density value. --> or density of center ??
         * @param trajectory Vector of trajectory points.
         */
        void next_Density(int i, int j, int k, double& rhocenter, vector<vector<int>>& trajectory);

        /**
         * @brief Adds surrounding density points to the trajectory.
         * 
         * @param i Grid index along the first (x) axis.
         * @param j Grid index along the second (y) axis.
         * @param k Grid index along the third (z) axis.
         * @param equals Vector of equal points.
         * @param rhocenter Reference to current density value. --> or density of center ??
         */
        void addSurroundingDensity(int i, int j, int k, vector<vector<int>>& equals, double& rhocenter);

        /**
         * @brief Sets boundary values to the closest interior value. Used by gradient and laplacian.
         *
         * @param nBound Thickness of the boundary layer.
         */
        void reset_Boundary(int nBound);


    public:
        /**
         * @brief Resets all grid values to zero and resizes _V to match the domain size.
         */
        void reset();

        /**
         * @brief Default constructor.
         *
         * Initializes all attributes to default values (calls default constructor on object members, empty vectors.).
         */
        Grid();

        /**
         * @brief Constructor from a domain.
         *
         * Builds a grid with the given domain.
         * 
         * @param d Domain reference.
         */
        Grid(const Domain& d);

        /**
         * @brief Constructor from a .cube file.
         * 
         * @param nameFile Input file stream opened on a .cube file.
         * @param Table PeriodicTable reference.
         */
        Grid(ifstream& nameFile, const PeriodicTable& Table);

        /**
         * @brief Initializes grid from a .cube file.
         *
         * @param nameFile Input file stream opened on a .cube file.
         * @param Table PeriodicTable reference.
         */
        void read_From_Cube(ifstream& nameFile, const PeriodicTable& Table);

        /**
         * @brief Returns the grid data as a 4D vector.
         * 
         * This method is slow since it copies a 4D vector. Avoid it if possible.
         */
        vector<vector<vector<vector<double>>>> get_values() const;

        /**
         * @brief Returns the domain object.
         */
        Domain get_domain() const;

        /**
         * @brief Returns the molecular/atomic structure on the grid.
         */
        Structure get_structure() const;

        /**
         * @brief Sets the domain.
         */
        void set_domain(const Domain& d);

        /**
         * @brief Sets the molecular/atomic structure on the grid.
         */
        void set_structure(const Structure& S);

        /**
         * @brief Sets the grid data.
         */
        void set_values(const vector<vector<vector<vector<double>>>>& U);

        /**
         * @brief Sets the value at grid indices (i,j,k,l).
         * 
         * @param rho Value to set
         * @param i Grid index i
         * @param j Grid index j
         * @param k Grid index k
         * @param l Value index l
         */
        void set_Vijkl(double rho, int i, int j, int k, int l);

        /**
         * @brief Overloads the addition operator for two grid objects.
         * 
         * @param g Grid to add.
         * @return Sum grid.
         */
        Grid operator+(const Grid& g); // ?? Should be a const method. Also: outside of class?

        /**
         * @brief Adds another grid to this one.
         * 
         * @param g Grid to add.
         */
        Grid add(const Grid& g); // ?? Should return Grid& ? Refactor: overload operator+= ?

        /**
         * @brief Overloads the multiplication operator to multiply two grids pointwise.
         * 
         * @param g Grid to multiply with.
         * @return Product grid.
         */
        Grid operator*(const Grid& g);

        /**
         * @brief Overloads the substraction operator to subtract two grids pointwise.
         
         * @param g Grid to subtract.
         * @return Difference grid.
         */
        Grid operator-(const Grid& g);

        /**
         * @brief Returns the sum of all grid values over all points.
         */
        double sum();

        /**
         * @brief Returns a Coulomb potential grid generated by a point charge.
         * 
         * @param q Charge value.
         * @param R Charge position.
         * @return Coulomb grid.
         */
        Grid coulomb_Grid(double q, vector<double> R);

        /**
         * @brief Integrates the values over the domain.
         * 
         * @return Integral value.
         */
        double integrate_Over_Dom();

        /**
         * @brief Computes the coefficients of the Laplacian operator in finite difference.
         *
         * @param nBound Boundary thickness.
         * @param fcx Laplacian component along the first (x) axis.
         * @param fcy Laplacian component along the second (y) axis.
         * @param fcz Laplacian component along the third (z) axis.
         * @param cc Output central coefficient.
         */
        void coefs_Laplacian(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz, double& cc) const;

        /**
         * @brief Returns a grid which values contain the Laplacian of this grid.
         * 
         * Note: the outer layers of thickness nBound will be zeroes as required by finite difference.
         * 
         * @param nBound Boundary thickness.
         * @return Laplacian grid.
         */
        Grid laplacian(int nBound) const;

        /**
         * @brief Computes the coefficients of the gradient operator in finite difference.
         * @param nBound Boundary thickness.
         * @param fcx Gradient component along the first (x) axis.
         * @param fcy Gradient component along the second (y) axis.
         * @param fcz Gradient component along the third (z) axis.
         */
        void coefs_Gradient(int nBound, vector<double>& fcx, vector<double>& fcy, vector<double>& fcz) const;

        /**
         * @brief Returns a grid which values contain the gradient of this grid.
         * 
         * Note: the outer layers of thickness nBound will be zeroes as required by finite difference.
         * Note 2: number of values per point will change to 4.
         * 
         * @param nBound Boundary thickness.
         * @return Gradient grid.
         */
        Grid gradient(int nBound) const;

        /**
         * @brief Returns the same structure with a finer grid.
         * 
         * Intermediate values are calculated by cubic interpolation.
         * 
         * @return Finer grid.
         */
        Grid finer_Grid();

        /**
         * @brief Returns the same structure with a coarser grid.
         * 
         * Values are averaged across the grid.
         * 
         * @return Coarser grid.
         */
        Grid coarser_Grid();

        /**
         * @brief Saves the grid to a .cube file.
         * 
         * @param name Output file stream opened on a .cube file.
         */
        void save(ofstream& name);

        /**
         * @brief Returns the value at grid indices (i,j,k,0).
         * 
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @return Value at (i,j,k,0).
         */
        double value(int i, int j, int k) const;

        /**
         * @brief Returns the value at grid indices (i,j,k,l).
         * 
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @param l Value index l.
         * @return Value at (i,j,k,l).
         */
        double value(int i, int j, int k, int l) const; //refactor : une seule fonction avec valeur l par defaut à 0 ?

        /**
         * @brief Advances to the next point in a trajectory.
         * 
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @param Current Reference to current value.
         * @param traj Trajectory vector.
         */
        void next(int i, int j, int k, double& Current, vector<vector<int>>& traj);

        /**
         * @brief Computes atom-attractor differences for a set of attractor points.
         * 
         * @param attract Vector of attractor points.
         * @return Vector of differences for each atom.
         */
        vector<double> atom_attract_diff(const vector<vector<int>>& attract);

        /**
         * @brief Adds surrounding equal points to the trajectory.
         * 
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @param equals Vector of equal points.
         * @param current Reference to current value.
         */
        void addSurroundingEqualPoints(int i,int j,int k, vector<vector<int>>& equals, double& current);

        /**
         * @brief Computes Atoms in Molecules (AIM) regions on the grid.
         *
         * Ref: A grid-based Bader analysis algorithm without lattice bias.
         * W. Tang, E. Sanville and G. Henkelman
         * J. Phys.: Condens. Matter (2009), 21 084204. DOI 10.1088/0953-8984/21/8/084204.
         *
         * @param nBound Boundary thickness.
         * @return AIM grid.
         */
        Grid aim_On_Grid(int nBound);

        /**
         * @brief Computes Atoms in Molecules (AIM) regions based on density.
         *
         * Ref: A grid-based Bader analysis algorithm without lattice bias.
         * W. Tang, E. Sanville and G. Henkelman
         * J. Phys.: Condens. Matter (2009), 21 084204. DOI 10.1088/0953-8984/21/8/084204.
         * 
         * @return Density-based AIM grid.
         */
        Grid aim_On_Grid_Density();

        /**
         * @brief Measures the size of the molecular structure and returns the largest sizes in each direction.
         * 
         * @param scale Scale factor.
         * @return Vector of largest sizes in each direction (x,y,z).
         */
        vector<double> sizeUpMol(double scale);

        /**
         * @brief Returns the electronic density at a point in space, weighted by the distance to grid points.
         * 
         * @param x First (x) coordinate.
         * @param y Second (y) coordinate.
         * @param z Third (z) coordinate.
         * @return Weighted electronic density.
         */
        double value(double x, double y, double z) const;
};


#endif //_CDFTT_GRID_H_INCLUDED
