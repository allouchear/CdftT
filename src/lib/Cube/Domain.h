#ifndef CDFTT_DOMAIN_H_INCLUDED
#define CDFTT_DOMAIN_H_INCLUDED

#include <array>
#include <fstream>
#include <vector>

#include <Common/Structure.h>


/**
 * @brief Domain class.
 *
 * This class is used to define a domain containing a structure or a molecule, in a 3D space.
 */
class Domain
{
    private:
        /** @brief Number of values stored per grid point. */
        int _Nval;

        /** @brief Number of grid points along the first (x) axis. */
        int _N1;

        /** @brief Number of grid points along the second (y) axis. */
        int _N2;

        /** @brief Number of grid points along the third (z) axis. */
        int _N3;

        /** @brief Origin coordinates (x,y,z). */
        std::array<double, 3> _origin;

        /** @brief Translation matrix (3x3). Maps grid indices to spatial coordinates. */
        std::array<std::array<double, 3>, 3> _T;

        /** @brief Inverse of the translation matrix. */
        std::array<std::array<double, 3>, 3> _inv_T;

        /** @brief Infinitesimal distance increment along the first (x) axis. */
        double _dx;

        /** @brief Infinitesimal distance increment along the second (y) axis. */
        double _dy;

        /** @brief Infinitesimal distance increment along the third (z) axis. */
        double _dz;

        /** @brief Infinitesimal volume increment (equals to _dx * _dy * _dz). */
        double _dv;

        //----------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        /** @brief Computes the infinitesimal distance and volume increments. */
        void computeInfinitesimalElements();

    
    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//
        
        /**
         * @brief Default Constructor.
         *
         * Sets all attributes to 0 except _Nval which is 1 by default.
         */
        Domain();
    
        /**
         * @brief Constructor from a .cube input.
         *
         * @param nameFile Input file stream opened on a .cube file.
         */
        Domain(std::ifstream& nameFile);

        /**
         * @brief Constructor with explicit grid sizes and maximum coordinates.
         *
         * Sets T to identity matrix.
         *
         * @param Nval Number of values per grid point.
         * @param N1 Number of grid points along the first (x) axis.
         * @param N2 Number of grid points along the second (y) axis.
         * @param N3 Number of grid points along the third (z) axis.
         * @param xmax Maximum coordinate along the first (x) axis.
         * @param ymax Maximum coordinate along the second (y) axis.
         * @param zmax Maximum coordinate along the third (z) axis.
         * @param T Translation matrix.
         */
        Domain(int Nval, int N1, int N2,int N3, double xmax, double ymax, double zmax);

        /**
         * @brief Constructor with explicit grid sizes and origin.
         *
         * Sets T to null matrix. If the origin pointer is NULL, sets origin to null vector.
         *
         * @param nVal Number of values per grid point.
         * @param n1 Number of grid points along the first (x) axis.
         * @param n2 Number of grid points along the second (y) axis.
         * @param n3 Number of grid points along the third (z) axis.
         * @param origin Reference to origin coordinates (3-element array).
         */
        Domain(int nVal, int n1, int n2, int n3, const std::array<double, 3>& origin);


        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//
        
        /**
         * @brief Returns the number of values per grid point.
         */
        int get_Nval() const;

        /**
         * @brief Returns the number of points along the first (x) axis.
         */
        int get_N1() const;
    
        /**
         * @brief Returns the number of points along the second (y) axis.
         */
        int get_N2() const;

        /**
         * @brief Returns the number of points along the third (z) axis.
         */
        int get_N3() const;
    
        /**
         * @brief Returns the coordinates of the origin.
         */
        const std::array<double, 3>& get_origin() const;
    
        /**
         * @brief Returns the translation matrix.
         */
        const std::array<std::array<double, 3>, 3>& get_T() const;
    
        /**
         * @brief Returns a given element of the translation matrix.
         *
         * @param i Row index.
         * @param j Column index.
         * @return Component Tij of the translation matrix.
         */
        double get_Tij(int i, int j) const;

        /**
         * @brief Returns the infinitesimal distance along the first (x) axis.
         */
        double get_dx() const;

        /**
         * @brief Returns the infinitesimal distance along the second (y) axis.
         */
        double get_dy() const;

        /**
         * @brief Returns the infinitesimal distance along the third (z) axis.
         */
        double get_dz() const;
        
        /**
         * @brief Returns the infinitesimal volume.
         */
        double get_dv() const;


        //----------------------------------------------------------------------------------------------------//
        // SETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Sets the number of values per grid point.
         *
         * @param N Number to set.
         */
        void set_Nval(int N);

        /**
         * @brief Sets the number of points along the first (x) axis.
         *
         * @param N Number to set.
         */
        void set_N1(int N);
    
        /**
         * @brief Sets the number of points along the second (y) axis.
         *
         * @param N Number to set.
         */
        void set_N2(int N);
    
        /**
         * @brief Sets the number of points along the third (z) axis.
         *
         * @param N Number to set.
         */
        void set_N3(int N);
        
        /**
         * @brief Sets a given element of the translation matrix.
         *
         * @param value Value to set.
         * @param i Row index.
         * @param j Column index.
         */
        void set_Tij(double value, int i, int j);

        /**
         * @brief Sets all the attributes.
         *
         * @param Nval Number of values per grid point.
         * @param N1 Number of grid points along the first (x) axis.
         * @param N2 Number of grid points along the second (y) axis.
         * @param N3 Number of grid points along the third (z) axis.
         * @param xmax Maximum coordinate along the first (x) axis.
         * @param ymax Maximum coordinate along the second (y) axis.
         * @param zmax Maximum coordinate along the third (z) axis.
         * @param T Translation matrix.
         */
        void set_all(int Nval, int N1, int N2, int N3, double xmax, double ymax, double zmax, const std::array<std::array<double, 3>, 3>& T);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the value along the first (x) axis of the point on indices i, j, k.
         *
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @return First (x) coordinate in space.
         */
        double x(int i, int j, int k) const;

        /**
         * @brief Returns the value along the second (y) axis of the point on indices i, j, k.
         *
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @return Second (y) coordinate in space.
         */
        double y(int i, int j, int k) const;

        /**
         * @brief Returns the value along the third (z) axis of the point on indices i, j, k.
         *
         * @param i Grid index i.
         * @param j Grid index j.
         * @param k Grid index k.
         * @return Third (z) coordinate in space.
         */
        double z(int i, int j, int k) const;

        /**
         * @brief Returns the closest index along the first (x) axis of the grid, for a given point in 3D space.
         *
         * @param x Point coordinate along the first (x) axis.
         * @param y Point coordinate along the second (y) axis.
         * @param z Point coordinate along the third (z) axis.
         * @return Closest grid index along the first (x) axis.
         */
        int i(double x, double y, double z) const;

        /**
         * @brief Returns the closest index along the second (y) axis of the grid, for a given point in 3D space.
         *
         * @param x Point coordinate along the first (x) axis.
         * @param y Point coordinate along the second (y) axis.
         * @param z Point coordinate along the third (z) axis.
         * @return Closest grid index along the second (y) axis.
         */
        int j(double x, double y, double z) const;

        /**
         * @brief Returns the closest index along the third (z) axis of the grid, for a given point in 3D space.
         *
         * @param x Point coordinate along the first (x) axis.
         * @param y Point coordinate along the second (y) axis.
         * @param z Point coordinate along the third (z) axis.
         * @return Closest grid index along the third (z) axis.
         */
        int k(double x, double y, double z) const;

        /**
         * @brief Loads data from a .cube file.
         *
         * Reads from a .cube file and initializes the number of atoms, the geometry, etc.
         *
         * @param nameFile Input file stream opened on a .cube file.
         */
        void readFromCube(std::ifstream& nameFile);

        /**
         * @brief Inverts the translation matrix and stores the result in _invT.
         */
        void inverse_T();

        /**
         * @brief Returns the maximum length betweens the atoms of a given Structure.
         *
         * Note that it is better to use Grid::sizeUpMol().
         *
         * @param S Structure reference containing atoms.
         * @param scale Scale factor (unused in current implementation).
         * @return Maximum inter-atomic distance.
         */
        double sizeUpMol(const Structure& S, double scale);


        //----------------------------------------------------------------------------------------------------//
        // OPERATOR OVERLOADS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Overloads the equality operator for two Domain objects.
         *
         * Two Domains are considered equal if they have the same number of values per grid point,
         * the same grid sizes, and identical translation matrices and origins.
         *
         * @param D Other Domain to compare this Domain with.
         * @return True if the two Domains are equal, false otherwise.
         */
        bool operator==(const Domain &D) const;

        /**
         * @brief Overloads the inequality operator for two Domain objects.
         *
         * Two Domains are considered unequal if they have different numbers of values per grid point,
         * different grid sizes, or different translation matrices or origins.
         *
         * @param D Other Domain to compare this Domain with.
         * @return True if the two Domains are equal, false otherwise.
         */
        bool operator!=(const Domain &D) const;
};

#endif // CDFTT_DOMAIN_H_INCLUDED
