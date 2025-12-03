#ifndef CDFTT_DOMAIN_H_INCLUDED
#define CDFTT_DOMAIN_H_INCLUDED

#include <fstream>
#include <vector>

#include <Common/Structure.h>

using std::ifstream;
using std::vector;

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
        double _O[3];

        /** @brief Translation matrix (3x3). Maps grid indices to spatial coordinates. */
        vector<vector<double>> _T;

        /** @brief Inverse of the translation matrix. */
        vector<vector<double>> _inv_T;

        /** @brief Infinitesimal distance increment along the first (x) axis. */
        double _dx;

        /** @brief Infinitesimal distance increment along the second (y) axis. */
        double _dy;

        /** @brief Infinitesimal distance increment along the third (z) axis. */
        double _dz;

        /** @brief Infinitesimal volume increment (equals to _dx * _dy * _dz). */
        double _dv;

    
    public:
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
        Domain(ifstream& nameFile);

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
        Domain(int Nval, int N1, int N2,int N3, double xmax, double ymax, double zmax, vector<vector<double>> T);

        /**
         * @brief Constructor with explicit grid sizes and origin.
         *
         * Sets T to null matrix. If the origin pointer is NULL, sets origin to null vector.
         *
         * @param i Number of values per grid point.
         * @param n Number of grid points along the first (x) axis.
         * @param m Number of grid points along the second (y) axis.
         * @param l Number of grid points along the third (z) axis.
         * @param O Pointer to origin coordinates (3-element array) or NULL.
         */
        Domain(int i, int n, int m, int l, double* O);
    
        /**
         * @brief Loads data from a .cube file.
         *
         * Reads from a .cube file and initializes the number of atoms, the geometry, etc.
         *
         * @param nameFile Input file stream opened on a .cube file.
         */
        void read_From_Cube(ifstream& nameFile);

        /**
         * @brief Returns the number of values per grid point.
         */
        int Nval() const;
        
        /**
         * @brief Sets the number of values per grid point.
         *
         * @param N Number to set.
         */
        void set_Nval(int N);

        /**
         * @brief Returns the number of points along the first (x) axis.
         */
        int get_N1() const;

        /**
         * @brief Sets the number of points along the first (x) axis.
         *
         * @param N Number to set.
         */
        void set_N1(int N);
    
        /**
         * @brief Returns the number of points along the second (y) axis.
         */
        int get_N2() const;
    
        /**
         * @brief Sets the number of points along the second (y) axis.
         *
         * @param N Number to set.
         */
        void set_N2(int N);

        /**
         * @brief Returns the number of points along the third (z) axis.
         */
        int get_N3() const;
    
        /**
         * @brief Sets the number of points along the third (z) axis.
         *
         * @param N Number to set.
         */
        void set_N3(int N);
    
        /**
         * @brief Returns the coordinates of the origin.
         */
        double* O() const;
    
        /**
         * @brief Returns the translation matrix.
         */
        vector<vector<double>> T() const;
    
        /**
         * @brief Returns a given element of the translation matrix.
         *
         * @param i Row index.
         * @param j Column index.
         * @return Component Tij of the translation matrix.
         */
        double Tij(int i, int j) const;
        
        /**
         * @brief Sets a given element of the translation matrix.
         *
         * @param v Value to set.
         * @param i Row index.
         * @param j Column index.
         */
        void set_T(double v, int i, int j);

        /**
         * @brief Returns the infinitesimal distance along the first (x) axis.
         */
        double dx() const;

        /**
         * @brief Returns the infinitesimal distance along the second (y) axis.
         */
        double dy() const;

        /**
         * @brief Returns the infinitesimal distance along the third (z) axis.
         */
        double dz() const;
        
        /**
         * @brief Returns the infinitesimal volume.
         */
        double get_dv() const;

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
        void set_all(int Nval, int N1, int N2, int N3, double xmax, double ymax, double zmax, vector<vector<double>> T);
};

#endif // CDFTT_DOMAIN_H_INCLUDED
