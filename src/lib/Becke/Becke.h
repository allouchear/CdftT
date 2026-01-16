#ifndef CDFTT_BECKE_INCLUDED
#define CDFTT_BECKE_INCLUDED

#include <iostream>
#include <functional>
#include <vector>

#include <Becke/GridPoints.h>
#include <Common/Descriptors.h>
#include <Common/Structure.h>
#include <Orbitals/Orbitals.h>
#include <Utils/FCHK.h>
#include <Utils/LOG.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/WFX.h>

/**
 * @brief Becke class.
 *
 * This class will be used to calculate partial charges.
 */
class Becke
{
    private:
        /** @brief Structure containing the molecule or the system. */
        Structure _molecule;

        /** @brief Orbitals object containing information on the orbitals. */
        Orbitals _orbitals;

        /** @brief Selected grid. */
        GridPoints _grid;

        /** @brief Storage of grid points for each atom: [atom][point][x,y,z,w]. The weights (w) are the product of the weight function (from fuzzy Voronoi decomposition of space) and the volume element. */
        std::vector<std::vector<std::vector<double>>> _grid_points;

        /** @brief Weights associated with each grid point. */
        std::vector<std::vector<double>> _grid_weights;

        /** @brief Volumes associated with each grid point. */
        std::vector<std::vector<double>> _grid_volumes;

        /** @brief Partial charges (per atom). */
        std::vector<double> _partial_charge;

        /** @brief Total energy. */
        double _energy;

        /** @brief Indicates whether multigrid mode is active (true) or not (false). */
        bool _multigrid;


    public:
        /** 
         * @brief Returns the total energy.
         */
        double get_Energy();

        /**
         * @brief Returns the partial charges (per atom).
         */
        std::vector<double> get_Partial_Charge();

        /**
         * @brief Prints partial charges to standard output.
         */
        void printCharges();

        /**
         * @brief Default constructor.
         *
         * This constructor is used to set all of the parameters to 0 or "None" value.
         */
        Becke();

        /**
         * @brief Constructor that sets up the structure.
         *
         * @param s Structure to use for the Becke calculation.
         */
        Becke(const Structure& s);

        /**
         * @brief Constructor that sets up the grid.
         *
         * @param g Grid to use for the Becke calculation.
         */
        Becke(const Grid& g);

        /**
         * @brief Constructor from a .wfx input.
         *
         * This constructor is used to set all parameters using a .wfx file.
         *
         * @param wfx WFX reader reference.
         * @param bin Binomial handler class reference.
         * @param table Periodic table reference.
         */
        Becke(WFX& wfx, Binomial& Bin, const PeriodicTable& Table);

        /**
         * @brief Constructor from a .fchk input.
         *
         * This constructor is used to set all parameters using a .fchk file.
         *
         * @param fchk FCHK reader reference.
         * @param bin Binomial handler class reference.
         * @param table Periodic table reference.
         */
        Becke(FCHK& fchk, Binomial& bin, const PeriodicTable& table);

        /**
         * @brief Constructor from a .molden or .gab input.
         *
         * This constructor is used to set all parameters using a .molden or .gab file.
         *
         * @param moldengab MOLDENGAB reader reference.
         * @param bin Binomial handler class reference.
         * @param table Periodic table reference.
         */
        Becke(MOLDENGAB& moldengab, Binomial& bin, const PeriodicTable& table);

        /**
         * @brief Constructor from a .log input.
         *
         * This constructor is used to set all parameters using a .log file.
         *
         * @param log LOG reader reference.
         * @param bin Binomial handler class reference.
         * @param table Periodic table reference.
         */
        Becke(LOG& log, Binomial& bin, const PeriodicTable& table);

        /**
         * @brief Default destructor.
         *
         * Not used explicitly.
         */
        ~Becke() {}

        /**
         * @brief Returns the Structure (molecule or system).
         *
         * @return Structure The molecule or system.
         */
        Structure get_struct() {return _molecule;}

        /**
         * @brief Returns the number of radial points for a given atomic number.
         *
         * @param Z Atomic number for which the number of radial points is required.
         * @return int Number of radial points.
         */
        int number_of_radial_points(int Z);

        /**
         * @brief Returns a grid for a given lebedev order.
         *
         * @param lebedev_order Lebedev order for angular quadrature.
         * @return GridPoints Angular grid.
         */
        GridPoints select_angular_grid(int lebedev_order);

        /**
         * @brief Constructs and updates _grid_points, _grid_weights and _grid_volumes.
         *
         * Constructs a Becke grid for each atom.
         *
         * @param kmax Indicates how fuzzy the Voronoi polyhedrons should be, with larger kmax values meaning that borders are fuzzier. (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         */
        void multicenter_grids(int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Merges all atomic grids into a single grid and recompute weights.
         *
         * @return Merged grid (columns: x,y,z,w).
         */
        std::vector<std::vector<double>> join_grids();

            //! A normal member taking two arguments and returning a double value.
            /*! \return The value of cutoff profiles. */
        double s(double mu, int k = 3); // ?

        /**
         * @brief Multicenter integration for functions of signature: double(const std::vector<GTF>&, double, double, double).
         *
         * @param f Function to integrate.
         * @param p Vector of GTF passed to f.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return double Value of the integral.
         */
        double multicenter_integration(std::function<double(const std::vector<GTF>&, double,double,double)> f, const std::vector<GTF>& p, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Multicenter integration for functions of signature: double(Orbitals&, int, int, double, double, double).
         *
         * @param f Function to integrate.
         * @param i ?
         * @param j ?
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return double Value of the integral.
         */
        double multicenter_integration(std::function<double(Orbitals&, int, int, double,double,double)>, int, int, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Multicenter integration for functions of signature: double(Orbitals&, int, int, double, double, double, int).
         *
         * @param f Function to integrate.
         * @param i ?
         * @param j ?
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @param spinType Spin type for the integral (default ALPHA).
         * @return double Value of the integral.
         */
        double multicenter_integration(std::function<double(Orbitals&, int, int, double, double, double, SpinType)> f, int i, int j, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5, SpinType spinType = SpinType::ALPHA);

        /**
         * @brief Multicenter integration from a density Grid.
         *
         * Creates a radial Becke grid from a density grid and interpolates the electronic density if the Becke grid points do not match the density grid points.
         *
         * @param g Density grid to use.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return double Value of the integral.
         */
        //! Create Becke grid from density grid
        /*! Creates a radial Becke grid from a density grid. Interpolates the electronic density if the points of Becke grid dont match the density grid*/
        double multicenter_integration(const Grid& g, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief TODO
         */
        double multicenter_integration(std::function<double(Orbitals&, int, int, double, double, double, SpinType, const std::array<double, 3>&, double)> f, int i, int j, SpinType spinType, const std::array<double, 3>& chargePosition, double charge, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Returns the table of integral values for a function of signature double(Orbitals&, double, double, double), evaluated on each grid (so on each atom).
         *
         * @param f Function to evaluate.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return Table of integral values per atom.
         */
        std::vector<double> multicenter_sub_integration(std::function<double(Orbitals&, double, double, double)> f, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
        
        /**
         * @brief ?
         *
         * @param g Density grid to use.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return Table of integral values per atom.
         */
        std::vector<double> multicenter_sub_integration(const Grid& g, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Calculates and updates the partial charge.
         *
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         */
        void partial_charge(int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Calculates and updates the partial charge from a density grid.
         *
         * @param g Density grid to use.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         */
        void partial_charge(const Grid& g, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Calculates and returns the ionic potential energy.
         */
        double ionic_potential(int i, int j, SpinType spinType, const std::array<double, 3>& chargePosition, double charge, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Electronic density value at point (x,y,z).
         *
         * @param Orb Orbitals object used to evaluate the density.
         * @param x X coordinate.
         * @param y Y coordinate.
         * @param z Z coordinate.
         * @return double Electronic density at the given point.
         */
        static double density(Orbitals&, double, double, double);

        /**
         * @brief Returns the overlap integral between two Gaussian-Type Functions (GTFs).
         *
         * @param g1 First GTF.
         * @param g2 Second GTF.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return double Value of the overlap integral.
         */
        double OverlapGTF(const GTF&, const GTF&, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Returns the product of many Gaussian-Type Functions (GTFs) evaluated at (x,y,z).
         *
         * @param p Vector of GTFs.
         * @param x X coordinate.
         * @param y Y coordinate.
         * @param z Z coordinate.
         * @return double Product value at the given point.
         */
        static double prodGTF(const std::vector<GTF>& p, double x, double y, double z);

        /**
         * @brief Returns the value of the overlap integral between two Contracted Gaussian-Type Functions (CGTFs).
         *
         * @param i Index of first CGTF.
         * @param j Index of second CGTF.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return double Value of the overlap integral.
         */
        double OverlapCGTF(int i, int j, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);

        /**
         * @brief Returns the product of two Contracted Gaussian-Type Functions (CGTFs) evaluated at a point (x,y,z).
         *
         * @param Orb Orbitals reference containing the CGTFs.
         * @param i Index of first CGTF.
         * @param j Index of second CGTF.
         * @param x X coordinate.
         * @param y Y coordinate.
         * @param z Z coordinate.
         * @return double Product value of the two CGTFs at (x,y,z).
         */
        static double CGTFstarCGTF(Orbitals& Orb, int i, int j, double x, double y, double z);

        /**
         * @brief Returns the overlap integral between two orbitals of indexes i and j.
         *
         * @param i Index of the first orbital.
         * @param j Index of the second orbital.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @param spinType SpinType for the integral (default ALPHA).
         * @return double Value of the overlap integral.
         */
        double overlap(int i, int j, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5, SpinType spinType = SpinType::ALPHA);

        /**
         * @brief Returns the value of the i-th orbital at point (x,y,z).
         *
         * @param orbitals Orbitals reference.
         * @param i Index of the chosen orbital.
         * @param x X coordinate.
         * @param y Y coordinate.
         * @param z Z coordinate.
         * @param spinType SpinType to consider (default ALPHA).
         * @return double Value of the chosen orbital at the point (x,y,z).
         */
        static double phi(Orbitals& orbitals, int i, double x, double y, double z, SpinType spinType = SpinType::ALPHA);

        /**
         * @brief Returns the product of two orbitals of indexes i and j at a point (x,y,z).
         *
         * @param Orb Orbitals reference.
         * @param i Index of the first orbital.
         * @param j Index of the second orbital.
         * @param x X coordinate.
         * @param y Y coordinate.
         * @param z Z coordinate.
         * @param spinType SpinType to consider (default ALPHA).
         * @return double Product value of the two orbitals at (x,y,z).
         */
        static double phiStarPhi(Orbitals& Orb, int i, int j, double x, double y, double z, SpinType spinType = SpinType::ALPHA);

        /**
         * @brief Returns the product of two orbitals of indexes i and j at a point (x,y,z) multiplied by the electrostatic potential V_ionic created by a point charge.
         *
         * @param[in] orbitals Orbitals reference.
         * @param[in] i Index of the first orbital.
         * @param[in] j Index of the second orbital.
         * @param[in] x X coordinate.
         * @param[in] y Y coordinate.
         * @param[in] z Z coordinate.
         * @param[in] spinType Spin type (ALPHA, BETA, ALPHA_BETA).
         * @param[in] position Position of the charge.
         * @param[in] charge Value of the charge.
         * @return Product value of the two orbitals at (x,y,z) multiplied by the electrostatic potential V_ionic.
         */
        static double phiStarVionicStarPhi(Orbitals& orbitals, int i, int j, double x, double y, double z, SpinType spinType, const std::array<double, 3>& chargePosition, double charge);

        /**
         * @brief Returns the HOMO energy.
         *
         * @return double HOMO energy.
         */
        double eHOMO()
        {
            _orbitals.HOMO();
            return _orbitals.eHOMO();
        }

        /**
         * @brief Returns the LUMO energy.
         *
         * @return double LUMO energy.
         */
        double eLUMO() {_orbitals.LUMO(); return _orbitals.eLUMO();}

        /**
         * @brief Returns partial charges and energy from a Grid.
         *
         * @param g Grid reference.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return Total energy (index 0) and partial charges (starting from index 1).
         */
        std::vector<double> PartialChargeAndEnergy(const Grid& g, int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);


        /**
         * @brief Returns partial charges and energy from a Grid.
         *
         * @param g Grid reference.
         * @param kmax Fuzzyness of the Voronoi polyhedrons (default 3).
         * @param lebedev_order Lebedev order for angular quadrature (default 41).
         * @param radial_grid_factor Radial grid multiplicative factor (default 5).
         * @return Total energy (first column) and partial charges (second column).
         */
        std::vector<std::vector<double>> PartialChargesAndEnergy(int kmax = 3, int lebedev_order = 41, int radial_grid_factor = 5);
};

#endif
