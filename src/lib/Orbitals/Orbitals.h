#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include <Basis/CGTF.h>
#include <Common/Descriptors.h>
#include <Common/PeriodicTable.h>
#include <Common/Structure.h>
#include <Cube/Grid.h>
#include <Utils/Enums.hpp>
#include <Utils/FCHK.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/WFX.h>



    //! An Orbitals class.
    /*! This class will be used to calculate descriptors. */
class Orbitals
{
    private:
        /** @brief Vector of CGTFs that compose the basis set. */
        std::vector<CGTF> _vcgtf;

        /** @brief Vector of CGTFs that compose the basis set without normalization. */
        std::vector<CGTF> _vcgtfUnnormalized;

        /** @brief Coefficients of the CGTF for the molecular orbitals. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals. The second dimension gives the coefficient of the i-th orbital. */
        std::vector<std::vector<std::vector<double>>> _coefficients;

        /** @brief Number of atomic orbitals. */
        int _numberOfAo;

        /** @brief Number of primitive functions. */
        int _numberOfGtf;

        /** @brief Number of molecular orbitals. */
        int _numberOfMo;

        /** @brief Number of alpha electrons. */
        int _numberOfAlphaElectrons;

        /** @brief Number of beta electrons. */
        int _numberOfBetaElectrons;

        /** @brief Number of atoms. */
        int _numberOfAtoms;

        /** @brief Table of primitive centers (atom indices) for each CGTF. */
        std::vector<int> _primitiveCenters;

        /** @brief Structure associated with the orbitals. */
        Structure _struct;

        /** @brief Table of atomic numbers for each atom. */
        std::vector<int> _atomicNumbers;

        /** @brief Table of atomic symbols for each atom. */
        std::vector<std::string> _symbol;

        /** @brief Table of molecular orbital energies. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals. The second dimension gives the energy of the i-th orbital. */
        std::vector<std::vector<double>> _orbitalEnergy;

        /** @brief Table of all fukui function values. The first index is for f-, the second for f+. The second dimension gives the fukui function value for the i-th atom. */
        std::vector<std::vector<double>> _all_f;

        /** @brief Orbital indexes (0-based) of the HOMO (index 0) and LUMO (index 1). */
        std::array<int, 2> _homoLumoIndexes;
        
        /** @brief Table of occupation numbers for each orbital. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals. The second dimension gives the occupation number of the i-th orbital. */
        std::vector<std::vector<double>> _occupationNumber;

        /** @brief Flag indicating how the spins are treated. If true, only alpha spin is considered (and is used to represent beta spin as well, i.e. the occupation number can be 2). Otherwise, alpha and beta spins are treated separately. */
        bool _alphaAndBeta;

        /** @brief Helper class to compute binomial coefficients. */
        Binomial _bino;

        /** @brief Descriptors for the orbitals. See @ref Descriptors for more information about the descriptors. */
        Descriptors _descriptors;

        double _energy;
        std::vector<double> _coordinates;
        bool _mixte;

        // debug
        std::vector<std::vector<double>> __debug_AOMatrix;
        std::vector<std::vector<std::vector<double>>> __debug_MOMatrix;
        double __debug_totalSumAO = 0.0;
        std::vector<double> __debug_totalSumMO = std::vector<double>(2, 0.0);


        //----------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the electronic density for a given spin-orbital with CGTFs already evaluated at a given point.
         *
         * @param[in] orbitalIndex The index of the orbital (0-based) for which to compute the electronic density.
         * @param[in] spinType The spin type of the orbital for which to compute the density.
         * @param[in] evaluatedCgtfs The evaluated CGTFs at the given point.
         *
         * @return double Electronic density for the specified spin-orbital at the given point.
         */
        double density(int orbitalIndex, SpinType spinType, const std::vector<double>& evaluatedCgtfs) const;

        /**
         * @brief Evaluates the CGTFs at the given point (x, y, z) and stores the results in the provided vector.
         * 
         * @param[in, out] evaluatedCgtfs A vector to store the evaluated CGTF values at the given point.
         * @param[in] coordinates An array containing the (x, y, z) coordinates at which to evaluate the CGTFs.
         */
        void evaluateCgtfsAtPoint(std::vector<double>& evaluatedCgtfs, const std::array<double, 3>& coordinates) const;

        /**
         * @brief Computes the square of the molecular orbital value (phi^2) for a given spin-orbital with CGTFs already evaluated at a given point.
         * 
         * @param[in] orbitalIndex The index of the orbital (0-based) for which to compute phi^2.
         * @param[in] spinType The spin type of the orbital for which to compute phi^2.
         * @param[in] evaluatedCgtfs The evaluated CGTFs at a given point.
         * 
         * @return double The value of phi^2 for the specified spin-orbital at a given point.
         */
        double phiSquared(int orbitalIndex, SpinType spinType, const std::vector<double>& evaluatedCgtfs) const;


    public :
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        //! A default constructor.
        /*! This constructor is used to set all of the parameters for Orbitals on 0 or "None" value. */
        Orbitals();

        //! A real constructor.
        /*! This constructor is used to add all of the parameters for Orbitals with the data in .wfx file. */
        Orbitals(WFX& wfxParser, const Binomial& bino, const PeriodicTable& periodicTable);

            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .fchk file. */
        Orbitals(FCHK& fchkParser, const Binomial& bino, const PeriodicTable& periodicTable);

            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .molden or .gab file. */
        Orbitals(MOLDENGAB& moldengabParser, const Binomial& bino, const PeriodicTable& periodicTable);
            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .log file. */
        Orbitals(LOG& logParser, const Binomial& bino, const PeriodicTable& periodicTable);


        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the table of all f values.
         */
        const std::vector<std::vector<double>>& get_all_f() const;

        /**
         * @brief Returns the table of CGTF which compose the Orbitals.
         */
        std::vector<CGTF> get_vcgtf() const;

        /**
         * @brief Returns the table of coefficients of all spin-orbitals which compose the Orbitals.
         */
        const std::vector<std::vector<std::vector<double>>>& get_coefficients() const;

        /**
         * @brief Returns the number of atomic orbitals.
         */
        int get_numberOfAo() const;

        /**
         * @brief Returns the number of molecular orbitals.
         */
        int get_numberOfMo() const;

        /**
         * @brief Returns the number of atoms.
         */
        int get_numberOfAtoms() const;

        /**
         * @brief Returns the table of primitive centers.
         */
        const std::vector<int>& get_primitiveCenters() const;

        /**
         * @brief Returns the Structure associated with the Orbitals.
         */
        const Structure& get_struct() const;

        /**
         * @brief Returns the table of atoms' symbols.
         */
        const std::vector<std::string>& get_symbol() const;

        /**
         * @brief Returns the table of molecular orbital's energy.
         */
        const std::vector<std::vector<double>>& get_orbitalEnergy() const;

        /**
         * @brief Returns the table of occupation numbers.
         */
        const std::vector<std::vector<double>>& get_occupationNumber() const;

        /**
         * @brief Returns false if alpha and beta spins are treated separately.
         */
        bool get_alphaAndBeta() const;

        /**
         * @brief Returns the Descriptors object associated with the Orbitals.
         */
        const Descriptors& get_descriptors() const;

        /**
         * @brief Returns the total energy associated with the Orbitals (ground state energy).
         */
        const double get_energy() const;


        //----------------------------------------------------------------------------------------------------//
        // SETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Sets the coefficients for all orbitals.
         * 
         * @param[in] coefficients The 3D vector of doubles where the coefficient values are stored. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals.
         */
        void set_coefficients(const std::vector<std::vector<std::vector<double>>>& coefficients);

        /**
         * @brief Sets the total energy associated with the Orbitals (ground state energy).
         * 
         * @param[in] energy The energy value to set, in Hartree.
         */
        void set_energy(double energy);

        /**
         * @brief Sets the orbital indexes (0-based) of the HOMO and LUMO.
         * 
         * @param[in] homoIndex The index of the HOMO.
         * @param[in] lumoIndex The index of the LUMO.
         */
        void set_homoLumoIndexes(int homoIndex, int lumoIndex);

        /**
         * @brief Sets the number of alpha electrons.
         *
         * @param[in] numberOfAlphaElectrons The number of alpha electrons.
         */
        void set_numberOfAlphaElectrons(int numberOfAlphaElectrons);

        /**
         * @brief Sets the number of beta electrons.
         *
         * @param[in] numberOfBetaElectrons The number of beta electrons.
         */
        void set_numberOfBetaElectrons(int numberOfBetaElectrons);

        /**
         * @brief Sets the number of atomic orbitals.
         *
         * @param[in] numberOfAo The number of atomic orbitals.
         */
        void set_numberOfAo(int numberOfAo);

        /**
         * @brief Sets the number of molecular orbitals.
         * 
         * @param[in] numberOfMo The number of molecular orbitals.
         */
        void set_numberOfMo(int numberOfMo);

        /**
         * @brief Sets the occupation numbers.
         * 
         * @param[in] occupationNumber The 2D vector of doubles where the occupation numbers are stored. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals.
         */
        void set_occupationNumber(const std::vector<std::vector<double>>& occupationNumber);

        /**
         * @brief Sets the orbital energies.
         * 
         * @param[in] orbitalEnergy The 2D vector of doubles where the energy values are stored. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals. The second dimension gives the energy of the i-th orbital.
         */
        void set_orbitalEnergy(const std::vector<std::vector<double>>& orbitalEnergy);

        /**
         * @brief Sets the table of CGTF which compose the Orbitals.
         * 
         * @param[in] vcgtf The vector of CGTF to set for the Orbitals.
         */
        void set_vcgtf(const std::vector<CGTF>& vcgtf);

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Evaluates the molecular orbitals at the given point (x, y, z) and stores the results in the provided 2D vector.
         * 
         * @param[out] evaluatedMOs A 2D vector to store the evaluated molecular orbital values at the given point. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals. The second dimension gives the value of the i-th orbital at the given point.
         * @param[in] coordinates An array containing the (x, y, z) coordinates at which to evaluate the molecular orbitals.
         */
        void evaluateAtPoint(std::vector<std::vector<double>>& evaluatedMOs, const std::array<double, 3>& coordinates) const;

        /**
         * @brief Returns the coefficients of the HOMO. The first dimension is for the spin and the second dimension is for the orbital (so it has a length of 1) to stay consistent with the _coefficients field.
         */
        std::vector<std::vector<std::vector<double>>> getHomoCoefficients();

        /**
         * @brief Returns the HOMO's energy.
         */
        std::vector<double> getHomoEnergy();

        /**
         * @brief Returns the index (0-based) of the HOMO.
         */
        int getHomoIndex();

        /**
         * @brief Returns the coefficients of the LUMO. The first dimension is for the spin and the second dimension is for the orbital (so it has a length of 1) to stay consistent with the _coefficients field.
         */
        std::vector<std::vector<std::vector<double>>> getLumoCoefficients();

        /**
         * @brief Returns the LUMO's energy.
         */
        std::vector<double> getLumoEnergy();

        /**
         * @brief Returns the index (0-based) of the LUMO.
         */
        int getLumoIndex();

        /**
         * @brief Returns the numbers of the occupied orbitals for alpha and beta spins.
         *
         * @return A 2D vector of integers corresponding to the number of the occupied orbitals. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals.
         */
        std::vector<std::vector<int>> getOccupiedOrbitalNumbers() const;

        /**
         * @brief Returns the numbers of the occupied and virtual orbitals for alpha and beta spins.
         * 
         * This method is slightly more efficient than calling getOccupiedOrbitalNumbers() and getVirtualOrbitalNumbers() separately, as it only loops once over the occupation numbers.
         * 
         * @param[out] occupiedOrbitalNumbers The 2D vector of integers where the number of the occupied orbitals are stored. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals.
         * @param[out] virtualOrbitalNumbers The 2D vector of integers where the number of the virtual orbitals are stored. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals.
         */
        void getOccupiedAndVirtualOrbitalNumbers(std::vector<std::vector<int>>& occupiedOrbitalNumbers, std::vector<std::vector<int>>& virtualOrbitalNumbers) const;

        /**
         * @brief Returns the i-th primitive center.
         */
        int getPrimitiveCenter(int i) const;

        /**
         * @brief Returns the numbers of the virtual orbitals for alpha and beta spins.
         * 
         * @return A 2D vector of integers corresponding to the number of the virtual orbitals. The first index is for alpha spin orbitals and the second index is for the beta spin orbitals.
         */
        std::vector<std::vector<int>> getVirtualOrbitalNumbers() const;
        
            //! A normal member taking no arguments and returning a void value.
            /*! Normalise all the CGTF which compose the Orbitals. */
        void NormaliseAllBasis();

            //! A normal member taking three arguments and returning a double value.
            /*! \return The value of ERI ???. */
        double ERIorbitals(Orbitals& q, Orbitals& r, Orbitals& s);

            //! A normal member taking three arguments and returning a double value.
            /*! \return The overlap between two orbitals i and j. */
        double overlap(const int i, const int j, const SpinType spinType);

            //! A normal member taking three arguments and returning a void value.
            /*! Print the value of an overlap between two orbitals i and j. */
        void printOverlap(const int i, const int j, const SpinType spinType);

        //! A normal member taking four arguments and returning a double value.
        /*! \return The overlap between three orbitals i, j, and k. */
        double Overlap3Orbitals(int i, int j, int k, int alpha=0);

            //! A normal member taking five arguments and returning a double value.
            /*! \return The overlap between four orbitals i, j, k, and l. */
        double Overlap4Orbitals(int i, int j, int k, int l, int alpha=0);

        /**
         * @brief Calculates the kinetic energy integral < psi | -1/2 nabla^2 | psi >.
         */
        double kinetic();

        /**
         * @brief Calculates the ion-electron integral for a given ion: < phi_i | V_ion | phi_j >.
         *
         * @param[in] chargePosition The position of the ion.
         * @param[in] charge The charge of the ion.
         * @return The lower triangular matrix < phi_i | V_ion | phi_j > (i >= j). The first index of this 3D vector gives the alpha spin matrix, the second gives the beta spin one.
         */
        std::vector<std::vector<std::vector<double>>> getIonicPotentialMatrix(const std::array<double, 3>& chargePosition, double charge, bool debug = false, bool printAOMatrix = false, bool printMOMatrix = false);

        /**
         * @brief Calculates the ion-electron integral for a given ion using a pseudo CGTF built from a constant unit GTF (exponant = 0, coefficient = 1).
         * 
         * @param[in] chargePosition The position of the ion.
         * @param[in] charge The charge of the ion.
         * @return The vector of < phi_i | V_ion | pseudoPhi >. The first index of this 2D vector gives the alpha spin values, the second gives the beta spin ones.
         */
        std::vector<std::vector<double>> getIonicPotentialVector_unitPseudoCgtf(const std::array<double, 3>& chargePosition, double charge, bool debug = false, bool printAOVector = false, bool printMOVector = false);

        /**
         * @brief Calculates the 3D matrix of triple orbital integrals. The first index corresponds to alpha spin orbitals and the second index corresponds to beta spin orbitals.
         * 
         * @param[in] showProgress If true, shows progress of the calculation.
         */
        std::vector<std::vector<std::vector<std::vector<double>>>> getTripleOrbitalIntegralMatrix(bool showProgress = false);

        //! A normal member taking no arguments and returning a double value.
        /*! \return The value of the integral of Orbitals * Orbitals. */
        double OrbstarOrb();

            //! A normal member taking three arguments and returning a double value.
            /*! \return The value of the integral of Orbitals * Orbitals(ix, iy, iz) * Orbitals. */
        double OrbxyzOrb(int ix, int iy, int iz);

            //! A normal member taking no arguments and returning a void value.
            /*! Normalise the basis.
             *  Can be delete. It was developpe for LCAO. */
        void normaliseBasis();

        /**
         * @brief Evaluates the molecular orbital values (phi) at a given point.
         * 
         * @param[in] coordinates An array containing the (x, y, z) coordinates
         * 
         * @return double The value of the molecular orbitals at the given point.
         */
        double func(const std::array<double, 3>& coordinates) const;

            //! A normal member taking no arguments and returning an int value.
            /*! \return The HOMO Molecular Orbital number. */
        void HOMO();

            //! A normal member taking no arguments and returning an int value.
            /*! \return The LUMO Molecular Orbital number. */
        void LUMO();

            //! A normal member taking no arguments and returning a std::vector<std::vector<double>> value.
            /*! \return The matrix of overlaps. */
        std::vector<std::vector<double>> get_S();

        //! A normal member taking no arguments and returning a std::vector<double> value.
        /*! \return The table of f value. */
        std::vector<double> get_f(int orb, int alpha=0);

            //! A normal member taking one argument and returning a void value.
            /*! Actualise _all_f value. */
        void get_f(int alpha=0);

            //! A normal member taking no arguments and returning a void value.
            /*! Actualise _all_f for HOMO and LUMO Orbitals. */
        void init_homoLumoIndexes();

            //! A normal member taking no arguments and returning a void value.
            /*! Print all the descriptors with the FMO method. */
        void printDescriptors();

            //! A normal member taking no arguments and returning a void value.
            /*! Change the HOMO on the orbital i and the LUMO on the orbitals j.
             *  Print all the descriptors with the FMO method. */
        void printDescriptors(int homoIndex, int lumoIndex);


        void makeDensityGrid(Grid& grid, const std::vector<std::vector<std::vector<double>>>& reducedDensityMatrix, bool showProgress = false) const;

        /**
         * @brief Creates a grid that stores electronic density (only for Ground State).
         * 
         * @param[in] domain The domain of the grid.
         * @param[in] showProgress Whether to display a progress bar during grid creation.
         */
        Grid makeGrid(const Domain& domain, bool showProgress = false);

        /**
         * @brief Computes the total electronic density for all orbitals at a given point (only for Ground State).
         *
         * @param[in] coordinates The (x, y, z) coordinates at which to compute the density.
         *
         * @return double Total electronic density at the given point.
         */
        double density(const std::array<double, 3>& coordinates) const;

        /**
         * @brief Computes the electronic density for a given spin-orbital at a given point (only for Ground State).
         *
         * @param[in] orbitalIndex The index of the orbital (0-based) for which to compute the electronic density.
         * @param[in] spinType The spin type of the orbital for which to compute the density.
         * @param[in] coordinates The (x, y, z) coordinates at which to compute the density.
         *
         * @return double Electronic density for the specified spin-orbital at the given point.
         */
        double density(int orbitalIndex, SpinType spinType, const std::array<double, 3>& coordinates) const;

        /**
         * @brief Computes the total electronic density for given orbitals at a given point (only for Ground State).
         *
         * @param[in] orbitalIndexes Vector of orbital indexes (0-based) for which to compute the density.
         * @param[in] orbitalSpins Vector of @ref SpinType for each orbital (must have the same length as orbitalIndexes).
         * @param[in] coordinates The (x, y, z) coordinates at which to compute the density.
         *
         * @return double Total electronic density for the given spin-orbitals at the given point.
         */
        double density(const std::vector<int>& orbitalIndexes, const std::vector<SpinType>& orbitalSpins, const std::array<double, 3>& coordinates) const;

        /**
         * @brief Computes the electronic density at coordinates (x,y,z).
         *
         * @param reducedDensityMatrix The matrix of coefficients for calculation of the density
         * @param spinType Spin type (alpha, beta, or both) for which to compute the density.
         * @param coordinates Array containing the (x, y, z) coordinates at which to evaluate the density.
         * 
         * @return double The electronic density evaluated at (x, y, z).
         */
        double density(const std::vector<std::vector<std::vector<double>>>& reducedDensityMatrix, SpinType spinType, const std::array<double, 3>& coordinates) const;

        /**
         * @brief Builds the Grid associated with the orbitals object.
         * 
         * @param[in] domain The domain of the grid.
         * @param[in] orbitalsNumbers The numbers of the orbitals to be included in the grid.
         * @param[in] orbitalsSpins The spins of the orbitals to be included in the grid.
         * @param[in] showProgress If true, shows progress of the grid creation.
         */
        Grid makeOrbGrid(const Domain& domain, const std::vector<int>& orbitalsNumbers, const std::vector<SpinType>& orbitalsSpins, bool showProgress = false);

            //! Electronic density
            /*! Calculates and returns the electronic density from molecular orbitals */
        std::vector<double> phis(const std::array<double, 3>& coordinates, const std::vector<int>& nums, const std::vector<SpinType>& typesSpin);
            
            //! Electron localisation function
            /*! Calculates and returns the ELF*/
        double ELF(const std::array<double, 3>& coordinates, double epsilon=2.87e-5);

            //! Make ELF grid
            /*! Make an ELF grid using Orbitals::ELF()*/
        /**
         * @brief Creates a grid that stores the Electron Localization Function (ELF) values.
         * 
         * @param[in] domain The domain of the grid.
         * @param[in] elfMethod The method to compute ELF values (default is ELFMethod::SAVIN).
         * @param[in] showProgress If true, shows progress of the grid creation.
         * 
         * @return Grid A grid containing the ELF values computed using the specified method.
         */
        Grid makeELFgrid(const Domain& domain, ELFMethod elfMethod = ELFMethod::SAVIN, bool showProgress = false);


            //! A normal member taking no arguments and returning a void value.
            /*! Sort the CGTF. 
             *  This method need to be delete later. 
             *  This was developped to solve a problem when you want save a .wfx file to .molden or .gab. */
        void Sorting();

            //! A normal member taking no arguments and returning a void value.
            /*! Denormalise all the CGTF which compose the Orbitals. */
        void DenormaliseAllBasis();

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file (tag = name.format). */
        void Save(const std::string& tag);

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file .wfx. */
        void Save_wfx(const std::string& tag);

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file .molden. */
        void Save_molden(const std::string& tag);

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file .gab. */
        void Save_gab(const std::string& tag);


        //----------------------------------------------------------------------------------------------------//
        // OPERATOR OVERLOADS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Overloads the output stream redirection operator for an Orbitals.
         *
         * Prints all the data of the Orbitals object.
         *
         * @param stream Output stream.
         * @param orbitals Orbitals to print.
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& stream, const Orbitals& orbitals);
};

    //! An operator member taking two arguments and returning a double value.
    /*! coord = (x,y,z). 
     *  \return The value of Orbitals * Orbitals at the point (x,y,z). */
double operator*(const Orbitals& a, const std::array<double, 3>& coordinates);

#endif
