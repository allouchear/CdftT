#ifndef CDFTT_ORBITALS_H_INCLUDED
#define CDFTT_ORBITALS_H_INCLUDED

#include <iostream>
#include <string>
#include <vector>

#include <Basis/CGTF.h>
#include <Common/Descriptors.h>
#include <Common/PeriodicTable.h>
#include <Common/Structure.h>
#include <Utils/Enums.hpp>
#include <Utils/FCHK.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/WFX.h>



    //! An Orbitals class.
    /*! This class will be used to calculate descriptors. */
class Orbitals
{
    private:
        std::vector<CGTF> _vcgtf;
        std::vector<std::vector<std::vector<double>>> _coefficients;
        int _numberOfAo;
        int _numberOfMo;
        int _numberOfAlphaElectrons;
        int _numberOfBetaElectrons;
        int _numberOfAtoms;
        std::vector<int> _primitiveCenters;
        Structure _struct;
        std::vector<int> _atomicNumbers;
        std::vector<std::string> _symbol;
        std::vector<std::vector<double>> _orbitalEnergy;
        std::vector<std::vector<double>> _all_f;
        std::vector<int> _numOrb;
        std::vector<std::vector<double>> _occupationNumber;
        bool _alphaAndBeta;
        Binomial _bino;
        Descriptors _descriptors;
        std::vector<CGTF> _vcgtfNonNormalise;
        int _numberOfGtf;
        double _energy;
        std::vector<double> _coordinates;
        bool _mixte;


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

            //! A default constructor.
            /*! This constructor is used to set all of the parameters for Orbitals on 0 or "None" value. */
        Orbitals();

            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .wfx file. */
        Orbitals(WFX& wfxParser, Binomial& bino, const PeriodicTable& periodicTable);

            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .fchk file. */
        Orbitals(FCHK& fchkParser, Binomial& bino, const PeriodicTable& periodicTable);

            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .molden or .gab file. */
        Orbitals(MOLDENGAB& moldengabParser, Binomial& bino, const PeriodicTable& periodicTable);
            //! A real constructor.
            /*! This constructor is used to add all of the parameters for Orbitals with the data in .log file. */
        Orbitals(LOG& logParser, Binomial& bino, const PeriodicTable& periodicTable);


        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

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
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the i-th primitive center.
         */
        int getPrimitiveCenter(int i) const;
        
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
         * @brief Calculates the ion-electron integral for a given ion : < phi_i | V_ion | phi_j >
         *
         * @param[in] chargePosition The position of the ion.
         * @param[in] charge The charge of the ion.
         * @return The matrix < phi_i | V_ion | phi_j > (the first index corresponds to alpha spin, the second to beta spin).
         */
        std::vector<std::vector<std::vector<double>>> getIonicPotentialMatrix(const std::array<double, 3>& chargePosition, double charge, bool debug = false);

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

            //! A normal member taking three arguments and returning a double value.
            /*! \return The value of the Orbitals at the point (x,y,z). */
        double func(double x, double y, double z) const;

            //! A normal member taking no arguments and returning an int value.
            /*! \return The HOMO Molecular Orbital number. */
        void HOMO();

            //! A normal member taking no arguments and returning an int value.
            /*! \return The LUMO Molecular Orbital number. */
        void LUMO();

            //! A normal member taking no arguments and returning a double value.
            /*! \return The HOMO's energy. */
        double eHOMO(int alpha=0) const {return _orbitalEnergy[alpha][_numOrb[0]];}

            //! A normal member taking no arguments and returning a double value.
            /*! \return The LUMO's energy. */
        double eLUMO(int alpha=0) const {return _orbitalEnergy[alpha][_numOrb[1]];}

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
        void HOMO_LUMO();

            //! A normal member taking two arguments and returning a void value.
            /*! Actualise _all_f for i and j Orbitals. */
        void HOMO_LUMO(int i, int j);

            //! A normal member taking no arguments and returning a void value.
            /*! Print all the descriptors with the FMO method. */
        void PrintDescriptors();

            //! A normal member taking no arguments and returning a void value.
            /*! Change the HOMO on the orbital i and the LUMO on the orbitals j.
             *  Print all the descriptors with the FMO method. */
        void PrintDescriptors(int i, int j);


        
            //! Make a grid of electronic density
            /*! creates a grid of electronic density. Values are calculated with Orbitals::density(x, y, z)*/
        Grid makeGrid(const Domain& d);

            //! Electronic density
            /*! Calculates and returns the electronic density from molecular orbitals */
        double density(double x, double y, double z);
            
            //! Make a grid of Molecular orbitals
            /*! Creates a grid of molecular orbitals. Values are calculated with Orbitals::phis()*/
        Grid makeOrbGrid(const Domain& d, const std::vector<int>& nums, const std::vector<int>& typesSpin);

            //! Electronic density
            /*! Calculates and returns the electronic density from molecular orbitals */
        std::vector<double> phis(double x, double y, double z, const std::vector<int>& nums, const std::vector<int>& typesSpin);
            
            //! Electron localisation function
            /*! Calculates and returns the ELF*/
        double ELF(const double& x, const double& y, const double& z, double epsilon=2.87e-5);

            //! Make ELF grid
            /*! Make an ELF grid using Orbitals::ELF()*/
        Grid makeELFgrid(const Domain& d,const double& epsilon=2.87e-5);


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
        void Save(std::string& tag);

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file .wfx. */
        void Save_wfx(std::string& tag);

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file .molden. */
        void Save_molden(std::string& tag);

            //! A normal member taking one argument and returning a void value.
            /*! Save all the data in Orbitals in a file .gab. */
        void Save_gab(std::string& tag);


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
        friend std::ostream& operator<<(std::ostream& stream, Orbitals& orbitals);
};

    //! An operator member taking two arguments and returning a double value.
    /*! coord = (x,y,z). 
     *  \return The value of Orbitals * Orbitals at the point (x,y,z). */
double operator*(const Orbitals& a, const std::vector<double>& coord);

#endif
