#ifndef CDFTT_CGTF_H_INCLUDED
#define CDFTT_CGTF_H_INCLUDED

#include <array>
#include <ostream>
#include <string>
#include <vector>

#include <Basis/GTF.h>


/**
 * @brief Contracted Gaussian-Type Function (CGTF) class.
 *
 * This class represents a contracted Gaussian-type function and is used in
 * LCAO (Linear Combination of Atomic Orbitals) constructions.
 */
class CGTF
{
    private:
        /** @brief Number of the center (atom index) associated with this CGTF. */ // ?
        int _num_center;

        /** @brief Number of Gaussian functions composing this CGTF. */
        int _numberOfFunctions;

        /** @brief Orbital shell type ("S", "P", "D", ...). */
        std::string _l_type;

        /** @brief Orbital format ("cartesian" or "spherical"). */
        std::string _l_format;

        /** @brief Global multiplicative factor for coefficients. */
        double _factor_coef;

        /** @brief Primitive GTFs composing the CGTF. */
        std::vector<GTF> _gtf;

        /** @brief Contraction coefficients for each primitive GTF. */
        std::vector<double> _coefficients;

        /** @brief Binomial handler class used by primitive GTFs. */
        Binomial _bino;

    public:
        /**
         * @brief Default constructor.
         *
         * Sets members to empty / zero values.
         */
        CGTF();

        /**
         * @brief Constructor from a std::vector of primitive GTF.
         *
         * @param gtfs std::vector of primitive GTF objects used to form the CGTF.
         */
        CGTF(std::vector<GTF> gtfs);

        /**
         * @brief Default destructor.
         *
         * Not used explicitly.
         */
        ~CGTF(){}
        
        /**
         * @brief Appends a coefficient for a primitive GTF.
         *
         * @param c Coefficient to append.
         */
        void setCoef(double);

        /**
         * @brief Sets the orbital format for this CGTF.
         *
         * @param format Format std::string (either "cartesian" or "spherical").
         */
        void setFormat(std::string);

        /**
         * @brief Returns the contraction coefficients.
         *
         * @return std::vector of coefficients for the primitive GTFs.
         */
        std::vector<double> coefficients() const {return _coefficients;}

        /**
         * @brief Returns the orbital type ("S", "P", "D", ...).
         *
         * @return Orbital type.
         */
        std::string Ltype() {return _l_type;}

        /**
         * @brief Returns the orbital format ("cartesian" or "spherical").
         *
         * @return Orbital format.
         */
        std::string Lformat() {return _l_format;}

        /**
         * @brief Returns the center (atom) index of this CGTF.
         *
         * @return Index of the center (atom) for this CGTF.
         */
        int NumCenter() {return _num_center;}

        /**
         * @brief Returns the global multiplicative factor for the coefficients.
         *
         * @return Global multiplicative factor for the coefficients.
         */
        double FactorCoef() {return _factor_coef;}

        /**
         * @brief Returns the number of primitive GTFs in this CGTF.
         *
         * @return Number of primitive GTFs composing the CGTF.
         */
        int numberOfFunctions() const
        {
            return _numberOfFunctions;
        }

        /**
         * @brief Returns the std::vector holding all primitive GTFs composing this CGTF.
         *
         * @return std::vector of primitive GTF objects.
         */
        std::vector<GTF> gtf() const
        {
            return _gtf;
        }

        /**
         * @brief Returns the Binomial object associated with this CGTF.
         *
         * @return Binomial handler class.
         */
        Binomial bino()
        {
            return _bino;
        }

        /**
         * @brief Computes the electron repulsion integral (ERI) involving four CGTFs.
         *
         * Computes the ERI (pq|rs) where p is this CGTF and q,r,s are provided.
         *
         * @param q Second CGTF.
         * @param r Third CGTF.
         * @param s Fourth CGTF.
         * @return Value of the ERI.
         */
        double ERICGTF(CGTF&, CGTF&, CGTF&);

        /**
         * @brief Normalizes the CGTF.
         *
         * Normalizes the CGTF. Each primitive radial function is normalized. Contraction coefficients
         * are then rescaled so that the CGTF has unit norm.
         */
        void normaliseCGTF();

        /**
         * @brief Denormalizes the primitive GTFs constituting the CGTF.
         *
         * This denormalizes the radial component of the constituent primitive GTFs.
         */
        void denormaliseCGTF();

        /**
         * @brief Returns the overlap integral between this CGTF and an other CGTF.
         *
         * @param right The other CGTF to compute the overlap integral with.
         * @return Value of the overlap integral between the two CGTFs.
         */
        double overlapCGTF(CGTF& right);

        /**
         * @brief Returns the overlap integral between this CGTF and two other ones.
         *
         * @param midle The middle CGTF.
         * @param right The right CGTF.
         * @return Value of the overlap integral between the three CGTFs.
         */
        double overlap3CGTF(CGTF& middle, CGTF& right);

        /**
         * @brief Returns the overlap integral between this CGTF and two other ones.
         *
         * @param midle The middle CGTF.
         * @param right The right CGTF.
         * @return Value of the overlap integral between the three CGTFs.
         */
        double overlap4CGTF(CGTF&, CGTF&, CGTF&);

        /**
         * @brief ?
         *
         * ?
         *
         * @param right The other CGTF.
         * @param ix ?
         * @param iy ?
         * @param iz ?
         * @return Integral value.
         */
        double CGTFxyzCGTF(CGTF&, int, int, int);

        /**
         * @brief Calculates the kinetic energy integral < thisCGTF | -1/2 nabla^2 | otherCGTF >.
         * 
         * @param[in] otherCGTF The other CGTF.
         */
        double kineticCGTF(const CGTF& otherCGTF);

        /**
         * @brief Calculates the ionic potential integral < thisCGTF | charge / |r - position| | otherCGTF >.
         * 
         * @param[in] otherCGTF The other CGTF.
         * @param[in] position Position of the ion.
         * @param[in] charge Charge of the ion.
         * 
         * @return The ionic potential integral.
         */
        double ionicPotentialCGTF(const CGTF &otherCGTF, const std::array<double, 3>& position, const double charge, bool debug = false);
        
            //! A normal member taking one argument and returning a double value.
            /*! \return ???. */
        double CGTFstarCGTF(CGTF&);

        //bool CGTFEqCGTF(CGTF&);

            //! A normal membre taking one argument and returning a void value.
            /*! Insert all the data in the CGTF. */
        void push_back(GTF&);

            //! A normal membre taking one argument and returning a void value.
            /*! Insert the num center in the CGTF. */
        void setNumCenter(int);

            //! A normal membre taking one argument and returning a void value.
            /*! Insert the shell type (S, P, D, ...) in the CGTF. */
        void setLtype(std::string);

            //! A normal membre taking one argument and returning a void value.
            /*! Insert the factor coefficient in the CGTF. */
        void setFactorCoef(double);

            //! A normal membre taking three arguments and returning a double value.
            /*! \return The value of the CGTF at the coordinates x,y,z. */
        double func(double x, double y, double z) const;
            
            //! Gradient of a CGTF
            /*! Computes the ith component of the gradient of CGTF*/
        double grad_CGTF(const double& x, const double& y, const double& z, int i);
};

/**
 * @brief Overloads the equality operator for two CGTF objects.
 *
 * Two CGTFs are considered equal if they contain the same set of primitive GTFs, independantly of the order.
 *
 * @param a First CGTF.
 * @param b Second CGTF.
 * @return True if the two CGTFs are equal, false otherwise.
 */
bool operator==(const CGTF& left, const CGTF& right);

/**
 * @brief Overloads the output stream redirection operator for a CGTF.
 *
 * Prints coefficients and primitive GTFs composing the CGTF.
 *
 * @param stream Output stream.
 * @param cgtf CGTF to print.
 * @return Reference to the output stream.
 */
std::ostream& operator<<(std::ostream &stream, const CGTF &cgtf);

/**
 * @brief Computes the product of many CGTFs evaluated at a given point (x,y,z).
 *
 * Evaluates each CGTF of the std::vector at the given (x,y,z) coordinates and
 * returns the product of the evaluations.
 *
 * @param cgtfs std::vector of CGTF objects.
 * @param coords 3-element coordinate std::vector (x,y,z).
 * @return Product of the CGTFs values at the given coordinates.
 */
double operator*(const std::vector<CGTF> &cgtfs, const std::vector<double> &coords);

#endif
