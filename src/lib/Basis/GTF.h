#ifndef CDFTT_GTF_H_INCLUDED
#define CDFTT_GTF_H_INCLUDED

#include <array>
#include <vector>

#include <Utils/Utils.h>


    //! GTF (Gaussian-Type Function) class.
    /*! This class will be used in the CGTF (Contracted Gaussian-Type Function) class. */

class GTF
{
    private:
        double _exponent;
        double _coefficient;
        std::array<double, 3> _coord;
        std::vector<int> _l;
        Binomial _bino;


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

            //! A default constructor.
            /*! This constructor is used to set all of the parameters for one GTF on 0 or "None" value. */
        GTF();

            //! A real constructor.
            /*! This constructor is used to add all of the parameters for one GTF. */
        GTF(const double exponent, const double coefficient, const std::array<double, 3>& coord, const std::vector<int>& l, const Binomial& binomial);


        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the exponent value.
         */
        double get_exponent() const;

        /**
         * @brief Returns the coefficient value.
         */
        double get_coefficient() const;

        /**
         * @brief Returns the coordinates of the atom associated with the GTF.
         */
        const std::array<double, 3>& get_coord() const;

        /**
         * @brief Returns the quantum numbers associated with the GTF.
         */
        const std::vector<int>& get_l() const;

        /**
         * @brief Returns the binomial helper object associated with the GTF.
         */
        Binomial& get_bino();


        //----------------------------------------------------------------------------------------------------//
        // SETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Sets the coefficient value.
         * @param[in] coefficient The new coefficient value.
         */
        void set_coefficient(const double coefficient);

        /**
         * @brief Sets the exponent value.
         * @param[in] exponent The new exponent value.
         */
        void set_exponent(const double exponent);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        void setL(const int index, const int value)
        {
            _l[index] = value;
        }


        //! A normal member taking one argument and returning a double value.
        /*! \return norme value. */
        double normeGTF();

            //! A normal member taking one argument and returning a double value.
            /*! \return norme value. */
        double normeGTF(GTF& q);

            //! A normal member taking one argument.
            /*! This member normalise the radial's composant of the GTF. */
        void normaliseRadialGTF();

            //! A normal member taking one argument.
            /*! This member unnormalise the radial's composant of the GTF. */
        void denormaliseRadialGTF();

            //! A normal member taking one argument.
            /*! This member normalise the GTF. */
        void normaliseGTF();

            //! A normal member taking one argument and returning a double value.
            /*! \return The GTF overlap. */
        double overlapGTF(GTF&);

            //! A normal member taking two arguments and returning a double value.
            /*! \return The GTF overlap. */
        double overlap3GTF(GTF&, GTF&);

            //! A normal member taking three arguments and returning a double value.
            /*! \return The GTF overlap. */
        double overlap4GTF(GTF&, GTF&, GTF&);

            //! A normal member taking one argument and returning a double value.
            /*! \return The result of an integral with two GTF. */
        double GTFstarGTF(GTF&);
        
            //! A normal member taking two arguments and returning a GTF value.
            /*! \return The result of an integral with three GTF. */
        double GTFstarGTFstarGTF(GTF&, GTF&);
        
            //! A normal member taking three arguments and returning a GTF value.
            /*! \return The result of an integral with four GTF. */
        double GTFstarGTFstarGTFstarGTF(GTF&, GTF&, GTF&);
        
            //! A normal member taking four arguments and returning a double value.
            /*! \return The result of an integral with three ???. */
        double GTFxyzGTF(GTF&, int, int, int);

        /**
         * @brief Calculates the kinetic energy integral < thisGTF | -1/2 nabla^2 | otherGTF >.
         *
         * Ref : Modern techniques in computational chemistry: Motecc -90
         * Enrico Clementi.
         * ESCOM, 1990. ISBN: 90-72199-07-3.
         * Page 392, Eq. (42).
         * 
         * @param[in] otherGTF The other GTF.
         */
        double kineticGTF(const GTF& otherGTF);

        /**
         * @brief Calculates the ionic potential integral < thisGTF | charge / |r - position| | otherGTF >.
         *
         * Ref : Modern techniques in computational chemistry: Motecc -90
         * Enrico Clementi.
         * ESCOM, 1990. ISBN: 90-72199-07-3.
         * Page 394, Eq. (48).
         *
         * @param[in] otherGTF The other GTF.
         * @param[in] position Position of the ion.
         * @param[in] charge Charge of the ion.
         *
         * @return The ionic potential integral.
         */
        double ionicPotentialGTF(const GTF& otherGTF, const std::array<double, 3>& position, const double charge, bool debug = false);
        
            //! A normal member taking three arguments and returning a double value.
            /*! \return The eri ???. */
        double ERIGTF(const GTF&, const GTF&,const GTF&);
            
        /**
         * @brief Evaluates the GTF at given coordinates.
         * 
         * @param[in] coordinates An array containing the (x, y, z) coordinates at which to evaluate the GTF.
         * 
         * @return double The value of the GTF at the given coordinates.
         */
        double func(const std::array<double, 3>& coordinates) const;
        
            //! An operator member taking one argument and returning a void value.
        void operator*=(double);
        
            //! An operator member taking one argument and returning a void value.
        void operator/=(double);

            //! A normal membre taking five arguments and returning a void value.
            /*! Insert all the data in the GTF. */
        void push_back(const double&, const double&, const std::array<double, 3>&, const std::vector<int>&, const Binomial&);
        
        /**
         * @brief Computes the i-th component of the gradient of the GTF at the given coordinates.
         * 
         * @param[in] coordinates The (x, y, z) coordinates at which to compute the gradient.
         * @param[in] i The index of the component of the gradient to compute (0 for x, 1 for y, 2 for z).
         * 
         * @return double The value of the i-th component of the gradient of the GTF at the given coordinates.
         */
        double grad_GTF(const std::array<double, 3>& coordinates, int i);

};

    //! An operator member taking two arguments and returning a bool value.
    /*! \return The bool value of an equality between two CGTF. */
bool operator==(GTF, GTF);

    //! An operator member taking two arguments and returning an std::ostream value.
    /*! Print all the data of one GTF */
std::ostream& operator<<(std::ostream&, const GTF& gtf);

/**
 * @brief Overloads the multiplication operator for a vector of GTFs and a set of coordinates, returning the product of the GTFs evaluated at those coordinates.
 * 
 * @param gtfs std::vector of GTF objects.
 * @param coordinates An array containing the (x, y, z) coordinates at which to evaluate the product of the GTFs.
 * 
 * @return double The product of the GTFs evaluated at the given coordinates.
 */
double operator*(const std::vector<GTF>& gtfs, const std::array<double, 3>& coordinates);

#endif
