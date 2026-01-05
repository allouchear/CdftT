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
        
            //! A normal member taking one argument and returning a double value.
            /*! \return The kinetic ???. */
        double kineticGTF(GTF&);
        
            //! A normal member taking three arguments and returning a double value.
            /*! \return The ionic potential ???. */
        double ionicPotentialGTF(const GTF&, const std::array<double, 3>&, double);
        
            //! A normal member taking three arguments and returning a double value.
            /*! \return The eri ???. */
        double ERIGTF(GTF&, GTF&, GTF&);
            
        double func(double x, double y, double z) const;
        
            //! An operator member taking one argument and returning a void value.
        void operator*=(double);
        
            //! An operator member taking one argument and returning a void value.
        void operator/=(double);

            //! A normal membre taking five arguments and returning a void value.
            /*! Insert all the data in the GTF. */
        void push_back(const double&, const double&, const std::array<double, 3>&, const std::vector<int>&, Binomial&);
        
            //! Gradient of GTF
            /*! Get the ith component of the gradient of a GTF*/
        double grad_GTF(const double& x, const double& y, const double& z, int i);

};

    //! An operator member taking two arguments and returning a bool value.
    /*! \return The bool value of an equality between two CGTF. */
bool operator==(GTF, GTF);

    //! An operator member taking two arguments and returning an std::ostream value.
    /*! Print all the data of one GTF */
std::ostream& operator<<(std::ostream&, GTF&);

    //! An operator member taking two arguments and returning a double value.
    /*! \return The double value of a product between a std::vector of GTF at the coordinates x,y,z*/
double operator*(const std::vector<GTF>&, const std::vector<double>&);

#endif
