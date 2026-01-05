#ifndef CDFTT_ELEMENT_H_INCLUDED
#define CDFTT_ELEMENT_H_INCLUDED

#include <string>
#include <vector>


/**
 * @brief Isotope class.
 *
 * This class is used in the Element class to complete the data on an element.
 */
class Isotope
{
    private:
        /** @brief Symbol of the isotope (e.g., "H", "D", "T"). */
        std::string _symbol;

        /** @brief Mass number of the isotope. */
        int _int_mass;

        /** @brief Exact mass of the isotope in atomic mass units. */
        double _real_mass;

        /** @brief Natural abundance (in percent) of the isotope. */
        double _abundance;


    public:
        /**
         * @brief Default constructor.
         *
         * Sets all numeric attributes to 0 and std::string attributes to "None".
         */
        Isotope();

        /**
         * @brief Constructor.
         * 
         * Initializes all data of one isotope.
         */
        Isotope(const std::string&, const int, const double, const double);

        /**
         * @brief Returns the isotope symbol.
         */
        std::string symbol() const
        {
            return _symbol;
        }

        /**
         * @brief Returns the mass number.
         */
        int int_mass() const
        {
            return _int_mass;
        }

        /**
         * @brief Returns the exact mass in atomic mass units.
         */
        double real_mass() const
        {
            return _real_mass;
        }

        /**
         * @brief Returns the natural abundance in percent.
         */
        double abundance() const
        {
            return _abundance;
        }
};

/**
 * @brief Element class.
 *
 * This class is used as a dictionary containing all the properties of an element.
 */
class Element
{
    private:
        /** @brief Full name of the element (e.g., "Hydrogen"). */
        std::string _name;

        /** @brief Chemical symbol of the element (e.g., "H"). */
        std::string _symbol;

        /** @brief Atomic number (Z) of the element. */
        int _atomicNumber;

        /** @brief Covalent radius of the element (in angstroms or project units). */
        double _covalentRadius;

        /** @brief Bond order adjusted radius. */
        double _bondOrderRadius;

        /** @brief Van der Waals radius. */
        double _vanDerWaalsRadius;

        /** @brief Generic atomic radius used by the project. */
        double _radius;

        /** @brief Maximum number of bonds allowed for this element. */
        int _maximumBondValence;

        /** @brief Atomic mass of the element. */
        double _mass;

        /** @brief Electronegativity (Pauling scale) of the element. */
        double _electronegativity;

        /** @brief List of isotopes associated with the element. */
        std::vector<Isotope> _isotope;


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Default constructor.
         *
         * Sets numeric attributes to 0 and std::string attributes to "None".
         */
        Element();

        /**
         * @brief Constructor.
         *
         * Initializes all data for an element.
         */
        Element(const std::string& name, const std::string& symbol, const int atomicNumber, const double covalentRadii, const double bondOrderRadii, const double vanDerWaalsRadii, const double radii, const int maximumBondValence, const double mass, const double electronegativity);


        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the name.
         */
        const std::string& get_name() const;

        /**
         * @brief Returns the symbol.
         */
        const std::string& get_symbol() const;

        /**
         * @brief Returns the atomic number.
         */
        int get_atomicNumber() const;

        /**
         * @brief Returns the covalent radius.
         */
        double get_covalentRadius() const;

        /**
         * @brief Returns the bond order radius.
         */
        double get_bondOrderRadius() const;

        /**
         * @brief Returns the van der Waals radius.
         */
        double get_vanDerWaalsRadius() const;

        /**
         * @brief Returns the generic radius.
         */
        double get_radius() const;

        /**
         * @brief Returns the maximum number of bonds allowed for this element.
         */
        int get_maximumBondValence() const;

        /**
         * @brief Returns the atomic mass in atomic mass units.
         */
        double get_mass() const;

        /**
         * @brief Returns the electronegativity (Pauling scale).
         */
        double get_electronegativity() const;


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the given isotope of the element.
         *
         * @param i Index of the isotope (starting from 1).
         * @return The requested Isotope object.
         */
        Isotope getIsotope(const int i) const;

        /**
         * @brief Adds an isotope to the element.
         *
         * @param isotope Isotope object to append to the internal isotope list.
         */
        void pushIsotope(const Isotope& isotope);
};


#endif
