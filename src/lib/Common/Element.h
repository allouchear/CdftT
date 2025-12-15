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
         * @brief Default destructor.
         *
         * Not used explicitly.
         */
        ~Isotope();

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
        int _atomic_number;

        /** @brief Covalent radius of the element (in angstroms or project units). */
        double _covalent_radii;

        /** @brief Bond order adjusted radius. */
        double _bond_order_radii;

        /** @brief Van der Waals radius. */
        double _van_der_waals_radii;

        /** @brief Generic atomic radius used by the project. */
        double _radii;

        /** @brief Maximum number of bonds allowed for this element. */
        int _maximum_bond_valence;

        /** @brief Atomic mass of the element. */
        double _mass;

        /** @brief Electronegativity (Pauling scale) of the element. */
        double _electronegativity;

        /** @brief List of isotopes associated with the element. */
        std::vector<Isotope> _isotope;


    public:
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

        /**
         * @brief Default destructor.
         *
         * Not used explicitly.
         */
        ~Element();

        /**
         * @brief Returns the name.
         */
        std::string name() const
        {
            return _name;
        }

        /**
         * @brief Returns the symbol.
         */
        std::string symbol() const
        {
            return _symbol;
        }

        /**
         * @brief Returns the atomic number.
         */
        int atomic_number() const
        {
            return _atomic_number;
        }

        /**
         * @brief Returns the covalent radius.
         */
        double covalent_radii() const
        {
            return _covalent_radii;
        }

        /**
         * @brief Returns the bond order radius.
         */
        double bond_order_radii() const
        {
            return _bond_order_radii;
        }

        /**
         * @brief Returns the van der Waals radius.
         */
        double van_der_waals_radii() const
        {
            return _van_der_waals_radii;
        }

        /**
         * @brief Returns the generic radius.
         */
        double radii() const
        {
            return _radii;
        }

        /**
         * @brief Returns the maximum number of bonds allowed for this element.
         */
        int maximum_bond_valence() const
        {
            return _maximum_bond_valence;
        }

        /**
         * @brief Returns the atomic mass in atomic mass units.
         */
        double mass() const
        {
            return _mass;
        }

        /**
         * @brief Returns the electronegativity (Pauling scale).
         */
        double electronegativity() const
        {
            return _electronegativity;
        }

        /**
         * @brief Returns the given isotope of the element.
         *
         * @param i Index of the isotope (starting from 1).
         * @return The requested Isotope object.
         */
        Isotope isotope(const int i) const
        {
            return _isotope[i - 1];
        }

        /**
         * @brief Adds an isotope to the element.
         *
         * @param isotope Isotope object to append to the internal isotope list.
         */
        void push_isotope(const Isotope& isotope)
        {
            _isotope.push_back(isotope);
        }
};


#endif
