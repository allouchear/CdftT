#ifndef CDFTT_ATOM_H_INCLUDED
#define CDFTT_ATOM_H_INCLUDED

#include <array>
#include <string>

#include <Common/PeriodicTable.h>
#include <Common/Element.h>


/**
 * @brief Atom class.
 *
 * This class will be used in the Structure class to create molecules or periodic structures.
 */
class Atom
{
    private:
        /** @brief Atom coordinates (x,y,z). */
        std::array<double, 3> _coordinates;

        /** @brief Gradient vector. */
        std::array<double, 3> _gradient;

        /** @brief Velocity vector. */
        std::array<double, 3> _velocity;

        /** @brief Element name of the atom. */
        std::string _name;

        /** @brief Chemical symbol of the atom (e.g. "H", "C"). */
        std::string _symbol;

        /** @brief Atomic number (Z) of the atom. */
        int _atomicNumber;

        /** @brief Partial charge of the atom. */
        double _charge;

        /** @brief Oxidation state of the atom. */
        double _charge0;

        /** @brief Hardness (eta) of the atom. */
        double _hardness; // eta

        /** @brief ? */
        double _width; // eta

        /** @brief Covalent radii of the atom. */
        double _covalentRadius;

        /** @brief Associated Element object (properties dictionary). */
        Element _e;

        /** @brief Molecular mechanics atom type. */
        std::string _mmType;

        /** @brief PDB atom type. */
        std::string _pdbType;

        /** @brief Residue name of the atom. */
        std::string _residueName;

        /** @brief Residue number of the atom. */
        int _residueNumber;

        /** @brief Index or count value. */
        int _N;


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Default constructor.
         *
         * Sets all numeric attributes to 0 and string attributes to "none".
         */
        Atom();

        /**
         * @brief Constructor.
         *
         * Creates an atom from the name of an Element and searching for it in PeriodicTable.
         *
         * @param PeriodicTable PeriodicTable used to find the element.
         * @param name Name or symbol of the element.
         */
        Atom(const PeriodicTable& PeriodicTable, const std::string& name);

        /**
         * @brief Constructor.
         *
         * Creates an atom from the atomic number of an Element and searching for it in PeriodicTable
         *
         * @param periodicTable PeriodicTable used to find the element
         * @param n Atomic number
         */
        Atom(const PeriodicTable& periodicTable, const int Z);

        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the coordinates of the atom.
         */
        const std::array<double, 3>& get_coordinates() const;

        /**
         * @brief Returns the gradient vector of the atom.
         */
        const std::array<double, 3>& get_gradient() const;

        /**
         * @brief Returns the velocity vector of the atom.
         */
        const std::array<double, 3>& get_velocity() const;

        /**
         * @brief Returns the name of the atom.
         */
        std::string get_name() const;

        /**
         * @brief Returns the chemical symbol of the atom.
         *
         * @return string Symbol of the atom.
         */
        std::string get_symbol() const;

        /**
         * @brief Returns the atomic number (Z) of the atom.
         *
         * @return int Atomic number (Z).
         */
        int get_atomicNumber() const;

        /**
         * @brief Returns the partial charge of the atom.
         *
         * @return double Partial charge.
         */
        double get_charge() const;

        /**
         * @brief Returns the oxidation state of the atom.
         *
         * @return double Oxidation state.
         */
        double get_charge0() const;

        /**
         * @brief Returns the hardness of the atom.
         *
         * @return double Hardness (eta).
         */
        double get_hardness() const;

        /**
         * @brief Returns the ? of the atom.
         *
         * @return double ?.
         */
        double get_width() const;

        /**
         * @brief Returns the covalent radii of the atom.
         *
         * @return double Covalent radii.
         */
        double get_covalentRadius() const;

        /**
         * @brief Returns the Element associated with the atom.
         *
         * Returns the element associated with the atom.
         * See @ref Element for full documentation of the Element class.
         *
         * @return Element The associated Element object.
         */
        Element get_element();


        //----------------------------------------------------------------------------------------------------//
        // SETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Sets the partial charge.
         */
        void set_charge(const double c);


        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the angle (in degrees) formed by this atom and two other atoms (a2, a3).
         *
         * @param a2 Second atom.
         * @param a3 Third atom.
         * @return double Angle in degrees.
         */
        double computeAngle(const Atom& a2, const Atom& a3) const;

        /**
         * @brief Computes the distance between this atom and another atom.
         *
         * @param a2 Other atom.
         * @return double Distance (in coordinates units).
         */
        double computeDistance(const Atom& a2) const;

        /**
         * @brief Computes the distance between this atom and a given coordinate.
         *
         * @param distantCoordinate Coordinate to compute the distance to.
         * @return double Distance (in coordinates units).
         */
        double computeDistance(const std::array<double, 3>& distantCoordinate) const;

        /**
         * @brief Computes the torsion (dihedral) angle (in degrees) defined by four atoms (this atom, a2, a3, a4).
         *
         * @param a2 Second atom.
         * @param a3 Third atom.
         * @param a4 Fourth atom.
         * @return double Torsion angle in degrees.
         */
        double computeTorsion(const Atom& a2, const Atom& a3, const Atom& a4) const;

        /**
         * @brief Sets the i-th coordinate to d.
         *
         * @param i Index of coordinate (0..2).
         * @param d Value to set.
         */
        void setCoordinate(const int i, const double d);

        /**
         * @brief Sets the i-th gradient component to d.
         *
         * @param i Index of gradient component (0..2).
         * @param d Value to set.
         */
        void setGradientComponent(const int i, const double d);

        /**
         * @brief Sets the i-th velocity component to d.
         *
         * @param i Index of velocity component (0..2).
         * @param d Value to set.
         */
        void setVelocityComponent(const int i, const double d);
};

#endif // CDFTT_ATOM_H_INCLUDED