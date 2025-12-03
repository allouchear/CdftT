#ifndef CDFTT_ATOM_H_INCLUDED
#define CDFTT_ATOM_H_INCLUDED

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
        double _coordinates[3];

        /** @brief Gradient vector. */
        double _gradient[3];

        /** @brief Velocity vector. */
        double _velocity[3];

        /** @brief Element name of the atom. */
        std::string _name;

        /** @brief Chemical symbol of the atom (e.g. "H", "C"). */
        std::string _symbol;

        /** @brief Atomic number (Z) of the atom. */
        int _atomic_number;

        /** @brief Partial charge of the atom. */
        double _charge;

        /** @brief Oxidation state of the atom. */
        double _charge_0;

        /** @brief Hardness (eta) of the atom. */
        double _hardness; // eta

        /** @brief ? */
        double _width; // eta

        /** @brief Covalent radii of the atom. */
        double _covalent_radii;

        /** @brief Associated Element object (properties dictionary). */
        Element _e;

        /** @brief Molecular mechanics atom type. */
        std::string _mm_Type;

        /** @brief PDB atom type. */
        std::string _pdb_Type;

        /** @brief Residue name of the atom. */
        std::string _residue_name;

        /** @brief Residue number of the atom. */
        int _residue_number;

        /** @brief Index or count value. */
        int _N;


    public:
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
         * @param Table PeriodicTable used to find the element.
         * @param name Name or symbol of the element.
         */
        Atom(const PeriodicTable& Table, const std::string& name);

        /**
         * @brief Constructor.
         *
         * Creates an atom from the atomic number of an Element and searching for it in PeriodicTable
         *
         * @param Table PeriodicTable used to find the element
         * @param n Atomic number
         */
        Atom(const PeriodicTable& Table, const int Z);

        /**
         * @brief Returns the coordinates of the atom.
         *
         * Returns the coordinates of the atom (3-element array).
         *
         * @return const double* 3-element array of coordinates.
         */
        const double* coordinates() const;

        /**
         * @brief Returns the gradient of the atom.
         *
         * Returns the gradient of the atom (3-element array).
         *
         * @return double* 3-element gradient array.
         */
        double* gradient();

        /**
         * @brief Returns the velocity of the atom.
         *
         * Returns the velocity of the atom (3-element array).
         *
         * @return double* 3-element velocity array.
         */
        double* velocity();

        /**
         * @brief Returns the name of the atom.
         *
         * @return string Name of the atom.
         */
        std::string name() const;

        /**
         * @brief Returns the chemical symbol of the atom.
         *
         * @return string Symbol of the atom.
         */
        std::string symbol() const;

        /**
         * @brief Returns the atomic number (Z) of the atom.
         *
         * @return int Atomic number (Z).
         */
        int atomic_number() const;

        /**
         * @brief Returns the partial charge of the atom.
         *
         * @return double Partial charge.
         */
        double charge() const;

        /**
         * @brief Returns the oxidation state of the atom.
         *
         * @return double Oxidation state.
         */
        double charge_0() const;

        /**
         * @brief Returns the hardness of the atom.
         *
         * @return double Hardness (eta).
         */
        double hardness() const;

        /**
         * @brief Returns the ? of the atom.
         *
         * @return double ?.
         */
        double width() const;

        /**
         * @brief Returns the covalent radii of the atom.
         *
         * @return double Covalent radii.
         */
        double covalent_radii() const;

        /**
         * @brief Returns the Element associated with the atom.
         *
         * Returns the element associated with the atom.
         * See @ref Element for full documentation of the Element class.
         *
         * @return Element The associated Element object.
         */
        Element element();

        /**
         * @brief Sets the i-th coordinate to d.
         *
         * @param i Index of coordinate (0..2).
         * @param d Value to set.
         */
        void set_coordinates(const int i, const double d);

        /**
         * @brief Sets the i-th gradient component to d.
         *
         * @param i Index of gradient component (0..2).
         * @param d Value to set.
         */
        void set_gradient(const int i, const double d);

        /**
         * @brief Sets the i-th velocity component to d.
         *
         * @param i Index of velocity component (0..2).
         * @param d Value to set.
         */
        void set_velocity(const int i, const double d);

        /**
         * @brief Sets the partial charge.
         *
         * @param c Charge value to set.
         */
        void set_charge(const double c);

        /**
         * @brief Default destructor.
         *
         * Not used explicitly.
         */
        ~Atom();

        /**
         * @brief Returns the distance between this atom and another atom.
         *
         * @param a2 Other atom.
         * @return double Distance (in coordinates units).
         */
        double get_distance(const Atom& a2) const;

        /**
         * @brief Returns the angle (in degrees) formed by this atom and two other atoms (a2, a3).
         *
         * @param a2 Second atom.
         * @param a3 Third atom.
         * @return double Angle in degrees.
         */
        double get_angle(const Atom& a2, const Atom& a3) const;

        /**
         * @brief Returns the torsion (dihedral) angle (in degrees) defined by four atoms (this atom, a2, a3, a4).
         *
         * @param a2 Second atom.
         * @param a3 Third atom.
         * @param a4 Fourth atom.
         * @return double Torsion angle in degrees.
         */
        double get_torsion(const Atom& a2, const Atom& a3, const Atom& a4) const;
};

#endif // CDFTT_ATOM_H_INCLUDED