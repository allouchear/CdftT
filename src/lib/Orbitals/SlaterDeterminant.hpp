#ifndef CDFTT_SLATERDETERMINANT_HPP_INCLUDED
#define CDFTT_SLATERDETERMINANT_HPP_INCLUDED

#include <iostream>
#include <utility>
#include <vector>

#include <Orbitals/Orbitals.h>
#include <Utils/Enums.hpp>


/** @brief SlaterDeterminant class.
 * 
 * Manages Slater determinants representing electronic configurations.
 */
class SlaterDeterminant
{
    private:
        /** @brief Type alias for the occupation of an orbital (orbital number, occupation number). */
        typedef std::pair<int, double> Occupation;

        /** @brief Occupied orbitals (1-based) and their occupation numbers. First index corresponds to alpha spin, second to beta spin. */
        std::vector<std::vector<Occupation>> _occupiedOrbitals;


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Default constructor.
         * 
         * Initializes empty vectors for occupied orbitals and occupation numbers.
         */
        SlaterDeterminant();

        /**
         * @brief Constructor for a ground state Slater determinant.
         * 
         * @param orbitals Orbitals instance from which to build the Slater determinant.
         */
        SlaterDeterminant(const Orbitals& orbitals);


        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Returns the occupied orbitals (1-based) and their occupation numbers. The first index corresponds to alpha spin, the second to beta spin.
         */
        const std::vector<std::vector<std::pair<int, double>>>& get_occupiedOrbitals() const;

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Updates the Slater determinant based on an electronic transition.
         * 
         * @param[in] initialOrbitalNumber The orbital number from which an electron is removed.
         * @param[in] initialSpin The spin type (alpha or beta) of the electron being removed.
         * @param[in] finalOrbitalNumber The orbital number to which an electron is added.
         * @param[in] finalSpin The spin type (alpha or beta) of the electron being added.
         * 
         * @return True if the update was successful, false otherwise (e.g., if the transition is invalid).
         */
        bool updateFromTransition(int initialOrbitalNumber, SpinType initialSpin, int finalOrbitalNumber, SpinType finalSpin);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the differences in occupied orbitals between two Slater determinants.
         * 
         * @param d_i First Slater determinant.
         * @param d_j Second Slater determinant.
         * 
         * @return A vector containing the differences in occupied orbitals for each spin type (first index corresponds to alpha spin, second to beta spin)
         */
        static std::vector<std::vector<std::pair<int, int>>> getDifferences(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j);

        /**
         * @brief Determines whether two Slater determinants are equivalent (i.e., have the same occupied orbitals and occupation numbers, without accounting for orbital ordering).
         * 
         * @param d_i First Slater determinant.
         * @param d_j Second Slater determinant.
         * 
         * @return True if the Slater determinants are equivalent, false otherwise.
         */
        static bool equivalent(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j);

        /**
         * @brief Computes the overlap between two Slater determinants.
         *
         * @param[in] d_i First Slater determinant.
         * @param[in] d_j Second Slater determinant.
         * 
         * @return The overlap < D_i | D_j > between the two Slater determinants.
         */
        static double overlap(const SlaterDeterminant& di, const SlaterDeterminant& dj);

        /**
         * @brief Computes the ionic potential matrix element between two Slater determinants.
         *
         * @param[in] d_i First Slater determinant.
         * @param[in] d_j Second Slater determinant.
         * @param[in] ionicMatrix Ionic matrix < phi_i | V_ion/electrons | phi_j > (the first index corresponds to alpha spin, the second to beta spin).
         * 
         * @return The ionic potential matrix element < D_i | V_ion/electrons | D_j >.
         */
        static double ionicPotential(const SlaterDeterminant& d_i, const SlaterDeterminant& d_j, const std::vector<std::vector<std::vector<double>>>& ionicMatrix);

        /**
         * @brief Parses a Slater determinant from a string representation.
         * 
         * @param[out] slaterDeterminant SlaterDeterminant object to populate from the string.
         * @param[in] str String representation of the Slater determinant.
         * 
         * @return True if parsing was successful, false otherwise.
         */
        static bool parseFromString(SlaterDeterminant& slaterDeterminant, const std::string& str);

        //----------------------------------------------------------------------------------------------------//
        // OPERATOR OVERLOADS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Overloads the equality operator for two SlaterDeterminant objects.
         * 
         * @param lhs Left-hand side SlaterDeterminant.
         * @param rhs Right-hand side SlaterDeterminant.
         * 
         * @return True if both Slater determinants are equal (have the same occupied orbitals and occupation numbers), false otherwise.
         */
        friend bool operator==(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs);

        /**
         * @brief Overloads the inequality operator for two SlaterDeterminant objects.
         *
         * @param lhs Left-hand side SlaterDeterminant.
         * @param rhs Right-hand side SlaterDeterminant.
         * 
         * @return True if both Slater determinants are different (have different occupied orbitals or occupation numbers), false otherwise.
         */
        friend bool operator!=(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs);

        /**
         * @brief Overloads the output stream redirection operator for a SlaterDeterminant.
         *
         * Prints the Slater determinant's occupied orbitals and occupation numbers.
         *
         * @param stream Output stream.
         * @param slaterDeterminant SlaterDeterminant to print.
         * 
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& stream, const SlaterDeterminant& slaterDeterminant);
};


#endif /* CDFTT_SLATERDETERMINANT_HPP_INCLUDED */