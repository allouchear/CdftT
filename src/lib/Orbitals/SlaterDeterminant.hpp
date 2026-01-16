#ifndef CDFTT_SLATERDETERMINANT_HPP_INCLUDED
#define CDFTT_SLATERDETERMINANT_HPP_INCLUDED

#include <vector>

#include <Orbitals/Orbitals.h>


/** @brief SlaterDeterminant class.
 * 
 * Manages Slater determinants representing electronic configurations.
 */
class SlaterDeterminant
{
    private:
        /** @brief Occupied orbitals and their occupation numbers. First index corresponds to alpha spin, second to beta spin. */
        std::vector<std::vector<std::pair<int, double>>> _occupiedOrbitals;


        //----------------------------------------------------------------------------------------------------//
        // STATIC FIELDS
        //----------------------------------------------------------------------------------------------------//

        /** @brief Static Orbitals instance shared among all SlaterDeterminant objects. */
        static Orbitals _s_orbitals_;

        /** @brief Flag indicating whether the static Orbitals instance has been set. */
        static bool _s_isOrbitalsSet_;


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
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Updates the Slater determinant based on an electronic transition.
         * 
         * @param[in] initialOrbitalNumber The orbital number from which an electron is removed.
         * @param[in] initialSpin The spin type (alpha or beta) of the electron being removed.
         * @param[in] finalOrbitalNumber The orbital number to which an electron is added.
         * @param[in] finalSpin The spin type (alpha or beta) of the electron being added.
         */
        void updateFromTransition(int initialOrbitalNumber, SpinType initialSpin, int finalOrbitalNumber, SpinType finalSpin);


        //----------------------------------------------------------------------------------------------------//
        // STATIC METHODS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Computes the differences in occupied orbitals between two Slater determinants.
         * 
         * @param di First Slater determinant.
         * @param dj Second Slater determinant.
         * @return A vector containing the differences in occupied orbitals for each spin type (first index corresponds to alpha spin, second to beta spin)
         */
        static std::vector<std::vector<std::pair<int, int>>> getDifferences(const SlaterDeterminant& di, const SlaterDeterminant& dj);

        /**
         * @brief Computes the ionic potential matrix element between two Slater determinants.
         * 
         * @param[in] di First Slater determinant.
         * @param[in] dj Second Slater determinant.
         * @param[in] ionicMatrix Ionic matrix < phi_i | V_ionic | phi_j > (the first index corresponds to alpha spin, the second to beta spin).
         * @return The ionic potential matrix element < Di | V_ionic | Dj >.
         */
        static double ionicPotential(const SlaterDeterminant& di, const SlaterDeterminant& dj, const std::vector<std::vector<std::vector<double>>>& ionicMatrix);

        //----------------------------------------------------------------------------------------------------//
        // OPERATOR OVERLOADS
        //----------------------------------------------------------------------------------------------------//

        /**
         * @brief Overloads the equality operator for two SlaterDeterminant objects.
         * 
         * @param lhs Left-hand side SlaterDeterminant.
         * @param rhs Right-hand side SlaterDeterminant.
         * @return True if both Slater determinants are equal (have the same occupied orbitals and occupation numbers), false otherwise.
         */
        friend bool operator==(const SlaterDeterminant& lhs, const SlaterDeterminant& rhs);

        /**
         * @brief Overloads the inequality operator for two SlaterDeterminant objects.
         *
         * @param lhs Left-hand side SlaterDeterminant.
         * @param rhs Right-hand side SlaterDeterminant.
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
         * @return Reference to the output stream.
         */
        friend std::ostream& operator<<(std::ostream& stream, const SlaterDeterminant& slaterDeterminant);
};


#endif /* CDFTT_SLATERDETERMINANT_HPP_INCLUDED */