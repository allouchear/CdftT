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
        /** @brief Occupied orbitals and their occupation numbers for the alpha spin. */
        std::vector<std::pair<int, double>> _alphaOccupation;

        /** @brief Occupied orbitals and their occupation numbers for the beta spin. */
        std::vector<std::pair<int, double>> _betaOccupation;

        
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
         * @brief Computes the ionic potential matrix element between two Slater determinants.
         * 
         * @param di First Slater determinant.
         * @param dj Second Slater determinant.
         * @param ionicMatrix Precomputed ionic potential matrix.
         * @return The ionic potential matrix element < Di | V_ion | Dj >.
         */
        static double ionicPotentialSlaterDeterminant(const SlaterDeterminant& di, const SlaterDeterminant& dj, const std::vector<std::vector<double>>& ionicMatrix);

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