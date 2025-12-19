#ifndef CDFTT_UTILS_H_INCLUDED
#define CDFTT_UTILS_H_INCLUDED

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>


//----------------------------------------------------------------------------------------------------//
// STRING MANAGEMENT FUNCTIONS
//----------------------------------------------------------------------------------------------------//

std::string to_lower(const std::string &str);
std::string to_upper(const std::string &str);
std::string trim_whitespaces(const std::string &str, bool leading, bool trailing);


//----------------------------------------------------------------------------------------------------//
// PRINT FUNCTIONS
//----------------------------------------------------------------------------------------------------//

void print_title(const std::string& title);
void print_error(const std::string& errorMessage);


//----------------------------------------------------------------------------------------------------//
// FILE PARSING FUNCTIONS
//----------------------------------------------------------------------------------------------------//

/** @brief Searches for the first occurrence of the specified parameter (eg: "RunType") in the input file and reads its value as a string.
 * 
 * The position of the input file stream is reset to the beginning of the file before searching.
 * The resulting value is trimmed of leading and trailing whitespaces.
 * 
 * @param[in] inputFile Input file stream to read from.
 * @param[in] parameter Parameter to search for.
 * @param[out] value Reference to a string where the read value will be stored.
 * @return True if the parameter was found and read successfully; false otherwise.
 */
bool readOneString(std::ifstream& inputFile, const std::string& parameter, std::string& value);

/** @brief Searches for the first occurrence of the specified parameter (eg: "RunType") in the input file and reads its value as the specified type.
 *
 * The position of the input file stream is reset to the beginning of the file before searching.
 *
 * @tparam T Type of the value to read (e.g., int, double, ...).
 * @param[in] inputFile Input file stream to read from.
 * @param[in] parameter Parameter to search for.
 * @param[out] value Reference to a type T variable where the read value will be stored.
 * @return True if the parameter was found and read successfully; false otherwise.
 */
template<typename T> bool readOneType(std::ifstream& inputFile, const std::string& tag, T& x);

/** @brief Searches for the first occurrence of the specified parameter (eg: "RunType") in the input file and reads its value as a list of the specified type.
 *
 * The position of the input file stream is reset to the beginning of the file before searching.
 * Accepted separators between list elements are spaces, commas, or semicolons.
 *
 * @tparam T Type of the values to read (e.g., int, double, ...).
 * @param[in] inputFile Input file stream to read from.
 * @param[in] parameter Parameter to search for.
 * @param[out] value Reference to a vector of type T where the read values will be stored.
 * @return True if the parameter was found and read successfully; false otherwise.
 */
template<typename T> bool readListType(std::ifstream& inputFile, const std::string& tag, std::vector<T>& x);

#include <Utils/Utils_fileParsing.tpp>


//----------------------------------------------------------------------------------------------------//
// MATRIX MANAGEMENT FUNCTIONS
//----------------------------------------------------------------------------------------------------//

bool diagonalisationOfATridiagonalMatrix(std::vector<double>& subDiagonal, std::vector<double>& eigenValues, std::vector<std::vector<double>>& eigenVectors);

/**
 * @brief Finds the eigenvalues and eigenvectors of a real symmetric matrix using the QL algorithm.
 * 
 * @param[in] matrixLowerTriangle Lower triangle of a real, symmetric matrix.
 * @param[out] eigenValues Vector where the computed eigenvalues will be stored.
 * @param[out] eigenVectors Matrix (as a vector of vectors) where the computed eigenvectors will be stored.
 * @return True if the computation was successful; false otherwise.
 */
bool findEigenValuesAndEigenVectorsOfSymmetricalMatrix(const std::vector<std::vector<double>>& matrixLowerTriangle, std::vector<double>& eigenValues, std::vector<std::vector<double>>& eigenVectors);

/**
 * @brief Reduces a real, symmetric matrix to tridiagonal form using Householder transformations.
 * 
 * @param[in,out] matrix Real symmetric matrix to be reduced (modified in place).
 * @param[out] diagonal Vector where the diagonal elements of the tridiagonal matrix will be stored.
 * @param[out] subDiagonal Vector where the sub-diagonal elements of the tridiagonal matrix will be stored.
 */
void reductionToTridiagonalMatrix(std::vector<std::vector<double>>& matrix, std::vector<double>& diagonal, std::vector<double>& subDiagonal);


//----------------------------------------------------------------------------------------------------//
// UTIL CLASSES
//----------------------------------------------------------------------------------------------------//

class CustomSizeData
{
    private:
        std::vector<double> _customSizeData;
    
    
    public:
        CustomSizeData();
        CustomSizeData(const std::vector<double>& customSizeData);

        double nx() const;

        double ny() const;

        double nz() const;

        double ox() const;

        double oy() const;

        double oz() const;

        double t11() const;

        double t12() const;

        double t13() const;

        double t21() const;

        double t22() const;

        double t23() const;

        double t31() const;

        double t32() const;

        double t33() const;

        double operator[](size_t i) const;
};


//! A factorial class.
/*! This class will be used in the GTF class for the normalisation. */
class Factorial
{
    private:
        std::vector<double> _tab;
    
    
    public:

            //! A default constructor.
            /*! This constructor create a table without a size (0). */
        Factorial();

            //! A real constructor.
            /*! This constructor is used to create a table from 0 to n factorial. */
        Factorial(int);

            //! A default desctructor.
            /*! We don't use it. */
        ~Factorial(){}

            //! A normal member taking one argument and returning a double value.
            /*! \return The n factorial value. */
        double factorial(int);

            //! A normal member taking one argument and returning a double value.
            /*! \return The n double factorial value. */
        double double_factorial(int);
};

//! A binomial class.
/*! This class will be used in the GTF class for some calculus. */
class Binomial
{
    private:
        std::vector<std::vector<double>> _tab;
        Factorial _fact;

    public:
        //! A real constructor.
        /*! This constructor is used to create a table from 0 to n binomial. */
        Binomial(int, Factorial &);

        //! A default constructor.
        /*! This constructor create a table without a size (0). */
        Binomial();

        //! A default desctructor.
        /*! We don't use it. */
        ~Binomial() {}

        //! A normal member taking two arguments and returning a double value.
        /*! \return The i,j binomial value. */
        double binomial(int, int);

        //! A normal member taking no arguments and returning a vector<vector<double>> value.
        /*! \return The binomial table. */
        std::vector<std::vector<double>> &tab()
        {
            return _tab;
        }

        //! A normal member taking no arguments and returning a Factorial value.
        /*! \return The factorial table. */
        Factorial& fact()
        {
            return _fact;
        }
};

    //! A method for approximation power (small number) taking two arguments and returning a double value.
    /*! \return The result */
double power(double, int);

    //! A method taking six arguments and returning a double value.
    /*! \return ???*/
double f(int, int, int, double, double, Binomial&);

    //! A method taking eight arguments and returning a double value.
    /*! \return ???*/
double Theta(int, int, int, int, double, double, double, Binomial&);

    //! A method taking one argument and returning an int value.
    /*! \return ???*/
int m1p(int);

    //! A method taking ten arguments and returning a double value.
    /*! \return ???*/
double A(int, int, int, int, int, double, double, double, double, Binomial&);

    //! A method taking ten arguments and returning a double value.
    /*! \return ???*/
double B(int,int,int , int , int , double , double , double , double , Factorial&);

    //! A method taking two arguments and returning a double value.
    /*! \return ???*/
double myGamma(int, Factorial&);

    //! A method taking three arguments and returning a double value.
    /*! \return ???*/
double F(int,double, Factorial&);

    //! A method taking three arguments and returning a vector<double> value.
    /*! \return ???*/
std::vector<double> getFTable(int, double, Factorial&);

    //! A method taking one argument and returning an int value.
    /*! \return The L type (S, Px, Py, ...) in the nomenclatur of .wfx file.*/
int getwfxType(const std::vector<int>& l);

    //! A method taking one argument and returning an int value.
    /*! \return The shell type (S, P, D, ...).*/
std::string getLType(const std::vector<int>& l);

#endif