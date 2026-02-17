#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// STRING MANAGEMENT FUNCTIONS
//----------------------------------------------------------------------------------------------------//

std::string int_to_string_withLeadingZeros(const int value, const int maxValue)
{
    // Count the number of digits in max
    int nbDigits = std::to_string(std::abs(maxValue)).length();
    
    // Format value with leading zeros
    std::stringstream sstreamWithLeadingZeros;
    sstreamWithLeadingZeros << std::setfill('0') << std::setw(nbDigits) << value;
    
    return sstreamWithLeadingZeros.str();
}

std::string to_lower(const std::string& str)
{
    std::string res = str;

    for (size_t i = 0; i < res.length(); i++)
    {
        res[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(res[i])));
    }

    return res;
}

std::string to_upper(const std::string& str)
{
    std::string res = str;

    for (size_t i = 0; i < res.length(); i++)
    {
        res[i] = static_cast<char>(std::toupper(static_cast<unsigned char>(res[i])));
    }

    return res;
}

std::string trim_whitespaces(const std::string& str, bool leading, bool trailing)
{
    std::string res = str;

    if (leading)
    {
        res = std::regex_replace(res, std::regex("^ +"), "");
    }

    if (trailing)
    {
        res = std::regex_replace(res, std::regex(" +$"), "");
    }

    return res;
}


//----------------------------------------------------------------------------------------------------//
// PRINT FUNCTIONS
//----------------------------------------------------------------------------------------------------//

void log(std::stringstream& messageStream, std::ostream& outputStream)
{
    std::vector<std::reference_wrapper<std::ostream>> outputStreams;

    outputStreams.emplace_back(outputStream);
    if (&outputStream != &std::cout)
    {
        outputStreams.emplace_back(std::cout);
    }

    for (auto& output : outputStreams)
    {
        if (output.get())
        {
            output.get() << messageStream.str();
        }
    }

    messageStream.clear();
    messageStream.str(std::string());
}

void print_title(const std::string& title)
{
    std::stringstream sstreamTitle;
    sstreamTitle << "| " << title << " |";

    std::string titleFullStr = sstreamTitle.str();
    int nbDashes = (int)(titleFullStr.length()) - 2;

    std::stringstream sstreamDashes;
    for (int i = 0; i < nbDashes; ++i)
    {
        sstreamDashes << "-";
    }
    std::string dashes = sstreamDashes.str();

    std::cout << '/' << dashes << '\\' << std::endl;
    std::cout << titleFullStr << std::endl;
    std::cout << '\\' << dashes << '/' << std::endl << std::endl;
}

void print_error(const std::string& errorMessage, std::ostream& outputStream)
{
    std::vector<std::reference_wrapper<std::ostream>> outputStreams;
    
    outputStreams.emplace_back(outputStream);
    if (&outputStream != &std::cerr)
    {
        outputStreams.emplace_back(std::cerr);
    }

    for (auto& output : outputStreams)
    {
        if (output.get())
        {
            output.get() << "!!!!!!!!!" << std::endl;
            output.get() << "! ERROR !" << std::endl;
            output.get() << "!!!!!!!!!" << std::endl << std::endl;

            output.get() << errorMessage << std::endl;
        }
    }
}


//----------------------------------------------------------------------------------------------------//
// FILE PARSING FUNCTIONS
//----------------------------------------------------------------------------------------------------//

bool readOneString(std::ifstream& inputFile, const std::string& tag, std::string& value)
{
    bool ok = false;

    if (tag.length() >= 1)
    {
        std::string TAG = to_upper(tag);

        inputFile.clear();
        inputFile.seekg(0);

        std::string t;
        std::string t2;

        value = "";

        while (!ok && !inputFile.eof())
        {
            std::getline(inputFile, t);
            if (inputFile.fail())
            {
                break;
            }

            // Trimming leading and trailing spaces
            t = trim_whitespaces(t, true, true);

            // Ignoring commment lines
            if (t[0] == '#')
            {
                continue;
            }

            t2 = to_upper(t);
            std::size_t pos = t2.find(TAG);

            // If tag not found, continue to next line
            if (pos == std::string::npos)
            {
                continue;
            }

            // Tag found, extracting value
            if (t2.find("=") != std::string::npos)
            {
                pos = t2.find("=");
                value = t.substr(pos + 1);
            }
            else // if no '=' is found, we assume a space separates tag and value
            {
                pos = t2.find(" ");
                value = t.substr(pos + 1);
            }

            if (value.length() > 0)
            {
                value = trim_whitespaces(value, true, true);
                ok = true;
            }
        }
    }

    return ok;
}


//----------------------------------------------------------------------------------------------------//
// MATRIX MANAGEMENT FUNCTIONS
//----------------------------------------------------------------------------------------------------//

bool diagonalisationOfATridiagonalMatrix(std::vector<double>& subDiagonal, std::vector<double>& eigenValues, std::vector<std::vector<double>>& eigenVectors)
{
    bool ok = true;
    int dimension = static_cast<int>(eigenValues.size());

    for (int i = 1; i < dimension; ++i)
    {
        subDiagonal[i - 1] = subDiagonal[i];
    }
    subDiagonal[dimension - 1] = 0.0;

    //for (int l = 0; ok && l < dimension; ++l)
    for (int l = 0; l < dimension; ++l)
    {
        int iteration = 0;

        int m;
        do
        {
            for (m = l; m < dimension - 1; ++m)
            {
                double dd = std::fabs(eigenValues[m]) + std::fabs(eigenValues[m + 1]);
                if (std::fabs(subDiagonal[m]) + dd == dd)
                {
                    break;
                }
            }

            if (m != l)
            {
                ++iteration;
                if (iteration == 30)
                {
                    ok = false;
                    break;
                }

                double g = (eigenValues[l + 1] - eigenValues[l]) / (2.0 * subDiagonal[l]);
                double r = std::sqrt((g * g) + 1.0);

                g = eigenValues[m] - eigenValues[l] + subDiagonal[l] / (g + (g < 0.0 ? -std::fabs(r) : std::fabs(r)));

                double s = 1.0;
                double c = 1.0;
                double p = 0.0;
                for (int i = m - 1; i >= l; --i)
                {
                    double f = s * subDiagonal[i];
                    double b = c * subDiagonal[i];

                    if (std::fabs(f) >= std::fabs(g))
                    {
                        c = g / f;
                        r = std::sqrt((c * c) + 1.0);

                        subDiagonal[i + 1] = f * r;

                        s = 1.0 / r;
                        c *= s;
                    }
                    else
                    {
                        s = f / g;
                        r = std::sqrt((s * s) + 1.0);

                        subDiagonal[i + 1] = g * r;

                        c = 1.0 / r;
                        s *= c;
                    }

                    g = eigenValues[i + 1] - p;
                    r = (eigenValues[i] - g) * s + 2.0 * c * b;
                    p = s * r;

                    eigenValues[i + 1] = g + p;
                    g = c * r - b;

                    for (int k = 0; k < dimension; ++k)
                    {
                        f = eigenVectors[k][i + 1];
                        eigenVectors[k][i + 1] = s * eigenVectors[k][i] + c * f;
                        eigenVectors[k][i] = c * eigenVectors[k][i] - s * f;
                    }
                }

                eigenValues[l] = eigenValues[l] - p;
                subDiagonal[l] = g;
                subDiagonal[m] = 0.0;
            }
        } while (m != l);
    }

    return ok;
}

bool findEigenValuesAndEigenVectorsOfSymmetricalMatrix(const std::vector<std::vector<double>>& matrixLowerTriangle, std::vector<double>& eigenValues, std::vector<std::vector<double>>& eigenVectors)
{
    bool ok = false;

    int dimension = static_cast<int>(matrixLowerTriangle.size());
    if (dimension > 1)
    {
        ok = true;

        // Building whole matrix from lower triangle
        eigenVectors = matrixLowerTriangle;
        for (int i = 0; i < dimension; ++i)
        {
            eigenVectors[i].resize(dimension);

            for (int j = i + 1; j < dimension; ++j)
            {
                eigenVectors[i][j] = eigenVectors[j][i];
            }
        }

        std::vector<double> subDiagonal;
        reductionToTridiagonalMatrix(eigenVectors, eigenValues, subDiagonal);

        ok = diagonalisationOfATridiagonalMatrix(subDiagonal, eigenValues, eigenVectors);
    }

    return ok;
}

void reductionToTridiagonalMatrix(std::vector<std::vector<double>>& matrix, std::vector<double>& diagonal, std::vector<double>& subDiagonal)
{
    // Ensure diagonal and subDiagonal have the correct size
    int dimension = static_cast<int>(matrix.size());
    diagonal.resize(dimension);
    subDiagonal.resize(dimension);

    for (int i = dimension - 1; i >= 1; --i)
    {
        int l = i - 1;

        double h = 0.0;
        double scale = 0.0;

        if (l > 0)
        {
            for (int k = 0; k <= l; ++k)
            {
                scale += std::fabs(matrix[i][k]);
            }

            if (scale == 0.0)
            {
                subDiagonal[i] = matrix[i][l];
            }
            else
            {
                for (int k = 0; k <= l; ++k)
                {
                    matrix[i][k] /= scale;
                    h += matrix[i][k] * matrix[i][k];
                }

                double f = matrix[i][l];
                double g = (f > 0.0) ? -std::sqrt(h) : std::sqrt(h);

                subDiagonal[i] = scale * g;
                h -= f * g;
                matrix[i][l] = f - g;
                
                f = 0.0;
                for (int j = 0; j <= l; ++j)
                {
                    matrix[j][i] = matrix[i][j] / h;

                    g = 0.0;
                    for (int k = 0; k <= j; ++k)
                    {
                        g += matrix[j][k] * matrix[i][k];
                    }
                    for (int k = j + 1; k <= l; ++k)
                    {
                        g += matrix[k][j] * matrix[i][k];
                    }

                    subDiagonal[j] = g / h;
                    f += subDiagonal[j] * matrix[i][j];
                }

                double hh = f / (h + h);
                for (int j = 0; j <= l; ++j)
                {
                    f = matrix[i][j];

                    g = subDiagonal[j] - hh * f;
                    subDiagonal[j] = g;
                    for (int k = 0; k <= j; ++k)
                    {
                        matrix[j][k] -= (f * subDiagonal[k] + g * matrix[i][k]);
                    }
                }
            }
        }
        else
        {
            subDiagonal[i] = matrix[i][l];
        }

        diagonal[i] = h;
    }

    diagonal[0] = 0.0;
    subDiagonal[0] = 0.0;
    for (int i = 0; i < dimension; ++i)
    {
        int l = i - 1;

        if (diagonal[i] != 0.0)
        {
            for (int j = 0; j <= l; ++j)
            {
                double g = 0.0;
                for (int k = 0; k <= l; ++k)
                {
                    g += matrix[i][k] * matrix[k][j];
                }
                for (int k = 0; k <= l; ++k)
                {
                    matrix[k][j] -= g * matrix[k][i];
                }
            }
        }

        diagonal[i] = matrix[i][i];
        matrix[i][i] = 1.0;
        for (int j = 0; j <= l; ++j)
        {
            matrix[j][i] = 0.0;
            matrix[i][j] = 0.0;
        }
    }
}

void sortEigenValuesAndEigenVectors(std::vector<double>& eigenValues, std::vector<std::vector<double>>& eigenVectors)
{
    size_t dimension = eigenValues.size();

    if (eigenVectors.size() != dimension)
    {
        print_error("Error in sortEigenValuesAndEigenVectors(): eigenVectors size does not match eigenValues size.");
        exit(1);
    }

    // Create a vector of indices
    std::vector<size_t> indices(dimension);
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices based on eigenvalues
    std::sort(indices.begin(), indices.end(), [&eigenValues](size_t a, size_t b) { return eigenValues[a] < eigenValues[b]; });

    // Create sorted copies
    std::vector<double> sortedEigenValues(dimension);
    std::vector<std::vector<double>> sortedEigenVectors(eigenVectors);

    // Reorganize eigenvalues and eigenvectors according to sorted indices
    for (size_t i = 0; i < dimension; ++i)
    {
        sortedEigenValues[i] = eigenValues[indices[i]];
        
        // Copy column indices[i] to column i (eigenvectors stored as columns)
        for (size_t k = 0; k < dimension; ++k)
        {
            sortedEigenVectors[k][i] = eigenVectors[k][indices[i]];
        }
    }

    // Replace with sorted data
    eigenValues = std::move(sortedEigenValues);
    eigenVectors = std::move(sortedEigenVectors);
}

//----------------------------------------------------------------------------------------------------//
// UTIL CLASSES
//----------------------------------------------------------------------------------------------------//

CustomSizeData::CustomSizeData()
{ }

CustomSizeData::CustomSizeData(const std::vector<double> &customSizeData) : _customSizeData(customSizeData)
{ }

double CustomSizeData::nx() const
{
    return _customSizeData[0];
}

double CustomSizeData::ny() const
{
    return _customSizeData[1];
}

double CustomSizeData::nz() const
{
    return _customSizeData[2];
}

double CustomSizeData::ox() const
{
    return _customSizeData[3];
}

double CustomSizeData::oy() const
{
    return _customSizeData[4];
}

double CustomSizeData::oz() const
{
    return _customSizeData[5];
}

double CustomSizeData::t11() const
{
    return _customSizeData[6];
}

double CustomSizeData::t12() const
{
    return _customSizeData[7];
}

double CustomSizeData::t13() const
{
    return _customSizeData[8];
}
double CustomSizeData::t21() const
{
    return _customSizeData[9];
}

double CustomSizeData::t22() const
{
    return _customSizeData[10];
}

double CustomSizeData::t23() const
{
    return _customSizeData[11];
}
double CustomSizeData::t31() const
{
    return _customSizeData[12];
}

double CustomSizeData::t32() const
{
    return _customSizeData[13];
}

double CustomSizeData::t33() const
{
    return _customSizeData[14];
}

double CustomSizeData::operator[](size_t i) const
{
    return _customSizeData[i];
}


Factorial::Factorial()
{
    _tab = std::vector<double>();
}

Factorial::Factorial(int n)
{
    _tab = std::vector<double>(n, 1);

    for (int i = 2; i < n; ++i)
    {
        for (int k = 0; k <= i / 2 - 1; ++k)
        {
            _tab[i] *= i - 2 * k;
        }
    }
}

double Factorial::factorial(int n)
{
    return n == 0 ? 1.0 : double_factorial(n) * double_factorial(n - 1);
}

double Factorial::double_factorial(int n)
{
    if (n<0)
        return 1;
    
    if (size_t(n) >= _tab.size())
    {
        double r;
        for (size_t i = _tab.size(); i < size_t(n); i++)
        {    
            r = 1;
            for (size_t k = 0; k <= i / 2 - 1; k++)
            {
                r *= i - 2 * k;
            }

            _tab.push_back(r);
        }
    }

    return _tab[n];
}

Binomial::Binomial()
{
    _fact = Factorial();
    _tab = std::vector<std::vector<double>>();
}

Binomial::Binomial(int i, Factorial& F) : _fact(F)
{
    std::vector<double> V(0);
    _tab = std::vector<std::vector<double>>(i, V);
    for (size_t k = 0; k < _tab.size(); k++)
    {
        _tab[k].resize(k+1);
        for (size_t l = 0; l <= k; l++)
        {
            _tab[k][l] = _fact.factorial(k) / _fact.factorial(l) / _fact.factorial(k - l);
        }
    }
}

double Binomial::binomial(int i, int j)
{
/*                                                    A debug
    if (i>_tab.size() || j>_tab[0].size())
    {
        cout<<"Test 1"<<endl;
        _tab=Binomial(i, j ,F).tab();
        cout<<"Test 2"<<endl;
    }
*/
    return _tab[i][j];
}

double power(double e, int n)
{
    double p = 1.0;
    if (std::fabs(e) < 1e-10)
    {
        p = (n == 0 ? 1.0 : 0.0);
    }
    else
    {
        for (int k = 1; k <= n; ++k)
        {
            p *= e;
        }
    }
    
    return p;
}

double f(int i, int l, int m, double A, double B, Binomial& binomial)
{
    double sum = 0.0;

    int jmin = (0 > i - m ? 0 : i - m);
    int jmax = (i < l ? i : l);

    for (int j = jmin; j <= jmax; ++j)
    {
        sum += binomial.binomial(l, j) * binomial.binomial(m, i - j) * power(-A, l - j) * power(-B, m - i + j);
    }

    return sum;
}

double Theta(int i,int r,int l1,int l2, double A, double B, double g, Binomial& Bi)
{
    return f(i,l1,l2,A,B,Bi)*Bi.fact().factorial(i)/Bi.fact().factorial(r)/Bi.fact().factorial(i-2*r)/pow(g,i-r);
}

int m1p(int i)
{
    return (i % 2 == 0) ? 1 : -1;
}

double A(int i, int r, int u, int l1, int l2, double A, double B, double C, double gamma, Binomial& binomial)
{
    return m1p(i + u) * f(i, l1, l2, A, B, binomial)
                      * binomial.fact().factorial(i)
                      * power(C, i - 2 * (r + u))
                      / (binomial.fact().factorial(r)
                         * binomial.fact().factorial(u)
                         * binomial.fact().factorial(i - 2 * r - 2 * u)
                         * power(4 * gamma, r + u));
}

double B(int i, int ip, int r, int rp, int u, double PQ, double d, double T1, double T2, Factorial& Fa)
{
    int ii = i + ip - 2 * r - 2 * rp;
    return m1p(ip + u) * T1 * T2 * Fa.factorial(ii)
                                 / Fa.factorial(u)
                                 / Fa.factorial(ii - 2 * u)
                                 * pow(PQ, ii - 2 * u)
                                 / (pow(4.0, i + ip - r - rp)
                                 * pow(d, ii - u));
}

/*
Old functions used to obtain F_nu(t) but the double factorial in the denominator of x variable in the while loop was diverging and we obtained infinite/nan values in the program.

double myGamma(int n, Factorial& Fa)
{
    return Fa.double_factorial(2 * n - 1) * sqrt(M_PI) / power(2, n);
}

double F(int n, double t, Factorial& factorial)
{
    double et = std::exp(-t);
    double twot = 2 * t;
    double T = 0.0;
    double x = 1.0;
    int i = 0;
    double DD = 1.0;
    double TMAX = 50.0;
    int MAXFACT = 200;
    double acc = 1e-16;

    if (std::fabs(t) <= acc)
    {
        T = 1.0 / (2 * n + 1);
    }
    else
    {
        if (t >= TMAX)
        {
            T = myGamma(n, factorial) / power(t, n) / 2.0 / std::sqrt(t);
        }
        else
        {
            while(std::fabs(x / T) > acc && (n + i) < MAXFACT)
            {
                x = factorial.double_factorial(2 * n - 1) / factorial.double_factorial(2 * (n + i + 1) - 1) * DD;
                T += x;
                i++;
                DD *= twot;
            }

            if (n + i >= MAXFACT)
            {
                std::cout << "Divergence in F, Ionic integrals" << std::endl;
                exit(1);
            }

            T *= et;
        }
    }

    return T;
}

std::vector<double> getFTable(int mMax, double t, Factorial& factorial)
{
    double tCritic = 30.0;

    std::vector<double> Fmt(mMax + 1);

    if (t > tCritic)
    {
        Fmt[0] = std::sqrt(M_PI / t) * 0.5;
        for (int m = 1; m <= mMax; ++m)
        {
            Fmt[m] = Fmt[m - 1] * (m - 0.5) / t;
        }
    }
    else
    {
        Fmt[mMax] = F(mMax, t, factorial);
        double expt = std::exp(-t);
        double twot = 2 * t;
        for (int m = mMax - 1; m >= 0; --m)
        {
            Fmt[m] = (twot * Fmt[m + 1] + expt) / (m * 2 + 1);
        }
    }

    return Fmt;
}
*/

double F(int nu, double t)
{
    // Exact t = 0 result
    double F = 1.0 / (2 * nu + 1);

    if (std::abs(t) > 1e-10)
    {
        // Initial value using error function for nu = 0
        F = 0.5 * std::sqrt(M_PI / t) * std::erf(std::sqrt(t));

        double expMinust = std::exp(-t);
        double twot = 2.0 * t;

        // Apply recurrence
        for (int n = 1; n <= nu; ++n)
        {
            F = ((2.0 * n + 1.0) * F - expMinust) / twot;
        }
    }

    return F;
}

std::vector<double> getFTable(int nuMax, double t, bool debug)
{
    // Check parameter validity
    if (nuMax < 0)
    {
        print_error("Error in getFTable(int nuMax, double t): nuMax must be non-negative.");
        exit(1);
    }

    if (t < 0.0)
    {
        print_error("Error in getFTable(int nuMax, double t): t must be non-negative.");
        exit(1);
    }

    if (debug)
    {
        std::cout << "t = " << t << std::endl;
    }

    std::vector<double> FTable(nuMax + 1);

    // Exact t = 0 result
    if (std::abs(t) < 1e-10)
    {
        for (int n = 0; n <= nuMax; ++n)
        {
            FTable[n] = 1.0 / (2 * n + 1);
        }
    }
    else
    {
        // Initial value using error function
        FTable[0] = 0.5 * std::sqrt(M_PI / t) * std::erf(std::sqrt(t));

        double expMinust = std::exp(-t);
        double twot = 2.0 * t;

        if (debug)
        {
            std::cout << "exp(-t) = " << expMinust << std::endl;
        }

        // Apply recurrence
        for (int n = 0; n < nuMax; ++n)
        {
            FTable[n + 1] = ((2.0 * n + 1.0) * FTable[n] - expMinust) / twot;
        }
    }

    return FTable;
}

int getwfxType(const std::vector<int>& l)
{
    int iType = 0;
    int shellType = l[0] + l[1] + l[2];

    if (shellType == 0) // S
    {
        iType = 1; // S
    }
    else if (shellType == 1) // P 
    {
        if (l[0] == 1)
        {
            iType = 2; // Px
        }
        else if (l[1] == 1)
        {
            iType = 3; // Py
        }
        else
        {
            iType = 4; // Pz
        }
    }
    else if (shellType == 2) // D
    {
        if (l[0] == 2)
        {
            iType = 5; // Dxx
        }
        else if (l[1] == 2)
        {
            iType = 6; // Dyy
        }
        else if (l[2] == 2)
        {
            iType = 7; // Dzz
        }
        else if (l[0] == 1 && l[1] == 1)
        {
            iType = 8; // Dxy
        }
        else if (l[0] == 1 && l[2] == 1)
        {
            iType = 9; // Dxz
        }
        else
        {
            iType = 10; // Dyz
        }
    }
    else if (shellType == 3) // F 
    {
        if (l[0] == 3)
        {
            iType = 11; // Fxxx
        }
        else if (l[1] == 3)
        {
            iType = 12; // Fyyy
        }
        else if (l[2] == 3)
        {
            iType = 13; // Fzzz
        }
        else if (l[0] == 2 && l[1] == 1)
        {
            iType = 14; // Fxxy
        }
        else if (l[0] == 2 && l[2] == 1)
        {
            iType = 15; // Fxxz
        }
        else if (l[1] == 2 && l[2] == 1)
        {
            iType = 16; // Fyyz
        }
        else if (l[0] == 1 && l[1] == 2)
        {
            iType = 17; // Fxyy
        }
        else if (l[0] == 1 && l[2] == 2)
        {
            iType = 18; // Fxzz
        }
        else if (l[1] == 1 && l[2] == 2)
        {
            iType = 19; // Fyzz
        }
        else
        {
            iType = 20; // Fxyz
        }
    }
    else if (shellType == 4) // G
    {
        if (l[0] == 4)
        {
            iType = 21; // Gxxxx
        }
        else if (l[1] == 4)
        {
            iType = 22; // Gyyyy
        }
        else if (l[2] == 4)
        {
            iType = 23; // Gzzzz
        }
        else if (l[0] == 3 && l[1] == 1)
        {
            iType = 24; // Gxxxy
        }
        else if (l[0] == 3 && l[2] == 1)
        {
            iType = 25; // Gxxxz
        }
        else if (l[0] == 1 && l[1] == 3)
        {
            iType = 26; // Gxyyy
        }
        else if (l[1] == 3 && l[2] == 1)
        {
            iType = 27; // Gyyyz
        }
        else if (l[0] == 1 && l[2] == 3)
        {
            iType = 28; // Gxzzz
        }
        else if (l[1] == 1 && l[2] == 3)
        {
            iType = 29; // Gyzzz
        }
        else if (l[0] == 2 && l[1] == 2)
        {
            iType = 30; // Gxxyy
        }
        else if (l[0] == 0 && l[2] == 2)
        {
            iType = 31; // Gxxzz
        }
        else if (l[1] == 2 && l[2] == 2)
        {
            iType = 32; // Gyyzz
        }
        else if (l[0] == 2 && l[1] == 1 && l[2] == 1)
        {
            iType = 33; // Gxxyz
        }
        else if (l[0] == 1 && l[1] == 2 && l[2] == 1)
        {
            iType = 34; // Gxyyz
        }
        else
        {
            iType = 35; // Gxyzz
        }
    }
    else // H and more
    {
        iType = 35;
        int L, ix, iy;

        for (L = 5; L <= 30; ++L)
        {
            for (ix = 0; ix < L; ++ix)
            {
                for (iy = 0; iy <= L - ix; ++iy)
                {
                    iType++;

                    if (l[0] == ix && l[1] == iy && l[2] == L - ix - iy)
                    {
                        return iType;
                    }
                }
            }
        }
    }

    return iType;
}

std::string getLType(const std::vector<int>& l)
{
    std::string LType = "None";
    int shellType = l[0] + l[1] + l[2];

    if (shellType == 0) // S
    {
        LType = "S"; // S
    }
    else if (shellType == 1) // P 
    {
        LType = "P";
    }
    else if (shellType == 2) // D
    {
        LType = "D";
    }
    else if (shellType == 3) // F 
    {
        LType = "F";
    }
    else if (shellType == 4) // G
    {
        LType = "G";
    }
    else // H and more
    {
        LType = std::to_string(shellType + int('H') - 5);
    }

    return LType;
}
