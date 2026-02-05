#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include "../Common/Atom.h"
#include "../Common/Constants.h"
#include "../Cube/Domain.h"


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Domain::Domain() :
    _Nval(1),
    _N1(0),
    _N2(0),
    _N3(0),
    _origin({ { 0.0, 0.0, 0.0 } }),
    _T({ { { 0.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0 } } }),
    _inv_T({ { { 1.0, 0.0, 0.0 },
               { 0.0, 1.0, 0.0 },
               { 0.0, 0.0, 1.0 } } }),
    _dx(0),
    _dy(0),
    _dz(0),
    _dv(0)
{ }

Domain::Domain(ifstream& nameFile) :
    _Nval(1),
    _N1(0),
    _N2(0),
    _N3(0),
    _origin({ { 0.0, 0.0, 0.0 } }),
    _T({ { { 0.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0 } } }),
    _inv_T({ { { 1.0, 0.0, 0.0 },
               { 0.0, 1.0, 0.0 },
               { 0.0, 0.0, 1.0 } } }),
    _dx(0),
    _dy(0),
    _dz(0),
    _dv(0)
{
    readFromCube(nameFile);
}

Domain::Domain(int Nval, int N1, int N2, int N3, double xmax, double ymax, double zmax) :
    _Nval(Nval),
    _N1(N1),
    _N2(N2),
    _N3(N3),
    _origin({ { - xmax, - ymax, - zmax } }),
    _T({ { { 2.0 * xmax / N1, 0.0, 0.0 },
           { 0.0, 2.0 * ymax / N2, 0.0 },
           { 0.0, 0.0, 2.0 * zmax / N3 } } }),
    _inv_T({{{1.0, 0.0, 0.0},
             {0.0, 1.0, 0.0},
             {0.0, 0.0, 1.0}}}),
    _dx(0.0),
    _dy(0.0),
    _dz(0.0),
    _dv(0.0)
{
    inverse_T();
    computeInfinitesimalElements();
}

Domain::Domain(int nVal, int n1, int n2, int n3, const std::array<double, 3>& origin) :
    _Nval(nVal),
    _N1(n1),
    _N2(n2),
    _N3(n3),
    _origin(origin),
    _T({ { { 0.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0 },
           { 0.0, 0.0, 0.0 } } }),
    _inv_T({ { { 1.0, 0.0, 0.0 },
               { 0.0, 1.0, 0.0 },
               { 0.0, 0.0, 1.0 } } }),
    _dx(0),
    _dy(0),
    _dz(0),
    _dv(0)
{
    inverse_T();
    computeInfinitesimalElements();
}


//----------------------------------------------------------------------------------------------------//
// PRIVATE METHODS
//----------------------------------------------------------------------------------------------------//

void Domain::computeInfinitesimalElements()
{
    for (int i = 0; i < 3; ++i)
    {
        _dx += _T[0][i] * _T[0][i];
        _dy += _T[1][i] * _T[1][i];
        _dz += _T[2][i] * _T[2][i];
    }

    _dx = std::sqrt(_dx);
    _dy = std::sqrt(_dy);
    _dz = std::sqrt(_dz);

    _dv = _dx * _dy * _dz;
}


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

int Domain::get_Nval() const
{
    return _Nval;
}

int Domain::get_N1() const
{
    return _N1;
}

int Domain::get_N2() const
{
    return _N2;
}

int Domain::get_N3() const
{
    return _N3;
}

const std::array<double, 3>& Domain::get_origin() const
{
    return _origin;
}

const std::array<std::array<double, 3>, 3>& Domain::get_T() const
{
    return _T;
}

double Domain::get_Tij(int i , int j) const
{
    return _T[i][j];
}

double Domain::get_dx() const
{
    return _dx;
}

double Domain::get_dy() const
{
    return _dy;
}

double Domain::get_dz() const
{
    return _dz;
}

double Domain::get_dv() const
{
    return _dv;
}


//----------------------------------------------------------------------------------------------------//
// SETTERS
//----------------------------------------------------------------------------------------------------//

void Domain::set_Nval(int N)
{
    _Nval = N;
}

void Domain::set_N1(int N)
{
    _N1 = N;
}

void Domain::set_N2(int N)
{
    _N2 = N;
}

void Domain::set_N3(int N)
{
    _N3 = N;
}

void Domain::set_Tij(double value, int i, int j)
{
    _T[i][j] = value;
}

void Domain::set_all(int Nval, int N1, int N2,int N3, double xmax, double ymax, double zmax, const array<array<double, 3>, 3>& T)
{
    _Nval = Nval;

    _N1 = N1;
    _N2 = N2;
    _N3 = N3;

    _origin = { { - xmax, - ymax, - zmax } };
    _T = T;

    inverse_T();
    computeInfinitesimalElements();
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

double Domain::x(int i, int j, int k) const
{
    return _origin[0] + _T[0][0] * i + _T[0][1] * j + _T[0][2] * k;
}

double Domain::y(int i, int j, int k) const
{
    return _origin[1] + _T[1][0] * i + _T[1][1] * j + _T[1][2] * k;
}

double Domain::z(int i, int j, int k) const
{
    return _origin[2] + _T[2][0] * i + _T[2][1] * j + _T[2][2] * k;
}

int Domain::i(double x, double y, double z) const
{
    return floor((x - _origin[0]) * _inv_T[0][0] + (y - _origin[1]) * _inv_T[0][1] + (z - _origin[2]) * _inv_T[0][2]);
}

int Domain::j(double x, double y, double z) const
{
    return floor((x - _origin[0]) * _inv_T[1][0] + (y - _origin[1]) * _inv_T[1][1] + (z - _origin[2]) * _inv_T[1][2]);
}

int Domain::k(double x, double y, double z) const
{
    return floor((x - _origin[0]) * _inv_T[2][0] + (y - _origin[1]) * _inv_T[2][1] + (z - _origin[2]) * _inv_T[2][2]);
}

void Domain::readFromCube(ifstream& nameFile)
{
    std::string s;
    getline(nameFile,s);
    std::stringstream ss(s);
    // Tokenize the input std::string by comma delimiter

    for(int i = 0; i < 3; ++i)
    {
        ss >> _origin[i];
    }

    _Nval = 1;
    ss >> _Nval;
    if(ss.fail())
    {
        _Nval = 1;
    }

    nameFile >> _N1;

    if(_N1 > 0)
    {
        for(int i = 0; i < 3; ++i)
        {
            nameFile >> _T[0][i];
        }

        nameFile >> _N2;

        for(int i = 0; i < 3; ++i)
        {
            nameFile >> _T[1][i];
        }

        nameFile >> _N3;

        for(int i = 0; i < 3; ++i)
        {
            nameFile >> _T[2][i];
        }

        _N1 = std::abs(_N1);
    }
    else
    {
        for(int i = 0; i < 3; ++i)
        {
            nameFile >> _T[0][i];
            _T[0][i] = _T[0][i] * Constants::ANGSTROM_TO_BOHR_RADIUS;
        }

        nameFile >> _N2;

        for(int i = 0; i < 3; ++i)
        {
            nameFile >> _T[1][i];
            _T[1][i] = _T[1][i] * Constants::ANGSTROM_TO_BOHR_RADIUS;
        }

        nameFile >> _N3;

        for(int i = 0; i < 3; ++i)
        {
            nameFile >> _T[2][i];
            _T[2][i] = _T[2][i] * Constants::ANGSTROM_TO_BOHR_RADIUS;
        }

        for(int i = 0; i < 3; ++i)
        {
            _origin[i] = _origin[i] * Constants::ANGSTROM_TO_BOHR_RADIUS;
        }
    }

    inverse_T();
    computeInfinitesimalElements();
}

void Domain::inverse_T()
{
    double t4, t6, t8, t10, t12, t14, t17;

    t4 = _T[0][0] * _T[1][1];
    t6 = _T[0][0] * _T[1][2];
    t8 = _T[0][1] * _T[1][0];
    t10 = _T[0][2] * _T[1][0];
    t12 = _T[0][1] * _T[2][0];
    t14 = _T[0][2] * _T[2][0];

    t17 = (t4 * _T[2][2] - t6 * _T[2][1] - t8 * _T[2][2] + t10 * _T[2][1] + t12 * _T[1][2] - t14 * _T[1][1]);

    if (abs(t17) < 1e-12)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                if (i == j)
                {
                    _inv_T[i][j] = 1;
                }
                else
                {
                    _inv_T[i][j] = 0;
                }
            }
        }
    }
    else
    {
        t17 = 1 / t17;

        _inv_T[0][0] = (_T[1][1] * _T[2][2] - _T[1][2] * _T[2][1]) * t17;
        _inv_T[0][1] = -(_T[0][1] * _T[2][2] - _T[0][2] * _T[2][1]) * t17;
        _inv_T[0][2] = -(-_T[0][1] * _T[1][2] + _T[0][2] * _T[1][1]) * t17;

        _inv_T[1][0] = -(_T[1][0] * _T[2][2] - _T[1][2] * _T[2][0]) * t17;
        _inv_T[1][1] = (_T[0][0] * _T[2][2] - t14) * t17;
        _inv_T[1][2] = -(t6 - t10) * t17;

        _inv_T[2][0] = -(-_T[1][0] * _T[2][1] + _T[1][1] * _T[2][0]) * t17;
        _inv_T[2][1] = -(_T[0][0] * _T[2][1] - t12) * t17;
        _inv_T[2][2] = (t4 - t8) * t17;
    }
}

double Domain::sizeUpMol(const Structure& S, double scale)
{
    double dmax = 0;

    std::vector<Atom> atoms = S.get_atoms();
    for (int i = 0; i < S.number_of_atoms(); i++)
    {
        for (int j = 0; j < S.number_of_atoms(); j++)
        {
            if (i == j)
            {
                continue;
            }

            double dtmp = atoms[i].computeDistance(atoms[j]);
            if (dtmp > dmax)
            {
                dmax = dtmp;
            }
        }
    }

    return dmax;
}


//----------------------------------------------------------------------------------------------------//
// OPERATOR OVERLOADS
//----------------------------------------------------------------------------------------------------//

bool Domain::operator==(const Domain& D) const
{
    return (D._Nval == _Nval
            && D._N1 == _N1
            && D._N2 == _N2
            && D._N3 == _N3
            && D._T == _T);
}

bool Domain::operator!=(const Domain& D) const
{
    return !(*this == D);
}


