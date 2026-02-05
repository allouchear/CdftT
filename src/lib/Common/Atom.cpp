#include <array>
#include <cmath>
#include <iostream>
#include <string>

#include "../Common/Atom.h"
#include "../Common/Constants.h"


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Atom::Atom():
    _coordinates({ 0.0, 0.0, 0.0 }),
    _gradient({ 0.0, 0.0, 0.0 }),
    _velocity({ 0.0, 0.0, 0.0 }),
    _name("none"),
    _symbol("none"),
    _atomicNumber(0),
    _charge(0.0),
    _charge0(0.0),
    _hardness(0.0),
    _width(0.0),
    _covalentRadius(0.0),
    _element(),
    _mmType("none"),
    _pdbType("none"),
    _residueName("none"),
    _residueNumber(0),
    _N(0)
{ }

Atom::Atom(const PeriodicTable& periodicTable, const std::string& name):
    _coordinates({ 0.0, 0.0, 0.0 }),
    _gradient({ 0.0, 0.0, 0.0 }),
    _velocity({ 0.0, 0.0, 0.0 }),
    _name("none"),
    _symbol("none"),
    _atomicNumber(0),
    _charge(0.0),
    _charge0(0.0),
    _hardness(0.0),
    _width(0.0),
    _covalentRadius(0.0),
    _element(),
    _mmType("none"),
    _pdbType("none"),
    _residueName("none"),
    _residueNumber(0),
    _N(0)
{
    _element = periodicTable.element(name);
    _name = _element.get_name();
    _symbol = _element.get_symbol();
    _atomicNumber = _element.get_atomicNumber();
    _covalentRadius=_element.get_covalentRadius();
}

Atom::Atom(const PeriodicTable& periodicTable, const int Z):
    _coordinates({ 0.0, 0.0, 0.0 }),
    _gradient({ 0.0, 0.0, 0.0 }),
    _velocity({ 0.0, 0.0, 0.0 }),
    _name("none"),
    _symbol("none"),
    _atomicNumber(0),
    _charge(0.0),
    _charge0(0.0),
    _hardness(0.0),
    _width(0.0),
    _covalentRadius(0.0),
    _element(),
    _mmType("none"),
    _pdbType("none"),
    _residueName("none"),
    _residueNumber(0),
    _N(0)
{
    _element = periodicTable.element(Z);
    _name = _element.get_name();
    _symbol = _element.get_symbol();
    _atomicNumber = _element.get_atomicNumber();
    _covalentRadius=_element.get_covalentRadius();
}


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

const std::array<double, 3>& Atom::get_coordinates() const
{
    return _coordinates;
}

const std::array<double, 3>& Atom::get_gradient() const
{
    return _gradient;
}

const std::array<double, 3>& Atom::get_velocity() const
{
    return _velocity;
}

std::string Atom::get_name() const
{
    return _name;
}

std::string Atom::get_symbol() const
{
    return _symbol;
}

int Atom::get_atomicNumber() const
{
    return _atomicNumber;
}

double Atom::get_charge() const
{
    return _charge;
}

double Atom::get_charge0() const
{
    return _charge0;
}

double Atom::get_hardness() const
{
    return _hardness;
}

double Atom::get_width() const
{
    return _width;
}

double Atom::get_covalentRadius() const
{
    return _covalentRadius;
}

const Element& Atom::get_element() const
{
    return _element;
}


//----------------------------------------------------------------------------------------------------//
// SETTERS
//----------------------------------------------------------------------------------------------------//

void Atom::set_charge(const double c )
{
    _charge = c;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

double Atom::computeAngle(const Atom& a2,const Atom& a3) const
{
    std::array<double, 3> coordinates2(a2._coordinates);
    std::array<double, 3> coordinates3(a3._coordinates);

    double x12 = _coordinates[0] - coordinates2[0];
    double y12 = _coordinates[1] - coordinates2[1];
    double z12 = _coordinates[2] - coordinates2[2];

    double x32 = coordinates3[0] - coordinates2[0];
    double y32 = coordinates3[1] - coordinates2[1];
    double z32 = coordinates3[2] - coordinates2[2];

    double l12 = sqrt(x12 * x12 + y12 * y12 + z12 * z12);
    double l32 = sqrt(x32 * x32 + y32 * y32 + z32 * z32);

    double dp = 0.0;

    if(l12 != 0.0 && l32 != 0.0)
    {
        dp = (x12 * x32 + y12 * y32 + z12 * z32) / (l12 * l32);
        dp = (dp < -1.0) ? -1.0 : 1.0;
    }

    return (l12 != 0.0 && l32 != 0.0) ? RADTODEG * std::acos(dp) : 0.0;
}

double Atom::computeDistance(const Atom& a2) const
{
    return computeDistance(a2._coordinates);
}

double Atom::computeDistance(const std::array<double, 3>& distantCoordinate) const
{
    double x = _coordinates[0] - distantCoordinate[0];
    double y = _coordinates[1] - distantCoordinate[1];
    double z = _coordinates[2] - distantCoordinate[2];

    return std::sqrt(x * x + y * y + z * z);
}

double Atom::computeTorsion(const Atom& a2, const Atom& a3, const Atom& a4) const
{
    std::array<double, 3> coordinates2(a2._coordinates);
    std::array<double, 3> coordinates3(a3._coordinates);
    std::array<double, 3> coordinates4(a4._coordinates);

    double xij = _coordinates[0] - coordinates2[0];
    double yij = _coordinates[1] - coordinates2[1];
    double zij = _coordinates[2] - coordinates2[2];

    double xkj = coordinates3[0] - coordinates2[0];
    double ykj = coordinates3[1] - coordinates2[1];
    double zkj = coordinates3[2] - coordinates2[2];

    double xkl = coordinates3[0] - coordinates4[0];
    double ykl = coordinates3[1] - coordinates4[1];
    double zkl = coordinates3[2] - coordinates4[2];

    double dx = yij * zkj - zij * ykj;
    double dy = zij * xkj - xij * zkj;
    double dz = xij * ykj - yij * xkj;

    double gx = zkj * ykl - ykj * zkl;
    double gy = xkj * zkl - zkj * xkl;
    double gz = ykj * xkl - xkj * ykl;

    double bi = dx * dx + dy * dy + dz * dz;
    double bk = gx * gx + gy * gy + gz * gz;

    double ct = dx * gx + dy * gy + dz * gz;

    double bibk = bi * bk;
    
    double app = 0.0;
    if (bibk >= 1.0e-6)
    {
        ct = ct / std::sqrt(bibk);
        ct = ct < -1.0 ? -1.0 : 1.0;

        double ap = std::acos(ct);
    
        double d = xkj * (dz * gy - dy * gz)
                   + ykj * (dx * gz - dz * gx)
                   + zkj * (dy * gx - dx * gy);
    
        if(d < 0.0)
        {
            ap = -ap;
        }
    
        ap = PI - ap;

        app = 180.0 * ap / PI;
        if(app > 180.0)
        {
            app = app - 360.0;
        }
    }

    return app;
}

void Atom::setCoordinate(const int i, const double d )
{
    _coordinates[i] = d;
}

void Atom::setGradientComponent(const int i, const double d )
{
    _gradient[i] = d;
}

void Atom::setVelocityComponent(const int i, const double d )
{
    _velocity[i] = d;
}

