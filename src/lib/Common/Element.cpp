#include <string>

#include "Element.h"


//----------------------------------------------------------------------------------------------------//
// CLASS ISOTOPE
//----------------------------------------------------------------------------------------------------//

Isotope::Isotope():
    _symbol("None"),
    _int_mass(0),
    _real_mass(0.0),
    _abundance(0.0)
{ }

Isotope::Isotope(const std::string& symbol, const int intMass, const double realMass, const double abundance):
    _symbol(symbol),
    _int_mass(intMass),
    _real_mass(realMass),
    _abundance(abundance)
{ }



//----------------------------------------------------------------------------------------------------//
// CLASS ELEMENT
//----------------------------------------------------------------------------------------------------//

//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Element::Element():
    _name("None"),
    _symbol("None"),
    _atomicNumber(0),
    _covalentRadius(0.0),
    _bondOrderRadius(0.0),
    _vanDerWaalsRadius(0.0),
    _radius(0.0),
    _maximumBondValence(0),
    _mass(0.0),
    _electronegativity(0.0) 
{
    _isotope.resize(1);
    _isotope[0]=Isotope();
}

Element::Element(const std::string& name, const std::string& symbol, const int atomicNumber, const double covalentRadii, const double bondOrderRadii, const double vanDerWaalsRadii, const double radii, const int maximumBondValence, const double mass, const double electronegativity):
    _name(name),
    _symbol(symbol),
    _atomicNumber(atomicNumber),
    _covalentRadius(covalentRadii),
    _bondOrderRadius(bondOrderRadii),
    _vanDerWaalsRadius(vanDerWaalsRadii),
    _radius(radii),
    _maximumBondValence(maximumBondValence),
    _mass(mass),
    _electronegativity(electronegativity)
{
    _isotope.resize(0);
}


//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

const std::string& Element::get_name() const
{
    return _name;
}

const std::string& Element::get_symbol() const
{
    return _symbol;
}

int Element::get_atomicNumber() const
{
    return _atomicNumber;
}

double Element::get_covalentRadius() const
{
    return _covalentRadius;
}

double Element::get_bondOrderRadius() const
{
    return _bondOrderRadius;
}

double Element::get_vanDerWaalsRadius() const
{
    return _vanDerWaalsRadius;
}

double Element::get_radius() const
{
    return _radius;
}

int Element::get_maximumBondValence() const
{
    return _maximumBondValence;
}

double Element::get_mass() const
{
    return _mass;
}

double Element::get_electronegativity() const
{
    return _electronegativity;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

Isotope Element::getIsotope(const int i) const
{
    return _isotope[i - 1];
}

void Element::pushIsotope(const Isotope& isotope)
{
    _isotope.push_back(isotope);
}