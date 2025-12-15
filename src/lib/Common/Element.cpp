#include<string>

#include"Element.h"


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

Isotope::~Isotope()
{ }


Element::Element():
    _name("None"),
    _symbol("None"),
    _atomic_number(0),
    _covalent_radii(0.0),
    _bond_order_radii(0.0),
    _van_der_waals_radii(0.0),
    _radii(0.0),
    _maximum_bond_valence(0),
    _mass(0.0),
    _electronegativity(0.0) 
{
    _isotope.resize(1);
    _isotope[0]=Isotope();
}

Element::Element(const std::string& name, const std::string& symbol, const int atomicNumber, const double covalentRadii, const double bondOrderRadii, const double vanDerWaalsRadii, const double radii, const int maximumBondValence, const double mass, const double electronegativity):
    _name(name),
    _symbol(symbol),
    _atomic_number(atomicNumber),
    _covalent_radii(covalentRadii),
    _bond_order_radii(bondOrderRadii),
    _van_der_waals_radii(vanDerWaalsRadii),
    _radii(radii),
    _maximum_bond_valence(maximumBondValence),
    _mass(mass),
    _electronegativity(electronegativity)
{
    _isotope.resize(0);
}

Element::~Element()
{ }

