#include<iostream>
#include<string>
#include<vector>
#include"Element.h"

Element::Element(string Name, string Symbol, int AtomicNumber, double CovalentRadii, double BondOrderRadii, 
	double VanDerWaalsRadii, double Radii, int MaximumBondValence, double Mass, 
	double Electronegativity) : _name(Name), _symbol(Symbol), _atomic_number(AtomicNumber), _covalent_radii(CovalentRadii), 
	_bond_order_radii(BondOrderRadii), 
	_van_der_waals_radii(VanDerWaalsRadii), _radii(Radii), _maximum_bond_valence(MaximumBondValence), _mass(Mass), _electronegativity(Electronegativity) {_isotope.resize(0);}

Element::Element() : _name("None"), _symbol("None"), _atomic_number(0), _covalent_radii(0), _bond_order_radii(0),
	_van_der_waals_radii(0), _radii(0), _maximum_bond_valence(0), _mass(0), _electronegativity(0) {}

Element::~Element(){}

Isotope::Isotope(string Symbol, int IntMass, double RealMass, double Abundance) : _symbol(Symbol), _int_mass(IntMass), _real_mass(RealMass), _abundance(Abundance) {}

Isotope::~Isotope(){}
