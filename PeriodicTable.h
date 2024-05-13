#ifndef PERIODICTABLE_H_INCLUDED
#define PERIODICTABLE_H_INCLUDED
#include<iostream>
#include<string>
#include<vector>

using namespace std;

class Isotope
{
	private:
		string _symbol;
		int _int_mass;
		double _real_mass;
		double _abundance;
	public:
		Isotope(){};
		Isotope(string, int, double, double);
		~Isotope();
		string symbol()
		{
			return _symbol;
		}

		int int_mass()
		{
			return _int_mass;
		}

		double real_mass()
		{
			return _real_mass;
		}

		double abundance()
		{
			return _abundance;
		}
};

class Element
{
	private:
		string _name;
		string _symbol;
		int _atomic_number;
		double _covalent_radii;
		double _bond_order_radii;
		double _van_der_waals_radii;
		double _radii;
		int _maximum_bond_valence;
		double _mass;
		double _electronegativity;
		vector<Isotope> _isotope;
	public:
		Element(){};
		Element(string, string, int, double, double, double, double, int, double, double);
		~Element();
		string name()
		{
			return _name;
		}

		string symbol()
		{
			return _symbol;
		}

		int atomic_number()
		{
			return _atomic_number;
		}

		double covalent_radii()
		{
			return _covalent_radii;
		}

		double bond_order_radii()
		{
			return _bond_order_radii;
		}

		double van_der_waals_radii()
		{
			return _van_der_waals_radii;
		}

		double radii()
		{
			return _radii;
		}

		int maximum_bond_valence()
		{
			return _maximum_bond_valence;
		}

		double mass()
		{
			return _mass;
		}

		double electronegativity()
		{
			return _electronegativity;
		}

		vector<Isotope> isotope()
		{
			return _isotope;
		}

		void push_isotope(Isotope ISO)
		{
			_isotope.push_back(ISO);
		}
};

class PeriodicTable
{
	private:
		vector<Element> _periodic_table;
	public:
		PeriodicTable();
		~PeriodicTable();
		Element element(int);
		Element element(string);
		void _add_element(Element);
		void _add_all_element();
		void _add_isotope(Isotope);
		void _add_all_isotope();
};

#endif
