#ifndef CDFTT_ELEMENT_H_INCLUDED
#define CDFTT_ELEMENT_H_INCLUDED
#include<iostream>
#include<string>
#include<vector>

using namespace std;

	//! An isotope class.
	/*! This class will be used in the Element class to complete the data on a element. */

class Isotope
{
	private:
		string _symbol;
		int _int_mass;
		double _real_mass;
		double _abundance;
	public:

			//! A default constructor.
			/*! In the case of a problem, this constructor create an isotope with all value on 0 and all string on "None". */

		Isotope();

			//! A real constructor.
			/*! This constructor is used to add all of the data of one isotope. */

		Isotope(const string&, const int&, const double&, const double&);

			//! A default desctructor.
			/*! We don't use it. */

		~Isotope();

			//! A normal member taking no arguments and returning a string value.
			/*! \return The isotope's symbol. */

		string symbol()
		{
			return _symbol;
		}

			//! A normal member taking no arguments and returning an integer value.
			/*! \return The approximate mass of an isotope. */

		int int_mass()
		{
			return _int_mass;
		}

			//! A normal member taking no arguments and returning an double value.
			/*! \return The real mass of an isotope. */

		double real_mass()
		{
			return _real_mass;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The abundace in percent of an isotope. */

		double abundance()
		{
			return _abundance;
		}
};

	//! An element class.
	/*! This class will be use as a dictionnary about all the properties of an element. */

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

			//! A default constructor.
			/*! In the case of a problem, this constructor create an element with all value on 0 and all string on "None". */

		Element();

			//! A real constructor.
			/*! This constructor is used to add all of the data of one element. */

		Element(const string&, const string&, const int&, const double&, const double&, const double&, const double&, const int&, const double&, const double&);

			//! A default destructor.
			/*! We don't use it. */

		~Element();

			//! A normal member taking no arguments and returning a string value.
			/*! \return The element's name. */

		string name()
		{
			return _name;
		}

			//! A normal member taking no arguments and returning a string value.
			/*! \return The element's symbol. */

		string symbol()
		{
			return _symbol;
		}

			//! A normal member taking no arguments and returning a integer value.
			/*! \return The atomic number of an element. */

		int atomic_number()
		{
			return _atomic_number;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The covalent radius of an element. */

		double covalent_radii()
		{
			return _covalent_radii;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The bond order radius of an element. */

		double bond_order_radii()
		{
			return _bond_order_radii;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The Van Der Waals radius of an element. */

		double van_der_waals_radii()
		{
			return _van_der_waals_radii;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The radius of an element. */

		double radii()
		{
			return _radii;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The maximum bond valence of an element. */

		int maximum_bond_valence()
		{
			return _maximum_bond_valence;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The mass of an element. */

		double mass()
		{
			return _mass;
		}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The electronegativity of an element. */

		double electronegativity()
		{
			return _electronegativity;
		}

			//! A normal member taking no arguments and returning an isotope value.
			/*! \return The isotope i of an element. 
				\sa Isotope class.
			*/

		Isotope isotope(const int& i)
		{
			return _isotope[i-1];
		}

			//! A normal member taking no arguments.
			/*! This add an isotope in the vector _isotope. */

		void push_isotope(const Isotope& ISO)
		{
			_isotope.push_back(ISO);
		}
};


#endif