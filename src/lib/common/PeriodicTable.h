#ifndef CDFTT_PERIODICTABLE_H_INCLUDED
#define CDFTT_PERIODICTABLE_H_INCLUDED
#include<iostream>
#include<vector>
#include"Element.h"

using namespace std;

	//! A periodic table class.
	/*! This class will be used in the program to have access of all the properties of an element. */

class PeriodicTable
{
	private:
		vector<Element> _periodic_table;
	public:

			//! A default constructor.
			/*! This create a periodic table. */

		PeriodicTable();

			//! A default desctructor.
			/*! We don't use it. */

		~PeriodicTable();

			//! A normal member taking one argument and returning an element value.
			/*! 
				\param i an integer argument (the atomic number).
			 	\return The element corresponding to the atomic number. */

		Element element(const int&);

			//! A normal member taking one argument and returning an element value.
			/*! 
				\param s a string argument (the symbol or the name).
			 	\return The element corresponding to the symbol or the name. */

		Element element(const string&);

			//! A normal member taking one argument.
			/*! Add an element.
				\param E an element argument. */

		void _add_element(const Element&);

			//! A normal member taking no arguments.
			/*! Add all elements and their isotopes  in the periodic table. */

		void _add_all_element();

			//! A normal member taking one argument.
			/*! Add an isotope in an element. 
				\param I an isotope argument. */

		void _add_isotope(Isotope);

			//! A normal member taking no arguments.
			/*! Add all isotopes for each elements.
				\sa _add_all_element. */

		void _add_all_isotope();
};

#endif
