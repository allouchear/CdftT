#ifndef CDFTT_PERIODICTABLE_H_INCLUDED
#define CDFTT_PERIODICTABLE_H_INCLUDED
#include<iostream>
#include<vector>
#include"Element.h"

using namespace std;

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
