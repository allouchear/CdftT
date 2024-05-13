#ifndef _CDFTT_ATOM_H_INCLUDED
#define _CDFTT_ATOM_H_INCLUDED
using namespace std;
#include "PeriodicTable.h"
#include <string>

class Atom
{
	double _coordinates[3];
	double _gradient[3];
	double _velocity[3];
	string _name;
	string _symbol;
	int _atomic_number;
	double _charge; // partial charge
	double _charge_0; // oxidation
	double _hardness; // eta
	double _width; // eta
	Element e;
	
	string _mm_Type;
	string _pdb_Type;
	string _residue_name;
	int _residue_number;
	int _N;
	public:
	Atom();
	Atom(const PeriodicTable& table, const string& name);
	~Atom();
	double get_distance(Atom& a2);
	double get_angle(Atom& a2, Atom& a3);
	double get_torsion(Atom& a2, Atom& a3, Atom& a4);
};

#endif //_CDFTT_ATOM_H_INCLUDED
