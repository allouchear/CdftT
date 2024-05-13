#ifndef _CDFTT_ATOM_H_INCLUDED
#define _CDFTT_ATOM_H_INCLUDED
using namespace std;
#include "PeriodicTable.h"
#include "Element.h"
#include <string>

	//! An Atom class.
	/*! This class will be used in the Structure class to create Structures and Molecules. */

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
	Element _e;
	
	string _mm_Type;
	string _pdb_Type;
	string _residue_name;
	int _residue_number;
	int _N;
	
	public:

	Atom();
		//!< Default constructor: 
		/*!<
		sets all attributes to 0 or "none"
		*/
	Atom(PeriodicTable& Table,const string& name);
		//!<Constructor: 
		/*!<
		creates an atom from the name of an Element and searching for it in PeriodicTable
		*/
	Atom(PeriodicTable& Table, const int& n);
		//!<Constructor:
		/*!<
		creates an atom from the atomic number of an Element and searching for it in PeriodicTable
		*/
	double* coordinates();
		//!< get() function
		/*!<
			Returns the coordinates of the atom
		*/
	double* gradient();
		//!< Get() function
		/*!<Returns the gradient of the atom
		*/
	double* velocity();
		//!< Get() function
		/*!< 
			Returns the velocity of the atom
		*/
	string name();
		//!< Get() function
		/*!<
			Returns the name of the atom
		*/
	string symbol();
		//!< Get() function 
		/*!<
			Returns the symbol of the atom
		*/
	int atomic_number();
		//!< Get() function
		/*!<
			Returns the atomic number of the atom
		*/
	double charge();
		//!< Get() function 
		/*!<
			Returns the partial charge of the atom
		*/
	double charge_0();
		//!< Get() function
		/*!<	
			Returns the oxidation of the atom
		*/
	double hardness();
		//!<Get() function
		/*!< 
			Returns the hardness of the atom
		*/
	double width();
		//!< Get() function
		/*!<
			Returns the width of the atom
		*/
	Element element();
		//!< Get() function
		/*!<	Returns the element associated with the atom */
		/*!<	The attributes of Element are properties of the atom	*/
		/*!<	Name			*/
		/*!<	Symbol			*/
		/*!<	Atomic Number		*/
		/*!<	Covalent Radius		*/
		/*!<	Bond order radii	*/
		/*!<	Van der Waals Radius	*/
		/*!<	Radius	*/
		/*!<	Maximum bond valence	*/
		/*!<	Mass	*/
		/*!<	electronegativity	*/
		
	~Atom();
		//!< Destructeur
		
	double get_distance(Atom& a2);
		//!< Renvoi la distance entre l'atom et un autre (a2)
		
	double get_angle(Atom& a2, Atom& a3);
		//!< Renvoi l'angle formée par l'atome et 2 autres(a2,a3)
		
	double get_torsion(Atom& a2, Atom& a3, Atom& a4);
		//!<Return
};

#endif //_CDFTT_ATOM_H_INCLUDED