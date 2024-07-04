#ifndef _CDFTT_ATOM_H_INCLUDED
#define _CDFTT_ATOM_H_INCLUDED
using namespace std;
#include <Common/PeriodicTable.h>
#include <Common/Element.h>
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
	double _covalent_radii;
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
	Atom(const PeriodicTable& Table,const string& name);
		//!<Constructor: 
		/*!<
		creates an atom from the name of an Element and searching for it in PeriodicTable
		*/
	Atom(const PeriodicTable& Table, const int& n);
		//!<Constructor:
		/*!<
		creates an atom from the atomic number of an Element and searching for it in PeriodicTable
		*/
	const double* coordinates() const;
		//!< get() function
		/*!<
			Returns the coordinates of the atom
		*/
		
	double* gradient();
		//!< Get() function
		/*!<
			Returns the gradient of the atom
		*/
	double* velocity();
		//!< Get() function
		/*!< 
			Returns the velocity of the atom
		*/
	string name() const;
		//!< Get() function
		/*!<
			Returns the name of the atom
		*/
	string symbol() const;
		//!< Get() function 
		/*!<
			Returns the symbol of the atom
		*/
	int atomic_number() const;
		//!< Get() function
		/*!<
			Returns the atomic number of the atom
		*/
	double charge() const;
		//!< Get() function 
		/*!<
			Returns the partial charge of the atom
		*/
	double charge_0() const;
		//!< Get() function
		/*!<	
			Returns the oxidation of the atom
		*/
	double hardness() const;
		//!<Get() function
		/*!< 
			Returns the hardness of the atom
		*/
	double width() const;
		//!< Get() function
		/*!<
			Returns the width of the atom
		*/
	double covalent_radii() const;
		//!< Get() function
		/*!<
			Returns the covalent radii of the atom
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
		
	void Set_coordinates(const int& i, const double& d);
		//!< Set() function
		/*!<
			Sets the ith coordinate to d
		*/
		
	void Set_gradient(const int& i, const double& d);
		//!< Set() function
		/*!<
			Sets the ith gradient component to d
		*/
	void Set_velocity(const int& i, const double& d);
		//!< Set() function
		/*!<
			Sets the ith velocity component to d
		*/
	void Set_charge(const double& c);
		//!< Set() function
		/*!<
			Sets the charge to c
		*/
	~Atom();
		//!< Destructeur
		
	double get_distance(const Atom& a2) const;
		//!< Renvoi la distance entre l'atom et un autre (a2)
		
	double get_angle(const Atom& a2,const Atom& a3) const;
		//!< Renvoi l'angle formée par l'atome et 2 autres(a2,a3)
		
	double get_torsion(const Atom& a2,const Atom& a3,const Atom& a4) const;
		//!<Return
};

#endif //_CDFTT_ATOM_H_INCLUDED
