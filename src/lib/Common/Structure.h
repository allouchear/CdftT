#ifndef CDFTT_STRUCTURE_H_INCLUDED
#define CDFTT_STRUCTURE_H_INCLUDED

#include<iostream>
#include<vector>
#include <Common/Atom.h>
#include <Utils/WFX.h>
#include <Utils/FCHK.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/LOG.h>

using namespace std;

	//! A structure class.
	/*! This class will be used to have molecules or periodic structures. */
class Structure
{
	private:
		vector<Atom> _atoms;
	public:

			//! A real constructor.
			/*! This constructor is used to add all of the data of all atoms used in our structure. */
		Structure(const vector<Atom>&);

			//! Default constructor.
			/*! In the case of a problem, this constructor create an element with all value on 0 and all string on "None". */
		Structure();

			//! .cube file reader
			/*! Reads .cube files to initialize the atoms in the structure and the structure itself.*/
			/*! note: uses Atom::Atom(const PeriodicTable& Table, const int& i), Table needs to be declared beforehand */
			/*! note 2: c++ sequential read requires that Domain::Domain(ifstream& nameFile) be called first */
		void read_From_Cube(ifstream& nameFile, int Natoms, const PeriodicTable& Table);
		
			//*Constructor
			/*!
				Calls read_From_Cube() to initialize the atoms of the structure
			*/
		Structure(ifstream& nameFile, const int Natoms, const PeriodicTable& Table);
			
			//! A default desctructor.
			/*! We don't use it. */
		~Structure(){}

			//! A normal member taking no arguments and returning an integer value.
			/*! \return The number of atoms in our structure. */
		int number_of_atoms() const
		{
			return _atoms.size();
		}

			//! A normal member taking one arguments and returning an atom value.
			/*! \return The atom i of our structure. */
		Atom atom(const int& i) const
		{
			return _atoms[i];
		}
			//! Get() function
			/*! returns _atoms as a vector of atoms */
		vector<Atom> atoms() const
		{
			return _atoms;
		}
		
			//! Operator +
			/*! Overload of + returns the structure with the biggest number of atoms */
		Structure operator+(const Structure& S) const;
		
			//! Operator -
			/*! Overload of - returns the structure with the biggest number of atoms */
		Structure add(const Structure& S);

			//! .wfx file reader
			/*! Reads .wfx files to initialize the atoms in the structure and the structure itself.*/
			/*! note: uses Atom::Atom(const PeriodicTable& Table, const int& i), Table needs to be declared beforehand */
		void read_from_wfx(WFX&, const PeriodicTable&);

			//*Constructor
			/*!
				Calls read_from_wfx() to initialize the atoms of the structure
			*/
		Structure(WFX&, const PeriodicTable&);

			//! .fchk file reader
			/*! Reads .fchk files to initialize the atoms in the structure and the structure itself.*/
			/*! note: uses Atom::Atom(const PeriodicTable& Table, const int& i), Table needs to be declared beforehand */
		void read_from_fchk(FCHK&, const PeriodicTable&);

			//*Constructor
			/*!
				Calls read_from_fchk() to initialize the atoms of the structure
			*/
		Structure(FCHK&, const PeriodicTable&);

			//! .molden and .gab file reader
			/*! Reads .molden and .gab files to initialize the atoms in the structure and the structure itself.*/
			/*! note: uses Atom::Atom(const PeriodicTable& Table, const int& i), Table needs to be declared beforehand */
		void read_from_moldengab(MOLDENGAB&, const PeriodicTable&);

			//*Constructor
			/*!
				Calls read_from_moldengab() to initialize the atoms of the structure
			*/
		Structure(MOLDENGAB&, const PeriodicTable&);

			//! .log file reader
			/*! Reads .log files to initialize the atoms in the structure and the structure itself.*/
			/*! note: uses Atom::Atom(const PeriodicTable& Table, const int& i), Table needs to be declared beforehand */
		void read_from_log(LOG&, const PeriodicTable&);

			//*Constructor
			/*!
				Calls read_from_log() to initialize the atoms of the structure
			*/
		Structure(LOG&, const PeriodicTable&);
};

#endif //CDFTT_STRUCTURE_H_INCLUDED
