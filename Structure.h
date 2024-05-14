#ifndef CDFTT_STRUCTURE_H_INCLUDED
#define CDFTT_STRUCTURE_H_INLCUDED

#include<iostream>
#include<vector>
#include"Atom.h"

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

			//! An default constructor.
			/*! In the case of a problem, this constructor create an element with all value on 0 and all string on "None". */

		Structure();

			//! .cub file reader
			/*! Reads .cub files to set values of _atoms' coordinates, velocity, ...*/
		void cube(ifstream& nameFile);
			
			//! A default desctructor.
			/*! We don't use it. */
		~Structure(){};

			//! A normal member taking no arguments and returning an integer value.
			/*! \return The number of atoms in our structure. */

		int number_of_atoms()
		{
			return _atoms.size();
		}

			//! A normal member taking one arguments and returning an atom value.
			/*! \return The atom i of our structure. */

		const Atom& atoms(const int& i)
		{
			return _atoms[i-1];
		}
		
};


#endif //CDFTT_STRUCTURE_H_INCLUDED
