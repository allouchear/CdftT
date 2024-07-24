#ifndef CDFTT_BECKE_INCLUDED
#define CDFTT_BECKE_INCLUDED

#include<iostream>
#include<functional>
#include <Utils/WFX.h>
#include <Utils/FCHK.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/LOG.h>
#include <Orbitals/Orbitals.h>
#include <Becke/GridPoints.h>
#include <Common/Descriptors.h>

using namespace std;

	//! A LOG class.
	/*! This class will be use to calculate partial charges. */
class Becke
{
	private:
		Structure _molecule;
		Orbitals _orbitals;
		GridPoints _grid;
		vector<vector<vector<double>>> _grid_points;
		vector<vector<double>> _grid_weights;
		vector<vector<double>> _grid_volumes;
		vector<double> _partial_charge;
		double _energy;
		bool _multigrid;
	public:
		double get_Energy();
		vector<double> get_Partial_Charge();
		void printCharges();

			//! A default constructor.
			/*! This constructor is used to set all of the parameters on 0 or "None" value. */
		Becke();

			//! A constructor taking one argument.
			/*! This constructor is used to setup the structure. */
		Becke(const Structure& s);

			//! A constructor taking one argument.
			/*! This constructor is used to setup the grid. */
		Becke(const Grid& g);

			//! A constructor taking one argument.
			/*! This constructor is used to set all of the parameters with a .wfx file. */
		Becke(WFX&, Binomial&, const PeriodicTable&);

			//! A constructor taking one argument.
			/*! This constructor is used to set all of the parameters with a .fchk file. */
		Becke(FCHK&, Binomial&, const PeriodicTable&);

			//! A constructor taking one argument.
			/*! This constructor is used to set all of the parameters with a .molden or a .gab file. */
		Becke(MOLDENGAB&, Binomial&, const PeriodicTable&);

			//! A constructor taking one argument.
			/*! This constructor is used to set all of the parameters with a .log file. */
		Becke(LOG&, Binomial&, const PeriodicTable&);

			//! A default desctructor.
			/*! We don't use it. */
		~Becke() {}

			//! A normal member taking no arguments and returning a Structure value.
			/*! \return The structure. */
		Structure str() {return _molecule;}

			//! A normal member taking one argument and returning an int value.
			/*! \return The number of radial points for an atomic number. */
		int number_of_radial_points(int);

			//! A normal member taking one argument and returning a GridPoints value.
			/*! \return The select grid for a Lmax. */
		GridPoints select_angular_grid(int);

			//! A normal member taking three arguments and returning a void value.
			/*! Construct and actualise _grid_points, _grid_weights, and _grid_volumes. 
			 * Construct a Becke grid for each atom.*/
		void multicenter_grids(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A normal member taking no arguments and returning a vector<vector<double>> value.
			/*! \return One grid (x,y,z,w).
			 * Merge all the grid in only one and recaculate the weights for each points. */
		vector<vector<double>> join_grids();

			//! A normal member taking two arguments and returning a double value.
			/*! \return The value of cutoff profiles. */
		double s(double, int);

			//! A normal member taking five arguments and returning a double value.
			/*! \return The value of an integral for a function like f(vector<GTF>, x, y, z). */
		double multicenter_integration(function<double(const vector<GTF>&, double,double,double)>, const vector<GTF>&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A normal member taking six arguments and returning a double value.
			/*! \return The value of an integral for a function like f(Orbitals, i, j, x, y, z). */
		double multicenter_integration(function<double(Orbitals&, int, int, double,double,double)>, int, int, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A normal member taking seven arguments and returning a double value.
			/*! \return The value of an integral for a function like f(Orbitals, i, j, x, y, z, spin). */
		double multicenter_integration(function<double(Orbitals&, int, int, double, double, double, int)> f, int i, int j, int kmax=3, int lebedev_order=41, int radial_grid_factor=5, int alpha=0);

			//! A normal member taking four arguments and returning a vector<double> value.
			/*! \return The table of values of an integral for a function like f(Orbitals, x, y, z) on each grid (so on each atom). */
		vector<double> multicenter_sub_integration(function<double(Orbitals&, double, double, double)> f, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A normal member taking three arguments and returning a void value.
			/*! Calculate and actualise _partial_charge. */
		void partial_charge(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A static member taking four arguments and returning a double value.
			/*! \return The value of the electronic density at the point (x,y,z). */
		static double density(Orbitals&, double, double, double);

			//! A normal member taking five arguments and returning a double value.
			/*! \return The value of an integral to have the overlap between two GTF. */
		double OverlapGTF(const GTF&, const GTF&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A static member taking four arguments and returning a double value.
			/*! \return The value of the product between many GTF at the point (x,y,z). */
		static double prodGTF(const vector<GTF>& p, double x, double y, double z);

			
		double multicenter_integration(const Grid& g, int kmax=3 , int lebedev_order=41, int radial_grid_factor=5);
		vector<double> multicenter_sub_integration(const Grid& g, int kmax=3 , int lebedev_order=41, int radial_grid_factor=5);
		void partial_charge(const Grid& g, int kmax=3 , int lebedev_order=41, int radial_grid_factor=5);


			//! A normal member taking five arguments and returning a double value.
			/*! \return The value of an integral to have the overlap between two CGTF i and j. */
		double OverlapCGTF(int i, int j, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);

			//! A static member taking six arguments and returning a double value.
			/*! \return The value of the product between two CGTF at the point (x,y,z). */
		static double CGTFstarCGTF(Orbitals& Orb, int i, int j, double x, double y, double z);

			//! A normal member taking six arguments and returning a double value.
			/*! \return The value of an integral to have the overlap between two Orbitals i and j. */
		double Overlap(int i, int j, int kmax=3, int lebedev_order=41, int radial_grid_factor=5, int alpha=0);

			//! A static member taking six arguments and returning a double value.
			/*! \return The value of the Orbital i at the point (x,y,z). */
		static double phi(Orbitals& Orb, int i, double x, double y, double z, int alpha=0);

			//! A static member taking six arguments and returning a double value.
			/*! \return The value of the product between two Orbitals i and j at the point (x,y,z). */
		static double phistarphi(Orbitals& Orb, int i, int j, double x, double y, double z, int alpha=0);

			//! A normal member taking no arguments and returning a double value.
			/*! \return The HOMO energy. */
		double eHOMO() {_orbitals.HOMO(); return _orbitals.eHOMO();}

			//! A normal member taking no arguments and returning a double value.
			/*! \return The LUMO energy. */
		double eLUMO() {_orbitals.LUMO(); return _orbitals.eLUMO();}


		vector<double> PartialChargeAndEnergy(const Grid& g, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);


			//! A normal member taking five arguments and returning a vector<double> value.
			/*! \return The energy (index 0) and all the partial charges (index 1). */
		vector<vector<double>> PartialChargesAndEnergy(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
};

#endif
