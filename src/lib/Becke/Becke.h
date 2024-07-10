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
		Becke();
		Becke(const Structure& s);
		Becke(const Grid& g);
		Becke(WFX&, Binomial&, const PeriodicTable&);
		Becke(FCHK&, Binomial&, const PeriodicTable&);
		Becke(MOLDENGAB&, Binomial&, const PeriodicTable&);
		Becke(LOG&, Binomial&, const PeriodicTable&);
		~Becke() {}
		Structure str() {return _molecule;}
		int number_of_radial_points(int);
		GridPoints select_angular_grid(int);
		void multicenter_grids(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		vector<vector<double>> join_grids();
		double s(double, int);
		double multicenter_integration(function<double(const vector<GTF>&, double,double,double)>, const vector<GTF>&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double multicenter_integration(function<double(Orbitals&, int, int, double,double,double)>, int, int, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double multicenter_integration(function<double(Orbitals&, int, int, double, double, double, int)> f, int i, int j, int kmax=3, int lebedev_order=41, int radial_grid_factor=5, int alpha=0);
		vector<double> multicenter_sub_integration(function<double(Orbitals&, double, double, double)> f, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		void partial_charge(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double density(Orbitals&, double, double, double);
		double OverlapGTF(const GTF&, const GTF&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodGTF(const vector<GTF>& p, double x, double y, double z);
		double overlapCGTF(const CGTF&, const CGTF&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodCGTF(const vector<CGTF>& p, double x, double y, double z);
		double multicenter_integration(const Grid& g, int kmax=3 , int lebedev_order=41, int radial_grid_factor=5);
		vector<double> multicenter_sub_integration(const Grid& g, int kmax=3 , int lebedev_order=41, int radial_grid_factor=5);
		void partial_charge(const Grid& g, int kmax=3 , int lebedev_order=41, int radial_grid_factor=5);
		double OverlapCGTF(int i, int j, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double CGTFstarCGTF(Orbitals& Orb, int i, int j, double x, double y, double z);
		double Overlap(int i, int j, int kmax=3, int lebedev_order=41, int radial_grid_factor=5, int alpha=0);
		static double phi(Orbitals& Orb, int i, double x, double y, double z, int alpha=0);
		static double phistarphi(Orbitals& Orb, int i, int j, double x, double y, double z, int alpha=0);

		double eHOMO() {_orbitals.HOMO(); return _orbitals.eHOMO();}
		double eLUMO() {_orbitals.LUMO(); return _orbitals.eLUMO();}
		vector<double> PartialChargeAndEnergy(const Grid& g, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		vector<vector<double>> PartialChargesAndEnergy(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double get_Energy();
		vector<double> get_Partial_Charge();
		void printCharges();
};
#endif
