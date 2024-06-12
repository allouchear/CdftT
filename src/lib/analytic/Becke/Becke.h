#ifndef CDFTT_BECKE_INCLUDED
#define CDFTT_BECKE_INCLUDED

#include<iostream>
#include<functional>
#include<analytic/Utils/WFX.h>
#include<analytic/Orbitals/Orbitals.h>
#include<analytic/Becke/GridPoints.h>
#include<common/Descriptors.h>

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
		Becke(WFX&, Binomial& Bin, const PeriodicTable&);
		~Becke() {}
		int number_of_radial_points(int);
		GridPoints select_angular_grid(int);
		void multicenter_grids(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		vector<vector<double>> join_grids();
		double s(double, int);
		double multicenter_integration(function<double(vector<GTF>, double,double,double)>, vector<GTF>, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double multicenter_integration(function<double(vector<CGTF>, double,double,double)>, vector<CGTF>, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double multicenter_integration(function<double(vector<LCAO>, double,double,double)>, vector<LCAO>, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		vector<double> multicenter_sub_integration(function<double(int, vector<double>, vector<LCAO>, double, double, double)> f, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		void partial_charge(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double density(int, vector<double>, vector<LCAO>, double, double, double);
		double overlapGTF(GTF, GTF, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodGTF(vector<GTF> p, double x, double y, double z);
		double overlapCGTF(CGTF, CGTF, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodCGTF(vector<CGTF> p, double x, double y, double z);
		double overlapLCAO(LCAO, LCAO, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodLCAO(vector<LCAO> p, double x, double y, double z);
		vector<double> PartialChargeAndEnergy(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
};

#endif