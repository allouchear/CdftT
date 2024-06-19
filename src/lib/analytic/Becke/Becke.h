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
		bool _alpha_and_beta;
	public:
		Becke();
		Becke(WFX&, Binomial& Bin, const PeriodicTable&);
		~Becke() {}
		Structure str() {return _molecule;}
		int number_of_radial_points(int);
		GridPoints select_angular_grid(int);
		void multicenter_grids(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		vector<vector<double>> join_grids();
		double s(double, int);
		double multicenter_integration(function<double(const vector<GTF>&, double,double,double)>, const vector<GTF>&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double multicenter_integration(function<double(const vector<CGTF>&, double,double,double)>, const vector<CGTF>&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		double multicenter_integration(function<double(const LCAO&, const vector<vector<double>>&)>, const LCAO&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5, int alpha=0);
		vector<double> multicenter_sub_integration(function<double(int, const vector<vector<double>>&, const LCAO&, const vector<vector<vector<double>>>&, double, double, double, bool)> f, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		void partial_charge(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double density(int, const vector<vector<double>>&, const LCAO&, const vector<vector<vector<double>>>&, double, double, double, bool);
		double overlapGTF(const GTF&, const GTF&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodGTF(const vector<GTF>& p, double x, double y, double z);
		double overlapCGTF(const CGTF&, const CGTF&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
		static double prodCGTF(const vector<CGTF>& p, double x, double y, double z);
		double overlapLCAO(const LCAO&, int kmax=3, int lebedev_order=41, int radial_grid_factor=5, int alpha=0);
		static double prodLCAO(const LCAO& p, const vector<vector<double>>& d);
		vector<vector<double>> PartialChargeAndEnergy(int kmax=3, int lebedev_order=41, int radial_grid_factor=5);
};

#endif