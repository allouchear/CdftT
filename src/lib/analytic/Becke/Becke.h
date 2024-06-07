#ifndef CDFTT_BECKE_INCLUDED
#define CDFTT_BECKE_INCLUDED

#include<iostream>
#include<functional>
#include<analytic/Utils/WFX.h>
#include<analytic/Orbitals/Orbitals.h>
#include<analytic/Becke/GridPoints.h>
#include<common/Structure.h>

using namespace std;

class Becke
{
	private:
		Structure _molecule;
		GridPoints _grid;
		vector<vector<vector<double>>> _grid_points;
		vector<vector<double>> _grid_weights;
		vector<vector<double>> _grid_volumes;
	public:
		Becke();
		Becke(WFX&, const PeriodicTable&);
		~Becke() {}
		int number_of_radial_points(int);
		GridPoints select_angular_grid(int);
		void multicenter_grids(int, int, int);
		vector<vector<double>> join_grids();
		double s(double, int);
		double multicenter_integration(function<double(vector<GTF>, double,double,double)>, vector<GTF>, double, double, double, int, int);
		double overlap(GTF, GTF);
		static double prod(vector<GTF> p, double x, double y, double z);

};

#endif