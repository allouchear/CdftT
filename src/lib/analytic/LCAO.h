#ifndef CDFTT_LCAO_H_INCLUDED
#define CDFTT_LCAO_H_INCLUDED

#include<iostream>
#include<analytic/CGTF.h>

using namespace std;

class LCAO
{
	private:
		int _numberOfFunctions;
		vector<double> _coefficient;
		vector<CGTF> _cgtf;
		Binomial _bino;
	public:
		LCAO();
		LCAO(vector<CGTF>);
		~LCAO(){}
		vector<CGTF> cgtf()
		{
			return _cgtf;
		}
		vector<double> coefficient()
		{
			return _coefficient;
		}
		int numberOfFunctions()
		{
			return _numberOfFunctions;
		}
		double ERICLCAO(LCAO&, LCAO&, LCAO&);
		void normaliseLCAO();
		double overlapLCAO(LCAO&);
		double overlap3LCAO(LCAO&, LCAO&);
		double overlap4LCAO(LCAO&, LCAO&, LCAO&);
		double kineticLCAO(LCAO&);
		double ionicPotentialLCAO(LCAO&, vector<double>, double);
		double LCAOstarLCAO(LCAO&);
		double LCAOxyzLCAO(LCAO&, int, int, int);
		//bool LCAOEqLCAO(LCAO&);
};

bool operator==(LCAO, LCAO);

#endif