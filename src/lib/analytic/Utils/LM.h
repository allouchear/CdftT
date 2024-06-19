#ifndef CDFTT_LM_H_INCLUDED
#define CDFTT_LM_H_INCLUDED

#include<iostream>
#include<vector>
#include<analytic/Utils/Utils.h>

using namespace std;

class LXYZ {

private:
	double _coef;
	vector<int> _l;

public:
	LXYZ();
	LXYZ(vector<int> L, double C);
	~LXYZ() {}
	double coef() {return _coef;}
	vector<int> l() {return _l;}
};

class Zlm {

private:
	int _l;
	int _m;
	int _numberOfCoefficients;
	vector<LXYZ> _lxyz;
	Binomial _bino;

public:
	Zlm();
	Zlm(int L, int M, Binomial& Bin);
	~Zlm() {}
	int numberOfCoefficients() {return _numberOfCoefficients;}
	vector<LXYZ> lxyz() {return _lxyz;}
	void setZlm0();
	void setCoefZlm();
};

void getlTable(int L, vector<int>& nCoefs, vector<vector<double>>& coefs, vector<vector<vector<int>>>& l, Binomial& Bin);

#endif