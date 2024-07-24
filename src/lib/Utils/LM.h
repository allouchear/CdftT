#ifndef CDFTT_LM_H_INCLUDED
#define CDFTT_LM_H_INCLUDED

#include<iostream>
#include<vector>
#include <Utils/Utils.h>

using namespace std;

	//! A LXYZ class.
	/*! This class will be used in the Zlm class. */
class LXYZ {

private:
	double _coef;
	vector<int> _l;

public:

		//! A default constructor.
		/*! This constructor is used to set all of the parameters on 0 or "None" value. */
	LXYZ();

		//! A constructor taking two arguments.
		/*! This constructor is used to set all of the parameters. */
	LXYZ(vector<int> L, double C);

		//! A default desctructor.
		/*! We don't use it. */
	~LXYZ() {}

		//! A normal member taking no arguments and returning a double value.
		/*! \return The coefficient. */
	double coef() {return _coef;}

		//! A normal member taking no arguments and returning a vector<int> value.
		/*! \return The table of l (lx, ly, lz). */
	vector<int> l() {return _l;}
};

	//! A Zlm class.
	/*! This class will be used in a function. */
class Zlm {

private:
	int _l;
	int _m;
	int _numberOfCoefficients;
	vector<LXYZ> _lxyz;
	Binomial _bino;

public:

		//! A default constructor.
		/*! This constructor use setZlm0(). */
	Zlm();

		//! A constructor taking three arguments.
		/*! This constructor is used to set all of the parameters. */
	Zlm(int L, int M, Binomial& Bin);

		//! A default desctructor.
		/*! We don't use it. */
	~Zlm() {}

		//! A normal member taking no arguments and returning an int value.
		/*! \return The number of coefficient. */
	int numberOfCoefficients() {return _numberOfCoefficients;}

		//! A normal member taking no arguments and returning a vector<LXYZ> value.
		/*! \return The table of LXYZ. */
	vector<LXYZ> lxyz() {return _lxyz;}

		//! A normal member taking no arguments and returning a void value.
		/*! Set all of the parameters on 0 or "None" value. */
	void setZlm0();

		//! A normal member taking no arguments and returning a void value.
		/*! Calcul and update the coefficient of each LXYZ. */
	void setCoefZlm();
};

	//! A function taking two arguments and returning a void value.
	/*! Actualise the different tables of coefficients and the l to construct our basis in Orbitals. 
	 * This function use the Zlm classe. */
void getlTable(int L, vector<int>& nCoefs, vector<vector<double>>& coefs, vector<vector<vector<int>>>& l, Binomial& Bin);

#endif