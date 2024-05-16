#ifndef CDFTT_CGTF_H_INCLUDED
#define CDFTT_CGTF_H_INCLUDED

#include<iostream>
#include<analytic/GTF.h>

using namespace std;

class CGTF
{
	private:
		int _num_center;
		int _numberOfFunctions;
		vector<GTF> _gtf;
		int _L;
		int _M;
		Binomial _bino;
	public:
		CGTF();
		CGTF(vector<GTF>);
		~CGTF(){}
		int numberOfFunctions()
		{
			return _numberOfFunctions;
		}
		vector<GTF> gtf()
		{
			return _gtf;
		}
		double ERICGTF(CGTF&, CGTF&, CGTF&);
		void normaliseCGTF();
		double overlapCGTF(CGTF&);
		double overlap3CGTF(CGTF&, CGTF&);
		double overlap4CGTF(CGTF&, CGTF&, CGTF&);
		double CGTFxyzCGTF(CGTF&, int, int, int);
		double kineticCGTF(CGTF&);
		double ionicPotentialCGTF(CGTF&, vector<double>, double);
		double CGTFstarCGTF(CGTF&);
		bool CGTFEqCGTF(CGTF&);
};

#endif
