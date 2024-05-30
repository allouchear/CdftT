#include<iostream>
#include<analytic/CGTF.h>

using namespace std;

CGTF::CGTF()
{
	_gtf.resize(0);
	_numberOfFunctions=0;
	_L=0;
	_M=0;
	_bino = Binomial();
}

CGTF::CGTF(vector<GTF> A) : _gtf(A)
{
	_numberOfFunctions=_gtf.size();
	_L=0;
	_M=0;
	_bino=A[0].bino();
}

double CGTF::ERICGTF(CGTF& q, CGTF& r, CGTF& s)
{
	int np,nq;
	int nr,ns;
	double sum = 0.0;

	for(np=0;np<_numberOfFunctions;np++)
		for(nq=0;nq<q.numberOfFunctions();nq++)
			for(nr=0;nr<r.numberOfFunctions();nr++)
				for(ns=0;ns<s.numberOfFunctions();ns++)
					sum += _gtf[np].ERIGTF(q.gtf()[nq],r.gtf()[nr],s.gtf()[ns]); 

	return sum;
}

void CGTF::normaliseCGTF()
{
	int n,np;
	double sum=0.0;
	for(n=0 ; n<_numberOfFunctions ; n++)
		sum += _gtf[n].coefficient()*_gtf[n].coefficient()* _gtf[n].GTFstarGTF(_gtf[n]);

	for(n=0;n<_numberOfFunctions-1 ; n++)
		for(np=n+1; np<_numberOfFunctions; np++)
			sum += 2*_gtf[n].coefficient()*_gtf[np].coefficient()*_gtf[n].GTFstarGTF(_gtf[np]);

	if(sum>1.e-20)
	{
		sum = sqrt(sum);
		for(n=0 ; n<_numberOfFunctions ; n++)
			_gtf[n]/=sum;
	}
	else
	{
		cout<<"A Contacted Gaussian Type function is nul"<<endl;
		exit(1);
	}
}

double CGTF::overlapCGTF(CGTF& right)
{
	double sum=0.0;

	//cout<<"Number of Function GTF in CGTF = "<<_numberOfFunctions<<endl;

	for(int n=0;n<_numberOfFunctions;n++)
		for(int np=0;np<_numberOfFunctions;np++)
			sum += _gtf[n].overlapGTF(right._gtf[np]);

	//cout<<"Test Overlap in CGTF "<<endl;

	return sum;
}

double CGTF::overlap3CGTF(CGTF& midle, CGTF& right)
{
	double sum=0.0;
	int n;
	int np;
	int ns;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<midle.numberOfFunctions();np++)
			for(ns=0;ns<right.numberOfFunctions();ns++)
				sum += _gtf[n].overlap3GTF(midle.gtf()[np],right.gtf()[ns]);

	return sum;
}

double CGTF::overlap4CGTF(CGTF& B, CGTF& C, CGTF& D)
{
	double sum=0.0;
	int np;
	int nq;
	int nr;
	int ns;

	for(np=0;np<_numberOfFunctions;np++)
		for(nq=0;nq<B.numberOfFunctions();nq++)
			for(nr=0;nr<C.numberOfFunctions();nr++)
				for(ns=0;ns<D.numberOfFunctions();ns++)
					sum += _gtf[np].overlap4GTF(B.gtf()[nq],C.gtf()[nr],D.gtf()[ns]);

	return sum;
}

double CGTF::CGTFxyzCGTF(CGTF& right, int ix, int iy, int iz)
{
	double sum=0.0;
	int n;
	int ns;
	vector<double> C(3,0);
	vector<int> l {ix, iy, iz};
	GTF m(0.0, 1.0, C, l, _bino);

	for(n=0;n<_numberOfFunctions;n++)
		for(ns=0;ns<right.numberOfFunctions();ns++)
			sum += _gtf[n].overlap3GTF(m,right.gtf()[ns]);

	return sum;
}

double CGTF::kineticCGTF(CGTF& right)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<right.numberOfFunctions();np++)
			sum += _gtf[n].kineticGTF(right.gtf()[np]);


	return sum;
}

double CGTF::ionicPotentialCGTF(CGTF& right, vector<double> C, double Z)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<right.numberOfFunctions();np++)
			sum += _gtf[n].ionicPotentialGTF(right.gtf()[np], C, Z); 

	return sum;
}

double CGTF::CGTFstarCGTF(CGTF& right)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<right._numberOfFunctions;np++)
			sum += _gtf[n].GTFstarGTF(right.gtf()[np]);

	return sum;
}
/*
bool CGTF::CGTFEqCGTF(CGTF& t2)
{
	int i;
	int c;
	if(_numberOfFunctions != t2.numberOfFunctions()) return false;
	for(i=0;i<3;i++)
		if(_gtf[0].l()[i] != t2.gtf()[0].l()[i]) return false;
	for(i=0;i<_numberOfFunctions;i++)
	{
		if(fabs(_gtf[i].exposant()-t2.gtf()[i].exposant())>1e-10) return false;
		if(fabs(_gtf[i].coefficient()-t2.gtf()[i].coefficient())>1e-10) return false;
		for(c=0;c<3;c++)
			if(fabs(_gtf[i].coord()[c]-t2.gtf()[i].coord()[c])>1e-10) return false;
	}
	return true;
}
*/

void CGTF::push_back(GTF& gtf)
{
	_gtf.push_back(gtf);
	_numberOfFunctions++;
}

bool operator==(CGTF a, CGTF b)
{
	size_t i,j;
	size_t c=0;
	if(a.gtf().size()!=b.gtf().size())
		return false;
	for(i=0; i<a.gtf().size(); i++)
		for(j=0; j<a.gtf().size(); j++)
			if(a.gtf()[i]==b.gtf()[j])
				c++;
	if(c==a.gtf().size())
		return true;
	else
		return false;
}

ostream& operator<<(ostream& flux, CGTF& cgtf)
{
	for(int i=0; i<cgtf.numberOfFunctions(); i++)
	{
		if(i>0)
			flux<<setw(20)<<" ";
		flux<<cgtf.gtf()[i]<<endl;
	}
	return flux;
}