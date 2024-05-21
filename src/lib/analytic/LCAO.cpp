#include<iostream>
#include<analytic/LCAO.h>

using namespace std;

LCAO::LCAO()
{
	_cgtf.resize(0);
	_numberOfFunctions=0;
	_bino=Binomial();
}

LCAO::LCAO(vector<CGTF> A) : _cgtf(A)
{
	_numberOfFunctions=_cgtf.size();
	_bino=A[0].bino();
}

double LCAO::ERICLCAO(LCAO& q, LCAO& r, LCAO& s)
{
	int np,nq;
	int nr,ns;
	double sum = 0.0;

	for(np=0;np<_numberOfFunctions;np++)
	for(nq=0;nq<q.numberOfFunctions();nq++)
	for(nr=0;nr<r.numberOfFunctions();nr++)
	for(ns=0;ns<s.numberOfFunctions();ns++)
		sum += _cgtf[np].ERICGTF(q.cgtf()[nq],r.cgtf()[nr],s.cgtf()[ns]); 

	return sum;
}

void LCAO::normaliseLCAO()
{
	int n,ns,ng,np;
	double sum=0.0;
	for(n=0 ; n<_numberOfFunctions ; n++)
		for(ns=0; ns<_cgtf[n].numberOfFunctions(); ns++)
			sum += _cgtf[n].gtf()[ns].coefficient()*_cgtf[n].gtf()[ns].coefficient()* _cgtf[n].CGTFstarCGTF(_cgtf[n]);

	for(n=0;n<_numberOfFunctions-1 ; n++)
		for(np=n+1; np<_numberOfFunctions; np++)
			for(ns=0 ; ns<_cgtf[n].numberOfFunctions() ; ns++)
				for(ng=0; ng<_cgtf[np].numberOfFunctions() ; ng++)
					sum += 2*_cgtf[n].gtf()[ns].coefficient()*_cgtf[np].gtf()[ng].coefficient()*_cgtf[n].CGTFstarCGTF(_cgtf[np]);

	if(sum>1.e-20)
	{
		sum = sqrt(sum);
		for(n=0 ; n<_numberOfFunctions ; n++)
			for(ng=0 ; ng<_cgtf[n].numberOfFunctions() ; ng++)
			_cgtf[n].gtf()[ng]/=sum;
	}
	else
	{
		cout<<"A Contacted Gaussian Type function is nul"<<endl;
		exit(1);
	}
}

double LCAO::overlapLCAO(LCAO& right)
{
	double sum=0.0;
	int n;
	int np;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<_numberOfFunctions;np++)
			sum += _cgtf[n].overlapCGTF(right.cgtf()[np]);

	return sum;
}

double LCAO::overlap3LCAO(LCAO& midle, LCAO& right)
{
	double sum=0.0;
	int n;
	int np;
	int ns;

	for(n=0;n<_numberOfFunctions;n++)
	for(np=0;np<midle.numberOfFunctions();np++)
	for(ns=0;ns<right.numberOfFunctions();ns++)
			sum += _cgtf[n].overlap3CGTF(midle.cgtf()[np],right.cgtf()[ns]);

	return sum;
}

double LCAO::overlap4LCAO(LCAO& B, LCAO& C, LCAO& D)
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
			sum += _cgtf[np].overlap4CGTF(B.cgtf()[nq],C.cgtf()[nr],D.cgtf()[ns]);

	return sum;
}

double LCAO::kineticLCAO(LCAO& right)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<right.numberOfFunctions();np++)
			sum += _cgtf[n].kineticCGTF(right.cgtf()[np]);


	return sum;
}

double LCAO::ionicPotentialLCAO(LCAO& right, vector<double> C, double Z)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<right.numberOfFunctions();np++)
			sum += _cgtf[n].ionicPotentialCGTF(right.cgtf()[np], C, Z); 

	return sum;
}

double LCAO::LCAOstarLCAO(LCAO& right)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<right.numberOfFunctions();np++)
			sum += _cgtf[n].CGTFstarCGTF(right.cgtf()[np]);

	return sum;
}


double LCAO::LCAOxyzLCAO(LCAO& right, int ix, int iy, int iz)
{
	double sum=0.0;
	int n;
	int ns;
	vector<double> C(3,0);
	vector<int> l {ix, iy, iz};
	GTF m1(0.0, 1.0, C, l, _bino);
	vector<GTF> mbis (1,m1);
	CGTF m2(mbis);

	for(n=0;n<_numberOfFunctions;n++)
		for(ns=0;ns<right.numberOfFunctions();ns++)
				sum += _cgtf[n].gtf()[ns].overlap3GTF(m2.gtf()[0],right.cgtf()[n].gtf()[ns]);

	return sum;
}

/*
bool LCAO::LCAOEqLCAO(LCAO& t2)
{
	int i;
	int c;
	if(_numberOfFunctions != t2.numberOfFunctions()) return false;
	for(i=0;i<3;i++)
		if(_cgtf[0].gtf().l()[i] != t2.cgtf()[0].gtf().l()[i]) return false;
	for(i=0;i<_numberOfFunctions;i++)
	{
		if(fabs(_cgtf[i].gtf().exposant()-t2.cgtf().gtf()[i].exposant())>1e-10) return false;
		if(fabs(_cgtf[i].gtf().coefficient()-t2.cgtf()[i].gtf().coefficient())>1e-10) return false;
		for(c=0;c<3;c++)
			if(fabs(_cgtf[i].gtf().coord()[c]-t2.cgtf()[i].gtf().coord()[c])>1e-10) return false;
	}
	return true;
}
*/
bool operator==(LCAO a, LCAO b)
{
	size_t i,j;
	size_t c=0;
	if(a.cgtf().size()!=b.cgtf().size())
		return false;
	for(i=0; i<a.cgtf().size(); i++)
		for(j=0; j<a.cgtf().size(); j++)
			if(a.cgtf()[i]==b.cgtf()[j])
				c++;
	if(c==a.cgtf().size())
		return true;
	else
		return false;
}
