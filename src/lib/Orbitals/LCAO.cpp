#include<iostream>
#include<iomanip>
#include <Utils/LM.h>
#include <Orbitals/LCAO.h>

using namespace std;

LCAO::LCAO()
{
	_cgtf=vector<CGTF> ();
	_numberOfFunctions=0;
	_bino=Binomial();
}

LCAO::LCAO(vector<CGTF> A) : _cgtf(A)
{
	_numberOfFunctions=_cgtf.size();
	_bino=A[0].bino();
}

LCAO::LCAO(FCHK& fchk, Binomial& Bin)
{
	_bino=Bin;
	int lmax = fchk.HighestAngularMomentum();
	int nShells = fchk.NumberOfContractedShells();
	int llmax = (lmax+1)*(lmax+2)/2;
	vector<int> numAtoms = fchk.ShellToAtomMap();
	vector<int> nPrimitivesByShell = fchk.NumberOfPrimitivesPerShell();
	vector<int> nCoefs (llmax);
	vector<int> shellTypes = fchk.ShellTypes();
	vector<double> contractionsCoefsSP = fchk.spContractionCoefficients();
	vector<double> coordinatesForShells = fchk.CoordinatesForShells();
	vector<double> primitiveExponents = fchk.PrimitiveExponents();
	vector<vector<double>> coefs (llmax, vector<double> (llmax));
	vector<vector<vector<int>>> l (3, vector<vector<int>> (llmax, vector<int> (llmax)));

	int NOrb = 0;
	for(int nS=0;nS<nShells;nS++) 
	{
		if(shellTypes[nS]<-1)
			NOrb += 2*abs(shellTypes[nS])+1; /* Spherical D, F, G, ...*/
		else if(shellTypes[nS]==-1)
			NOrb +=  4; /* This a SP.*/
		else
			NOrb +=  (shellTypes[nS]+1)*(shellTypes[nS]+2)/2; /* Cartesian S,P,D,F,G,..*/
	}

	vector<CGTF> temp (NOrb);
	int kOrb = 0;
	int kPrimitive = 0;
	for(int nS = 0;nS<nShells; nS++)
	{
		int nM = 0;
		/* printf("begin primitive nS = %d\n",nS);*/
		if(shellTypes[nS]<-1)
			nM = 2*abs(shellTypes[nS])+1; /* Sperical D, F, G, ...*/
		else if(shellTypes[nS]==-1)
			nM = 1; /* This a SP. Make S before */
		else
			nM = (shellTypes[nS]+1)*(shellTypes[nS]+2)/2;

		/* printf("nM = %d\n",nM);*/
		if(shellTypes[nS]==-1)
			getlTable(0, nCoefs, coefs, l, _bino); /* This a SP. Make S before */
		else
			getlTable(shellTypes[nS], nCoefs, coefs, l, _bino); 
		/* printf("end getlTable\n");*/
		for(int m=0;m<nM;m++)
		{
			int ip,j,n;
			/* printf("P : m = %d nCoef = %d nPrim = %d\n",m,nCoefs[m],nPrimitivesByShell[nS]);*/
			temp[kOrb]= CGTF ();
			temp[kOrb].setNumCenter(numAtoms[nS]-1);
  			j = -1;
	 		for(ip=0;ip<nPrimitivesByShell[nS];ip++)
 				for(n=0;n<nCoefs[m];n++)
	 			{
	 		   		j++;
	   				vector<double> coord_ (3);
	   				vector<int> l_ (3);
	   				for(int c=0;c<3;c++)
	   				{
	   					coord_[c] = coordinatesForShells[c+nS*3];
						l_[c] = l[c][m][n];
	   				}
	   				GTF gtf (primitiveExponents[kPrimitive+ip], contractionsCoefsSP[kPrimitive+ip]*coefs[m][n], coord_, l_, _bino);
	   				temp[kOrb].push_back(gtf);
	 			}
			kOrb++;
		}
		if(shellTypes[nS]==-1) /* This a SP. Now make P*/
		{
			getlTable(-1, nCoefs, coefs, l, _bino);
			nM = 3;
			for(int m=0;m<nM;m++)
			{
				int ip,j,n;
				/* printf("P : m = %d nCoef = %d nPrim = %d\n",m,nCoefs[m],nPrimitivesByShell[nS]);*/
				temp[kOrb]= CGTF ();
				temp[kOrb].setNumCenter(numAtoms[nS]-1);
          			j = -1;
	 			for(ip=0;ip<nPrimitivesByShell[nS];ip++)
 					for(n=0;n<nCoefs[m];n++)
	 				{
	 		   			j++;
	   					vector<double> coord_ (3);
	   					vector<int> l_ (3);
	   					for(int c=0;c<3;c++)
	   					{
	   						coord_[c] = coordinatesForShells[c+nS*3];
	   						l_[c] = l[c][m][n];
	   					}
	   				GTF gtf (primitiveExponents[kPrimitive+ip], contractionsCoefsSP[kPrimitive+ip]*coefs[m][n], coord_, l_, _bino);
	   				temp[kOrb].push_back(gtf);
	 			}
				kOrb++;
			}
		}
		/* printf("end primitive nS = %d\n",nS);*/
		kPrimitive += nPrimitivesByShell[nS];
	}
	_cgtf=temp;
	_numberOfFunctions=NOrb;

	if(_cgtf.size() != (size_t)NOrb)
	{
		cout<<"Error : Their is "<<_cgtf.size()<<" CGTF for "<<NOrb<<" basis in file."<<endl;
		cout<<"Please check your file."<<endl;
		exit(1);
	}
}

double LCAO::ERILCAO(LCAO& q, LCAO& r, LCAO& s)
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
	cout<<"Norme LCAO = "<<sum<<endl;
}

double LCAO::overlapLCAO()
{
	double sum=0.0;
	int n;
	int np;

#ifdef ENABLE_OMP
#pragma omp parallel for private(n,np) reduction(+:sum)
#endif
	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<_numberOfFunctions;np++)
			sum += _cgtf[n].overlapCGTF(_cgtf[np]);

	//cout<<"Test Overlap in LCAO "<<endl;

	return sum;
}

double LCAO::overlap3LCAO(const vector<double>& c1, const vector<double>& c2, const vector<double>& c3)
{
	double sum=0.0;
	int n;
	int np;
	int ns;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<_numberOfFunctions;np++)
			for(ns=0;ns<_numberOfFunctions;ns++)
				sum += c1[n]*c2[np]*c3[ns]*_cgtf[n].overlap3CGTF(_cgtf[np],_cgtf[ns]);

	return sum;
}

double LCAO::overlap4LCAO(const vector<double>& c1, const vector<double>& c2, const vector<double>& c3, const vector<double>& c4)
{
	double sum=0.0;
	int np;
	int nq;
	int nr;
	int ns;

	for(np=0;np<_numberOfFunctions;np++)
		for(nq=0;nq<_numberOfFunctions;nq++)
			for(nr=0;nr<_numberOfFunctions;nr++)
				for(ns=0;ns<_numberOfFunctions;ns++)
					sum += c1[np]*c2[nq]*c3[nr]*c4[ns]*_cgtf[np].overlap4CGTF(_cgtf[nq],_cgtf[nr],_cgtf[ns]);

	return sum;
}

double LCAO::kineticLCAO(LCAO& right)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfFunctions;n++)
		for(np=0;np<numberOfFunctions();np++)
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

double LCAO::func(const vector<double>& c, double x, double y, double z) const
{
	double r=0.0;
	if(c.size()!=_cgtf.size())
	{
		cout<<"Error, their is "<<c.size()<<" coefficients for "<<_cgtf.size()<<" CGTF."<<endl;
		cout<<"Please, check the code or your file !"<<endl;
		exit(1);
	}
	for(size_t i=0; i<c.size(); i++)
	{
		if(abs(c[i])>1e-10)
			r+=c[i] * _cgtf[i].func(x,y,z);
	}

	return r;
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
		if(fabs(_cgtf[i].gtf().coefficients()-t2.cgtf()[i].gtf().coefficients())>1e-10) return false;
		for(c=0;c<3;c++)
			if(fabs(_cgtf[i].gtf().coord()[c]-t2.cgtf()[i].gtf().coord()[c])>1e-10) return false;
	}
	return true;
}
*/

void LCAO::push_back(const vector<CGTF>& cgtf)
{
	_cgtf=cgtf;
	_numberOfFunctions=_cgtf.size();

}
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

ostream& operator<<(ostream& flux, const LCAO& lcao)
{
	flux<<std::scientific;
	flux<<std::setprecision(10);
	flux<<std::setw(20);
	flux<<left<<setw(20)<<"CoefLCAO"<<setw(20)<<"Exp"<<setw(5)<<"Lx"<<setw(5)<<"Ly"<<setw(5)
	<<"Lz"<<setw(20)<<"x"<<setw(20)<<"y"<<setw(20)<<"z"<<endl;
	for(int i=0; i<lcao.numberOfFunctions(); i++)
		flux<<left<<lcao.cgtf()[i]<<endl;

	return flux;
}

double operator*(const LCAO& a, const vector<vector<double>>& m)
{
	double r=1.0;
	for(size_t i=1; i<m.size(); i++)
		r*=a.func(m[i],m[0][0],m[0][1],m[0][2]);
	
	return r;
}