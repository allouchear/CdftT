#include<iostream>
#include<iomanip>
#include<analytic/Orbitals/Orbitals.h>
#include<analytic/Utils/LM.h>

using namespace std;

Orbitals::Orbitals()
{
	_all_f = vector<vector<double>> ();
	_vcgtf = vector<CGTF> ();
	_symbol = vector<string> ();
	_orbital_energy = vector<vector<double>> ();
	_primitive_centers = vector<int> ();
	_atomic_numbers = vector<int> ();
	_numOrb = vector<int> ();
	_numberOfAo=0;
	_numberOfMo=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_number_of_atoms=0;
	_occupation_number=vector<vector<double>> ();
	_alpha_and_beta=false;
	_bino=Binomial ();
}

Orbitals::Orbitals(WFX& wfx, Binomial& Bin, const PeriodicTable& Table)
{
	int i,j;
	_bino=Bin;
	GTF gtf;
	_vcgtf = vector<CGTF> (wfx.Number_of_Primitives());
	vector<vector<double>> Coord(wfx.Number_of_Nuclei(), vector<double> (0));
	
	for(i=0; i<wfx.Number_of_Nuclei(); i++)
		for(j=i*3; j<3*(1+i); j++)
			Coord[i].push_back(wfx.Nuclear_Cartesian_Coordinates()[j]);

	gtf=GTF();
	for(j=0; j<wfx.Number_of_Primitives(); j++)
	{			
			gtf.push_back(wfx.Primitive_Exponents()[j], 1.0, Coord[wfx.Primitive_Centers()[j]-1], setLxyz(wfx.Primitive_Types()[j]), Bin);
			_vcgtf[j].push_back(gtf);
			_vcgtf[j].setCoef(1.0);
	}

	_coefficients=vector<vector<vector<double>>> (2);
	_coefficients[0]=vector<vector<double>> (wfx.Molecular_Orbital_Primitive_Coefficients()[0].size());
	_coefficients[1]=vector<vector<double>> (wfx.Molecular_Orbital_Primitive_Coefficients()[1].size());
	for(size_t i=0; i<wfx.Molecular_Orbital_Primitive_Coefficients()[0].size(); i++)
		_coefficients[0][i]=wfx.Molecular_Orbital_Primitive_Coefficients()[0][i].Coefficients();
	for(size_t i=0; i<wfx.Molecular_Orbital_Primitive_Coefficients()[1].size(); i++)
		_coefficients[1][i]=wfx.Molecular_Orbital_Primitive_Coefficients()[1][i].Coefficients();

	_primitive_centers=wfx.Primitive_Centers();
	_atomic_numbers=wfx.Atomic_Number();
	_numberOfMo=wfx.Number_of_Occupied_Molecular_Orbital();
	if(!wfx.AlphaAndBeta())
		_numberOfMo/=2;
	_number_of_alpha_electrons=wfx.Number_of_Alpha_Electrons();
	_number_of_beta_electrons=wfx.Number_of_Beta_Electrons();
	_number_of_atoms=wfx.Number_of_Nuclei();
	_orbital_energy=wfx.Molecular_Orbital_Energies();
	_symbol=wfx.Nuclear_Names();
	_numOrb=vector<int> (2,0);
	_occupation_number=wfx.Molecular_Orbital_Occupation_Numbers();
	_alpha_and_beta=wfx.AlphaAndBeta();
	_descriptors=Descriptors(wfx, Table);
	_numberOfAo=_vcgtf.size();
}

Orbitals::Orbitals(FCHK& fchk, Binomial& Bin, const PeriodicTable& Table)
{
	_numberOfMo=fchk.NumberOfBasisFunctions();

	_bino=Bin;
	int lmax = fchk.HighestAngularMomentum();
	int nShells = fchk.NumberOfContractedShells();
	int llmax = (lmax+1)*(lmax+2)/2;
	vector<int> numAtoms = fchk.ShellToAtomMap();
	vector<int> nPrimitivesByShell = fchk.NumberOfPrimitivesPerShell();
	vector<int> nCoefs (llmax);
	vector<int> shellTypes = fchk.ShellTypes();
	vector<double> contractionsCoefs = fchk.ContractionCoefficients();
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
	_vcgtf = vector<CGTF> (NOrb);
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
			_vcgtf[kOrb]= CGTF ();
			_vcgtf[kOrb].setNumCenter(numAtoms[nS]-1);
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
	   				GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
	   				_vcgtf[kOrb].push_back(gtf);
	 				_vcgtf[kOrb].setCoef(contractionsCoefs[kPrimitive+ip]*coefs[m][n]);
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
				_vcgtf[kOrb]= CGTF ();
				_vcgtf[kOrb].setNumCenter(numAtoms[nS]-1);
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
	   				GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
	   				_vcgtf[kOrb].push_back(gtf);
	   				_vcgtf[kOrb].setCoef(contractionsCoefsSP[kPrimitive+ip]*coefs[m][n]);
	 			}
				kOrb++;
			}
		}
		/* printf("end primitive nS = %d\n",nS);*/
		kPrimitive += nPrimitivesByShell[nS];
	}

	_numberOfAo=_vcgtf.size();

	if(_numberOfAo != _numberOfMo)
	{
		cout<<"Error : Their is "<<_vcgtf.size()<<" CGTF for "<<_numberOfMo<<" basis in file."<<endl;
		cout<<"Please check your file."<<endl;
		exit(1);
	}

	_number_of_alpha_electrons=fchk.NumberOfAlphaElectrons();
	_number_of_beta_electrons=fchk.NumberOfBetaElectrons();
	_number_of_atoms=fchk.NumberOfAtoms();
	_coefficients=vector<vector<vector<double>>> (2, vector<vector<double>> ());
	int nOrb_alpha=fchk.AlphaMOCoefficients().size()/fchk.AlphaOrbitalEnergies().size();
	int nOrb_beta=fchk.BetaMOCoefficients().size()/fchk.BetaOrbitalEnergies().size();

	_coefficients[0]=vector<vector<double>> (nOrb_alpha, vector<double> (_numberOfMo));
	_coefficients[1]=vector<vector<double>> (nOrb_beta, vector<double> (_numberOfMo));

	for(int i=0; i<nOrb_alpha; i++)
		for(int j=0; j<_numberOfMo; j++)
			_coefficients[0][i][j]=fchk.AlphaMOCoefficients()[i*_numberOfMo+j];
	for(int i=0; i<nOrb_beta; i++)
		for(int j=0; j<_numberOfMo; j++)		
			_coefficients[1][i][j]=fchk.BetaMOCoefficients()[i*_numberOfMo+j];
	_orbital_energy=vector<vector<double>> (2, vector<double> ());
	_orbital_energy[0]=fchk.AlphaOrbitalEnergies();
	_orbital_energy[1]=fchk.BetaOrbitalEnergies();

	_occupation_number=vector<vector<double>> (2);
	_alpha_and_beta=fchk.AlphaAndBeta();

	if(_alpha_and_beta)
	{
		_occupation_number[0]=vector<double> (_numberOfMo);
		for(int i=0; i<_numberOfMo; i++)
			_occupation_number[0][i]=fchk.AlphaOccupation()[i]+fchk.BetaOccupation()[i];
	}
	else
	{
		_occupation_number[0]=fchk.AlphaOccupation();
		_occupation_number[1]=fchk.BetaOccupation();
	}

	_atomic_numbers=fchk.AtomicNumbers();
	_symbol=vector<string> (_number_of_atoms);
	
	for(int i=0; i<_number_of_atoms; i++)
		_symbol[i]=Table.element(_atomic_numbers[i]).name();

	_numOrb=vector<int> (2,0);
	_descriptors=Descriptors(fchk, Table);
	_primitive_centers=vector<int> ();
	NormaliseAllBasis();
}

double Orbitals::ERIorbitals(Orbitals& q, Orbitals& r, Orbitals& s)
{
	int np,nq;
	int nr,ns;
	double sum = 0.0;

	for(np=0;np<_numberOfAo;np++)
		for(nq=0;nq<q._numberOfAo;nq++)
			for(nr=0;nr<r._numberOfAo;nr++)
				for(ns=0;ns<s._numberOfAo;ns++)
					sum += _vcgtf[np].ERICGTF(q.vcgtf()[nq],r.vcgtf()[nr],s.vcgtf()[ns]); 

	return sum;
}

double Orbitals::Overlap(int i, int j, int alpha)
{
	double sum=0.0;

#ifdef ENABLE_OMP
#pragma omp parallel for reduction(+:sum)
#endif
	for(size_t m=0; m<_coefficients[alpha][i].size(); m++)
		for(size_t n=0; n<_coefficients[alpha][j].size(); n++)
			sum+=_coefficients[alpha][i][m]*_coefficients[alpha][j][n]*_vcgtf[m].overlapCGTF(_vcgtf[n]);

	return sum;
}

void Orbitals::PrintOverlap(int i, int j, int alpha)
{
	cout<<"Overlap <"<<i<<"|"<<j<<"> = "<<Overlap(i, j, alpha)<<endl;
}

double Orbitals::Overlap3Orbitals(int i, int j, int k, int alpha)
{
	double sum=0.0;
	int n;
	int np;
	int ns;

	for(n=0;n<_numberOfAo;n++)
		for(np=0;np<_numberOfAo;np++)
			for(ns=0;ns<_numberOfAo;ns++)
				sum += _coefficients[alpha][i][n]*_coefficients[alpha][j][np]*_coefficients[alpha][k][ns]*_vcgtf[n].overlap3CGTF(_vcgtf[np],_vcgtf[ns]);

	return sum;
}

double Orbitals::Overlap4Orbitals(int i, int j, int k, int l, int alpha)
{
	double sum=0.0;
	int np;
	int nq;
	int nr;
	int ns;

	for(np=0;np<_numberOfAo;np++)
		for(nq=0;nq<_numberOfAo;nq++)
			for(nr=0;nr<_numberOfAo;nr++)
				for(ns=0;ns<_numberOfAo;ns++)
					sum += _coefficients[alpha][i][np]*_coefficients[alpha][j][nq]*_coefficients[alpha][k][nr]*_coefficients[alpha][l][ns]*_vcgtf[np].overlap4CGTF(_vcgtf[nq],_vcgtf[nr],_vcgtf[ns]);

	return sum;
}

double Orbitals::kinetic()
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfAo;n++)
		for(np=0;np<_numberOfAo;np++)
			sum += _vcgtf[n].kineticCGTF(_vcgtf[np]);


	return sum;
}

double Orbitals::ionicPotential(vector<double> C, double Z)
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfAo;n++)
		for(np=0;np<_numberOfAo;np++)
			sum += _vcgtf[n].ionicPotentialCGTF(_vcgtf[np], C, Z); 

	return sum;
}

double Orbitals::OrbstarOrb()
{
	int n;
	int np;
	double sum=0.0;

	for(n=0;n<_numberOfAo;n++)
		for(np=0;np<_numberOfAo;np++)
			sum += _vcgtf[n].CGTFstarCGTF(_vcgtf[np]);

	return sum;
}

double Orbitals::OrbxyzOrb(int ix, int iy, int iz)
{
	double sum=0.0;
	int n;
	int ns;
	vector<double> C(3,0);
	vector<int> l {ix, iy, iz};
	GTF m1(0.0, 1.0, C, l, _bino);
	vector<GTF> mbis (1,m1);
	CGTF m2(mbis);

	for(n=0;n<_numberOfAo;n++)
		for(ns=0;ns<_numberOfAo;ns++)
				sum += _vcgtf[n].gtf()[ns].overlap3GTF(m2.gtf()[0],_vcgtf[n].gtf()[ns]);

	return sum;
}

void Orbitals::NormaliseAllBasis()
{
	int k;

	for(k=0;k<_numberOfAo;k++)
		_vcgtf[k].normaliseCGTF();
}

double Orbitals::func(double x, double y, double z) const
{
	double r=0.0;
	int n;

	if(_alpha_and_beta)
		n=1;
	else
		n=2;

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<_numberOfMo; j++)
		{
			if(_coefficients[i][j].size()!=_vcgtf.size())
			{
				cout<<"Error, their is "<<_coefficients[i][j].size()<<" coefficients for "<<_vcgtf.size()<<" CGTF."<<endl;
				cout<<"Please, check the code or your file !"<<endl;
				exit(1);
			}
			for(int k=0; k<_numberOfMo; k++)
			{
				if(abs(_coefficients[i][j][k])>1e-10)
					r+=_coefficients[i][j][k] * _vcgtf[k].func(x,y,z);
			}
		}
	}

	return r;
}

void Orbitals::HOMO()
{
	if(_number_of_alpha_electrons>=_number_of_beta_electrons)
		_numOrb[0] = _number_of_alpha_electrons-1;
	else
		_numOrb[0] = _number_of_beta_electrons-1;
}

void Orbitals::LUMO()
{
	_numOrb[1] = _numOrb[0]+1;
}

vector<vector<double>> Orbitals::get_S()
{
	int i,j;
	vector<vector<double>> S (_numberOfAo, vector<double> (_numberOfAo,0.0));

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,j)
#endif
	for(i=0; i<_numberOfAo; i++)
		for(j=i; j<_numberOfAo; j++)
			S[i][j]=S[j][i]=_vcgtf[i].overlapCGTF(_vcgtf[j]);

	return S;
}

vector<double> Orbitals::get_f(int orb, int alpha)
{
	int i;
	vector<vector<double>> S =get_S();
	size_t nu,xi;
	vector<double> f(_number_of_atoms,0.0);

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,nu,xi)
#endif
	for(i=0; i<_number_of_atoms; i++)
		for(nu=0; nu<_coefficients[alpha][orb].size(); nu++)
		{
			if(i+1 == _primitive_centers[nu])
				f[i]+=_coefficients[alpha][orb][nu]*_coefficients[alpha][orb][nu];

			for(xi=0; xi<_coefficients[alpha][orb].size(); xi++)
				if(xi!=nu && i+1 == _primitive_centers[nu])
					f[i]+=_coefficients[alpha][orb][xi]*_coefficients[alpha][orb][nu];
		}

	return f;
}

void Orbitals::get_f(int alpha)
{
	vector<vector<double>> S=get_S();
	int i;
	size_t j,nu,xi;
	vector<double> V(_number_of_atoms,0.0);
	vector<vector<double>> f(_numOrb.size(), V);

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,j,nu,xi)
#endif
	for(i=0; i<_number_of_atoms; i++)
		for(j=0; j<_numOrb.size(); j++)
			for(nu=0; nu<_coefficients[alpha][_numOrb[j]].size(); nu++)
			{
				if(i+1 == _primitive_centers[nu])
					f[j][i]+=_coefficients[alpha][_numOrb[j]][nu]*_coefficients[alpha][_numOrb[j]][nu];

				for(xi=0; xi<_coefficients[alpha][_numOrb[j]].size(); xi++)
					if(xi!=nu && i+1 == _primitive_centers[nu])
						f[j][i]+=_coefficients[alpha][_numOrb[j]][xi]*_coefficients[alpha][_numOrb[j]][nu]*S[xi][nu];
			}

	_all_f = f;
}

void Orbitals::HOMO_LUMO()
{
	HOMO();
	LUMO();
}

void Orbitals::HOMO_LUMO(int i, int j)
{
	_numOrb[0]=i;
	_numOrb[1]=j;
}

void Orbitals::PrintDescriptors()
{
	HOMO_LUMO();
	get_f();
	_descriptors.set_mu_fk_data(_all_f, eHOMO(), eLUMO());
	_descriptors.compute_all();
	cout<<_descriptors<<endl;
}

void Orbitals::PrintDescriptors(int i, int j)
{
	HOMO_LUMO(i,j);
	get_f();
	_descriptors.set_mu_fk_data(_all_f, eHOMO(), eLUMO());
	_descriptors.compute_all();
	cout<<_descriptors<<endl;
}

double operator*(const Orbitals& a, const vector<double>& coord)
{
	double r=1.0;
	for(size_t i=1; i<coord.size(); i++)
		r*=a.func(coord[0],coord[1],coord[2]);
	
	return r;
}

ostream& operator<<(ostream& flux, Orbitals& Orb)
{
	flux<<scientific;
	flux<<setprecision(10);
	flux<<setw(20);
	flux<<left<<setw(20)<<"Coef CGTF"<<setw(20)<<"Coef GTF"<<setw(20)<<"Exp"<<setw(5)<<"Lx"<<setw(5)<<"Ly"<<setw(5)
	<<"Lz"<<setw(20)<<"x"<<setw(20)<<"y"<<setw(20)<<"z"<<endl;
	for(int i=0; i<Orb.NumberOfAo(); i++)
		flux<<left<<Orb.vcgtf()[i]<<endl;

	return flux;
}