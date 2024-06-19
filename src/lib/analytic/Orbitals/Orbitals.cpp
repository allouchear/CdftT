#include<iostream>
#include<iomanip>
#include<analytic/Orbitals/Orbitals.h>
#include<analytic/Utils/LM.h>

using namespace std;

Orbitals::Orbitals()
{
	_all_f = vector<vector<double>> ();
	_lcao = LCAO ();
	_symbol = vector<string> ();
	_orbital_energy = vector<vector<double>> ();
	_primitive_centers = vector<int> ();
	_atomic_numbers = vector<int> ();
	_numOrb = vector<int> ();
	_numberOfFunctions=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_number_of_atoms=0;
	_occupation_number=vector<vector<double>> ();
}

Orbitals::Orbitals(WFX& wfx, Binomial& Bin, const PeriodicTable& Table)
{
	int i,j;

	GTF gtf;
	vector<CGTF> vcgtf (wfx.Number_of_Primitives());
	_lcao = LCAO ();
	vector<vector<double>> Coord(wfx.Number_of_Nuclei(), vector<double> (0));
	
	for(i=0; i<wfx.Number_of_Nuclei(); i++)
		for(j=i*3; j<3*(1+i); j++)
			Coord[i].push_back(wfx.Nuclear_Cartesian_Coordinates()[j]);

	gtf=GTF();
	for(j=0; j<wfx.Number_of_Primitives(); j++)
	{			
			gtf.push_back(wfx.Primitive_Exponents()[j], 1.0, Coord[wfx.Primitive_Centers()[j]-1], setLxyz(wfx.Primitive_Types()[j]), Bin);
			vcgtf[j].push_back(gtf);
	}
	_lcao.push_back(vcgtf);

	_coefficients=vector<vector<vector<double>>> (2);
	_coefficients[0]=vector<vector<double>> (wfx.Molecular_Orbital_Primitive_Coefficients()[0].size());
	_coefficients[1]=vector<vector<double>> (wfx.Molecular_Orbital_Primitive_Coefficients()[1].size());

	for(size_t i=0; i<wfx.Molecular_Orbital_Primitive_Coefficients()[0].size(); i++)
		_coefficients[0][i]=wfx.Molecular_Orbital_Primitive_Coefficients()[0][i].Coefficients();

	for(size_t i=0; i<wfx.Molecular_Orbital_Primitive_Coefficients()[1].size(); i++)
		_coefficients[1][i]=wfx.Molecular_Orbital_Primitive_Coefficients()[1][i].Coefficients();

	_primitive_centers=wfx.Primitive_Centers();
	_atomic_numbers=wfx.Atomic_Number();
	_numberOfFunctions=wfx.Number_of_Occupied_Molecular_Orbital();
	if(!wfx.AlphaAndBeta())
		_numberOfFunctions/=2;
	_number_of_alpha_electrons=wfx.Number_of_Alpha_Electrons();
	_number_of_beta_electrons=wfx.Number_of_Beta_Electrons();
	_number_of_atoms=wfx.Number_of_Nuclei();
	_orbital_energy=wfx.Molecular_Orbital_Energies();
	_symbol=wfx.Nuclear_Names();
	_numOrb=vector<int> (2,0);
	_occupation_number=wfx.Molecular_Orbital_Occupation_Numbers();
	_descriptors=Descriptors(wfx, Table);
}

Orbitals::Orbitals(FCHK& fchk, Binomial& Bin, const PeriodicTable& Table)
{
	_numberOfFunctions=fchk.NumberOfBasisFunctions();
	_lcao = LCAO (fchk, Bin);
	_number_of_alpha_electrons=fchk.NumberOfAlphaElectrons();
	_number_of_beta_electrons=fchk.NumberOfBetaElectrons();
	_number_of_atoms=fchk.NumberOfAtoms();
	_coefficients=vector<vector<vector<double>>> (2, vector<vector<double>> ());

	int nOrb_alpha=fchk.AlphaMOCoefficients().size()/fchk.AlphaOrbitalEnergies().size();
	int nOrb_beta=fchk.BetaMOCoefficients().size()/fchk.BetaOrbitalEnergies().size();
	int np=fchk.NumberOfPrimitivesShells();

	for(int i=0; i<nOrb_alpha; i++)
		for(int j=0; j<np; j++)
			_coefficients[0][i][j]=fchk.AlphaMOCoefficients()[i*np+j];
	for(int i=0; i<nOrb_beta; i++)
		for(int j=0; j<np; j++)		
			_coefficients[1][i][j]=fchk.BetaMOCoefficients()[i*np+j];

	_orbital_energy=vector<vector<double>> (2, vector<double> ());
	_orbital_energy[0]=fchk.AlphaOrbitalEnergies();
	_orbital_energy[1]=fchk.BetaOrbitalEnergies();
	_occupation_number=vector<vector<double>> (2);
	_occupation_number[0]=fchk.AlphaOccupation();
	_occupation_number[1]=fchk.BetaOccupation();
	_atomic_numbers=fchk.AtomicNumbers();
	_symbol=vector<string> (_number_of_atoms);
	
	for(int i=0; i<_number_of_atoms; i++)
		_symbol[i]=Table.element(_atomic_numbers[i]).name();

	_numOrb=vector<int> (2,0);
	_descriptors=Descriptors(fchk, Table);
	_primitive_centers=vector<int> ();
}

double Orbitals::Overlap(int i, int j, int alpha)
{
	double sum=0.0;

	for(size_t m=0; m<_coefficients[alpha][i].size(); m++)
		for(size_t n=0; n<_coefficients[alpha][j].size(); n++)
			sum+=_coefficients[alpha][i][m]*_coefficients[alpha][j][n]*_lcao.overlapLCAO();

	return sum;
}

void Orbitals::PrintOverlap(int i, int j, int alpha)
{
	cout<<"OverlapLCAO <"<<i<<"|"<<j<<"> = "<<Overlap(i, j, alpha)<<endl;
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

vector<vector<double>> Orbitals::get_S() const
{
	int i,j,n;
	vector<CGTF> cgtf =_lcao.cgtf();
	n=cgtf.size();
	vector<vector<double>> S (n, vector<double> (n,0.0));
	

#ifdef ENABLE_OMP
#pragma omp parallel for private(i,j)
#endif
	for(i=0; i<n; i++)
		for(j=i; j<n; j++)
			S[i][j]=S[j][i]=cgtf[i].overlapCGTF(cgtf[j]);

	return S;
}

vector<double> Orbitals::get_f(int orb, int alpha) const
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