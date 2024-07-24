#include<iostream>
#include<iomanip>
#include <Orbitals/Orbitals.h>
#include <Utils/LM.h>

using namespace std;

Orbitals::Orbitals()
{
	_struct=Structure();
	_all_f = vector<vector<double>> ();
	_vcgtf = vector<CGTF> ();
	_vcgtf_non_normalise = vector<CGTF> ();
	_symbol = vector<string> ();
	_orbital_energy = vector<vector<double>> ();
	_primitive_centers = vector<int> ();
	_atomic_numbers = vector<int> ();
	_numOrb = vector<int> ();
	_numberOfAo=0;
	_numberOfMo=0;
	_number_of_gtf=0;
	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;
	_number_of_atoms=0;
	_energy=0.0;
	_occupation_number=vector<vector<double>> ();
	_alpha_and_beta=false;
	_bino=Binomial ();
}

Orbitals::Orbitals(WFX& wfx, Binomial& Bin, const PeriodicTable& Table)
{
	_struct=Structure(wfx, Table);
	int i,j;
	_bino=Bin;
	_coordinates=wfx.Nuclear_Cartesian_Coordinates();
	GTF gtf;
	_vcgtf = vector<CGTF> (wfx.Number_of_Primitives());
	_vcgtf_non_normalise = vector<CGTF> ();
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
			_vcgtf[j].setFactorCoef(1.0);
			_vcgtf[j].setNumCenter(wfx.Primitive_Centers()[j]);
			_vcgtf[j].setLtype(getLType(_vcgtf[j].gtf()[0].l()));
			_vcgtf[j].setFormat("Cart");
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
	_number_of_gtf=wfx.Number_of_Primitives();
	if(!wfx.AlphaAndBeta())
		_numberOfMo/=2;
	_number_of_alpha_electrons=wfx.Number_of_Alpha_Electrons();
	_number_of_beta_electrons=wfx.Number_of_Beta_Electrons();
	_number_of_atoms=wfx.Number_of_Nuclei();
	_orbital_energy=wfx.Molecular_Orbital_Energies();
	_symbol=wfx.Nuclear_Names();
	_numOrb=vector<int> (2,0);
	_energy=wfx.Energy();
	_occupation_number=wfx.Molecular_Orbital_Occupation_Numbers();
	_alpha_and_beta=wfx.AlphaAndBeta();
	_descriptors=Descriptors(wfx, Table);
	_numberOfAo=_vcgtf.size();
	_vcgtf_non_normalise=_vcgtf;
}

Orbitals::Orbitals(FCHK& fchk, Binomial& Bin, const PeriodicTable& Table)
{
	_struct=Structure(fchk, Table);
	_vcgtf_non_normalise = vector<CGTF> ();
	_numberOfMo=fchk.NumberOfBasisFunctions();
	_coordinates=fchk.CurrentCartesianCoordinates();
	_bino=Bin;
	_number_of_gtf=0;
	_energy=fchk.TotalEnergy();
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
	string format;

	for(int nS = 0;nS<nShells; nS++)
	{
		int nM = 0;
		/* printf("begin primitive nS = %d\n",nS);*/
		if(shellTypes[nS]<-1)
		{
			nM = 2*abs(shellTypes[nS])+1; /* Sperical D, F, G, ...*/
			format="Sphe";
		}
		else if(shellTypes[nS]==-1)
		{
			nM = 1; /* This a SP. Make S before */
			format="Cart";
		}
		else
		{
			nM = (shellTypes[nS]+1)*(shellTypes[nS]+2)/2;
			format="Cart";
		}

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
			_vcgtf[kOrb].setNumCenter(numAtoms[nS]);
			_vcgtf[kOrb].setFactorCoef(1.0);
			_vcgtf[kOrb].setFormat(format);
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
	 				_vcgtf[kOrb].setLtype(getLType(l_));
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
				_vcgtf[kOrb].setNumCenter(numAtoms[nS]);
				_vcgtf[kOrb].setFactorCoef(1.0);
				_vcgtf[kOrb].setFormat(format);
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
	   				_vcgtf[kOrb].setLtype(getLType(l_));
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
		_occupation_number[1]=_occupation_number[0];
	}
	else
	{
		_occupation_number[0]=fchk.AlphaOccupation();
		_occupation_number[1]=fchk.BetaOccupation();
	}

	_atomic_numbers=fchk.AtomicNumbers();
	_symbol=vector<string> (_number_of_atoms);
	
	for(int i=0; i<_number_of_atoms; i++)
		_symbol[i]=Table.element(_atomic_numbers[i]).symbol();

	_numOrb=vector<int> (2,0);
	_descriptors=Descriptors(fchk, Table);

	for(size_t i=0; i<_vcgtf.size(); i++)
		for(int j=0; j<_vcgtf[i].numberOfFunctions(); j++)
			_primitive_centers.push_back(_vcgtf[i].NumCenter());

	_vcgtf_non_normalise=_vcgtf;

	NormaliseAllBasis();

	for(size_t i=0; i<_vcgtf.size(); i++)
		_number_of_gtf+=_vcgtf[i].numberOfFunctions();
}

Orbitals::Orbitals(MOLDENGAB& moldengab, Binomial& Bin, const PeriodicTable& Table)
{
	_struct=Structure(moldengab, Table);
	_vcgtf_non_normalise = vector<CGTF> ();
	_number_of_gtf=0;
	_energy=0;
	_coefficients=vector<vector<vector<double>>> (2);
	_coefficients[0]=moldengab.AlphaMOCoefs();
	_coefficients[1]=moldengab.BetaMOCoefs();

	_number_of_atoms=moldengab.NumberOfAtoms();
	_coordinates=vector<double> (_number_of_atoms*3);
	for(int i=0; i<_number_of_atoms; i++)
	{
		int k=0;
		for(int j=i*3; j<(i+1)*3; j++)
		{
			_coordinates[j]=moldengab.Coordinates()[i][k];
			k++;
		}
	}
	_numberOfMo=moldengab.NumberOfMO();
	_alpha_and_beta=moldengab.AlphaAndBeta();

	if(!_alpha_and_beta)
		_numberOfMo/=2;

	_number_of_alpha_electrons=0;
	_number_of_beta_electrons=0;

	if(_alpha_and_beta)
	{
		for(size_t i=0; i<moldengab.AlphaOccupation().size(); i++)
			if(moldengab.AlphaOccupation()[i]==1)
				_number_of_alpha_electrons++;
		_number_of_beta_electrons=_number_of_alpha_electrons;
	}

	else
	{
		for(size_t i=0; i<moldengab.AlphaOccupation().size(); i++)
			if(moldengab.AlphaOccupation()[i]==1)
				_number_of_alpha_electrons++;

		for(size_t i=0; i<moldengab.BetaOccupation().size(); i++)
			if(moldengab.BetaOccupation()[i]==1)
				_number_of_beta_electrons++;
	}

	_bino=Bin;
		
	_atomic_numbers=moldengab.AtomicNumbers();
	_symbol=moldengab.Symbol();
	_orbital_energy=vector<vector<double>> (2);
	_orbital_energy[0]=moldengab.AlphaEnergies();
	_orbital_energy[1]=moldengab.BetaEnergies();	
	_numOrb=vector<int> (2,0);
	_occupation_number=vector<vector<double>> (2);
	_occupation_number[0]=moldengab.AlphaOccupation();
	_occupation_number[1]=moldengab.BetaOccupation();
	_descriptors=Descriptors(moldengab, Table);

	vector<int> Ltypes = moldengab.Ltypes();
	int nShells = Ltypes.size();

	int lmax=0;
	for(int i=0; i<nShells; i++)
		if(lmax<abs(Ltypes[i]))
			lmax=abs(Ltypes[i]);

	int llmax = (lmax+1)*(lmax+2)/2;
	vector<int> nPrimitivesByShell = moldengab.NumberOfGtf();
	vector<int> nCoefs (llmax);
	
	vector<double> FactorCoefs = moldengab.FactorCoefficients();
	vector<double> CgtfCoefs = moldengab.CgtfCoefficients();
	vector<vector<double>> coordinatesForShells;
	vector<int> NatBasis = moldengab.NatBasis();

	for(int i=0; i<_number_of_atoms; i++)
		for(int j=0; j<NatBasis[i]; j++)
			coordinatesForShells.push_back(moldengab.Coordinates()[i]);

	vector<double> primitiveExponents = moldengab.Exposants();
	vector<vector<double>> coefs (llmax, vector<double> (llmax));
	vector<vector<vector<int>>> l (3, vector<vector<int>> (llmax, vector<int> (llmax)));

	int NOrb = moldengab.NumberOfMOCoefs();
	_vcgtf = vector<CGTF> (NOrb);
	string format;

	int kOrb = 0;
	int kPrimitive = 0;
	for(int nS = 0;nS<nShells; nS++)
	{
		int nM = 0;

		if(Ltypes[nS]<-1)
		{
			nM = 2*abs(Ltypes[nS])+1; /* Sperical D, F, G, ...*/
			format="Sphe";
		}
		else if(Ltypes[nS]==-1)
		{
			nM = 1; /* This a SP. Make S before */
			format="Cart";
		}
		else
		{
			nM = (Ltypes[nS]+1)*(Ltypes[nS]+2)/2;
			format="Cart";
		}

		if(Ltypes[nS]==-1)
			getlTable(0, nCoefs, coefs, l, _bino); /* This a SP. Make S before */
		
		else
			getlTable(Ltypes[nS], nCoefs, coefs, l, _bino); 
	

		for(int m=0;m<nM;m++)
		{

			int ip,j,n;
			_vcgtf[kOrb]= CGTF ();
			_vcgtf[kOrb].setNumCenter(moldengab.Numcenter()[nS]);
			_vcgtf[kOrb].setLtype(moldengab.ShellTypes()[nS]);
			_vcgtf[kOrb].setFactorCoef(FactorCoefs[nS]);
			_vcgtf[kOrb].setFormat(format);
  			j = -1;

	 		for(ip=0;ip<nPrimitivesByShell[nS];ip++)
 				for(n=0;n<nCoefs[m];n++)
	 			{
	 		   		j++;
	   				vector<double> coord_ (3);
	   				vector<int> l_ (3);
	   				for(int c=0;c<3;c++)
	   				{
	   					coord_[c] = coordinatesForShells[kPrimitive+ip][c];
						l_[c] = l[c][m][n];
	   				}
	   				GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
	   				_vcgtf[kOrb].push_back(gtf);
	 				_vcgtf[kOrb].setCoef(FactorCoefs[kPrimitive+ip]*CgtfCoefs[kPrimitive+ip]*coefs[m][n]);
	 			}
			kOrb++;
		}
		if(Ltypes[nS]==-1) /* This a SP. Now make P*/
		{
			getlTable(-1, nCoefs, coefs, l, _bino);
			nM = 3;
			for(int m=0;m<nM;m++)
			{
				int ip,j,n;

				_vcgtf[kOrb]= CGTF ();
				_vcgtf[kOrb].setNumCenter(moldengab.Numcenter()[nS]);
				_vcgtf[kOrb].setLtype(moldengab.ShellTypes()[nS]);
				_vcgtf[kOrb].setFactorCoef(FactorCoefs[nS]);
				_vcgtf[kOrb].setFormat(format);
      			j = -1;
	 			for(ip=0;ip<nPrimitivesByShell[nS];ip++)
 					for(n=0;n<nCoefs[m];n++)
	 				{
	 		   			j++;
	   					vector<double> coord_ (3);
	   					vector<int> l_ (3);
	   					for(int c=0;c<3;c++)
	   					{
	   						coord_[c] = coordinatesForShells[nS][c];
	   						l_[c] = l[c][m][n];
	   					}
	   				GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
	   				_vcgtf[kOrb].push_back(gtf);
	   				_vcgtf[kOrb].setCoef(FactorCoefs[kOrb]*CgtfCoefs[kPrimitive+ip]*coefs[m][n]);
	 				}
				kOrb++;
			}
		}
		kPrimitive += nPrimitivesByShell[nS];
	}

	_numberOfAo=_vcgtf.size();

	if(_numberOfAo != moldengab.NumberOfMOCoefs())
	{
		cout<<"Error : Their is "<<_vcgtf.size()<<" CGTF for "<<moldengab.NumberOfMOCoefs()<<" basis in file."<<endl;
		cout<<"Please check your file."<<endl;
		exit(1);
	}

	_vcgtf_non_normalise=_vcgtf;

	for(size_t i=0; i<_vcgtf.size(); i++)
		for(int j=0; j<_vcgtf[i].numberOfFunctions(); j++)
			_primitive_centers.push_back(_vcgtf[i].NumCenter());

	NormaliseAllBasis();

	for(size_t i=0; i<_vcgtf.size(); i++)
		_number_of_gtf+=_vcgtf[i].numberOfFunctions();

	Sorting();
}

Orbitals::Orbitals(LOG& log, Binomial& Bin, const PeriodicTable& Table)
{
	_struct=Structure(log, Table);
	_vcgtf_non_normalise = vector<CGTF> ();
	_number_of_gtf=0;
	_coefficients=vector<vector<vector<double>>> (2);
	_coefficients[0]=log.AlphaMOcoefs();
	_coefficients[1]=log.BetaMOcoefs();
	_energy=log.Energy();
	_number_of_alpha_electrons=log.NumberOfAlphaElectrons();
	_number_of_beta_electrons=log.NumberOfBetaElectrons();

	_number_of_atoms=log.NumberOfAtoms();
	_coordinates=vector<double> (_number_of_atoms*3);
	for(int i=0; i<_number_of_atoms; i++)
	{
		int k=0;
		for(int j=i*3; j<(i+1)*3; j++)
		{
			_coordinates[j]=log.Coordinates()[i][k];
			k++;
		}
	}
	_numberOfMo=log.NumberOfMO();
	_alpha_and_beta=log.AlphaAndBeta();

	_bino=Bin;
		
	_atomic_numbers=log.AtomicNumbers();
	_symbol=log.Symbol();
	_orbital_energy=vector<vector<double>> (2);
	_orbital_energy[0]=log.AlphaEnergy();
	_orbital_energy[1]=log.BetaEnergy();	
	_numOrb=vector<int> (2,0);
	_occupation_number=vector<vector<double>> (2);
	_occupation_number[0]=log.AlphaOccupation();
	_occupation_number[1]=log.BetaOccupation();
	_descriptors=Descriptors(log, Table);

	vector<int> Ltypes = log.Ltypes();
	int nShells = Ltypes.size();

	int lmax=0;
	for(int i=0; i<nShells; i++)
		if(lmax<abs(Ltypes[i]))
			lmax=abs(Ltypes[i]);

	int llmax = (lmax+1)*(lmax+2)/2;
	vector<int> nPrimitivesByShell = log.NumberOfGtf();
	vector<int> nCoefs (llmax);
	
	vector<double> FactorCoefs = log.FactorCoefficients();
	vector<double> CgtfCoefs = log.CgtfCoefficients();
	vector<double> CgtfSpCoefs = log.CgtfSpCoefficients();
	vector<vector<double>> coordinatesForShells;
	vector<int> NatBasis = log.NatBasis();

	for(int i=0; i<_number_of_atoms; i++)
		for(int j=0; j<NatBasis[i]; j++)
			coordinatesForShells.push_back(log.Coordinates()[i]);

	vector<double> primitiveExponents = log.Exposants();
	vector<vector<double>> coefs (llmax, vector<double> (llmax));
	vector<vector<vector<int>>> l (3, vector<vector<int>> (llmax, vector<int> (llmax)));

	int NOrb = _numberOfMo;
	_vcgtf = vector<CGTF> (NOrb);
	string format;

	int kOrb = 0;
	int kPrimitive = 0;
	for(int nS = 0;nS<nShells; nS++)
	{
		int nM = 0;

		if(Ltypes[nS]<-1)
		{
			nM = 2*abs(Ltypes[nS])+1; /* Sperical D, F, G, ...*/
			format="Sphe";
		}
		else if(Ltypes[nS]==-1)
		{
			nM = 1; /* This a SP. Make S before */
			format="Cart";
		}
		else
		{
			nM = (Ltypes[nS]+1)*(Ltypes[nS]+2)/2;
			format="Cart";
		}

		if(Ltypes[nS]==-1)
			getlTable(0, nCoefs, coefs, l, _bino); /* This a SP. Make S before */
		else
			getlTable(Ltypes[nS], nCoefs, coefs, l, _bino); 

		for(int m=0;m<nM;m++)
		{
			int ip,j,n;
			_vcgtf[kOrb]= CGTF ();
  			j = -1;
	 		for(ip=0;ip<nPrimitivesByShell[nS];ip++)
 				for(n=0;n<nCoefs[m];n++)
	 			{
	 		   		j++;
	   				vector<double> coord_ (3);
	   				vector<int> l_ (3);
	   				for(int c=0;c<3;c++)
	   				{
	   					coord_[c] = coordinatesForShells[kPrimitive+ip][c];
						l_[c] = l[c][m][n];
	   				}
	   				GTF gtf (primitiveExponents[kPrimitive+ip], 1.0, coord_, l_, _bino);
	   				_vcgtf[kOrb].push_back(gtf);
	 				_vcgtf[kOrb].setCoef(FactorCoefs[kPrimitive+ip]*CgtfCoefs[kPrimitive+ip]*coefs[m][n]);
	 				_vcgtf[kOrb].setNumCenter(log.NumCenter()[kPrimitive+ip]);
					_vcgtf[kOrb].setLtype(getLType(l_));
					_vcgtf[kOrb].setFormat(format);
					_vcgtf[kOrb].setFactorCoef(FactorCoefs[kPrimitive+ip]);
	 			}
			kOrb++;
		}
		if(Ltypes[nS]==-1) /* This a SP. Now make P*/
		{
			getlTable(-1, nCoefs, coefs, l, _bino);
			nM = 3;
			for(int m=0;m<nM;m++)
			{
				int ip,j,n;

				_vcgtf[kOrb]= CGTF ();
          			j = -1;
	 			for(ip=0;ip<nPrimitivesByShell[nS];ip++)
 					for(n=0;n<nCoefs[m];n++)
	 				{
	 		   			j++;
	   					vector<double> coord_ (3);
	   					vector<int> l_ (3);
	   					for(int c=0;c<3;c++)
	   					{
	   						coord_[c] = coordinatesForShells[nS][c];
	   						l_[c] = l[c][m][n];
	   					}
	   				GTF gtf (primitiveExponents[kPrimitive+ip], 1.0, coord_, l_, _bino);
	   				_vcgtf[kOrb].push_back(gtf);
	   				_vcgtf[kOrb].setCoef(FactorCoefs[kOrb]*CgtfSpCoefs[kPrimitive+ip]*coefs[m][n]);
	   				_vcgtf[kOrb].setNumCenter(log.NumCenter()[kPrimitive+ip]);
					_vcgtf[kOrb].setLtype(getLType(l_));
					_vcgtf[kOrb].setFormat(format);
					_vcgtf[kOrb].setFactorCoef(FactorCoefs[kPrimitive+ip]);
	 				}
				kOrb++;
			}
		}
		kPrimitive += nPrimitivesByShell[nS];
	}

	_numberOfAo=_vcgtf.size();

	if(_numberOfAo != _numberOfMo)
	{
		cout<<"Error : Their is "<<_vcgtf.size()<<" CGTF for "<<_numberOfMo<<" basis in file."<<endl;
		cout<<"Please check your file."<<endl;
		exit(1);
	}

	if(log.NumberOfBasisFunctions() < _numberOfMo)
	{
		for(int i=0; i<_numberOfMo-log.NumberOfBasisFunctions(); i++)
		{
			vector<double> v (log.NumberOfBasisFunctions(),0);
			_coefficients[0].push_back(v);
			_coefficients[1].push_back(v);
			_occupation_number[0].push_back(0);
			_occupation_number[1].push_back(0);
		}
	}

	_vcgtf_non_normalise=_vcgtf;

	for(size_t i=0; i<_vcgtf.size(); i++)
		for(int j=0; j<_vcgtf[i].numberOfFunctions(); j++)
			_primitive_centers.push_back(_vcgtf[i].NumCenter());

	NormaliseAllBasis();

	for(size_t i=0; i<_vcgtf.size(); i++)
		_number_of_gtf+=_vcgtf[i].numberOfFunctions();
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

void Orbitals::DenormaliseAllBasis()
{
	int k;

	for(k=0;k<_numberOfAo;k++)
		_vcgtf[k].denormaliseCGTF();
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
Structure Orbitals::get_struct()
{
	return _struct;
}
Grid Orbitals::makeGrid(const Domain& d)
{
	Grid g;
	g.set_str(_struct);
	g.set_dom(d);
	g.reset();
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
	for(int i=0;i<d.N1();i++)
	{
		for(int j=0;j<d.N2();j++)
		{
			for(int k=0;k<d.N3();k++)
			{
				double rho=density(d.x(i,j,k), d.y(i,j,k), d.z(i,j,k));
				g.set_Vijkl(rho,i,j,k,0);
			}
		}
	}
	return g;
}
double Orbitals::density(double x, double y, double z)
{
	double rho=0.0;
	int n;

	if(AlphaAndBeta())
		n=1;
	else
		n=2;

	vector<double> v(_vcgtf.size());
        for(size_t k=0; k<_vcgtf.size(); k++)
		v[k] = _vcgtf[k].func(x,y,z);

	for(int j=0; j<NumberOfMo(); j++)
	for(int i=0; i<n; i++)
	{
		if(OccupationNumber()[i][j]>1e-10)
		{
			double phi=0;
        		for(size_t k=0; k<_vcgtf.size(); k++)
				phi += _coefficients[i][j][k]*v[k];
			rho+=OccupationNumber()[i][j] * phi*phi;
		}
	}
	return rho;
}
Grid Orbitals::makeOrbGrid(const Domain& d, const vector<int>& nums, const vector<int>& typesSpin)
{
	Grid g;
	g.set_str(_struct);
	g.set_dom(d);
	g.reset();
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
	for(int i=0;i<d.N1();i++)
	{
		for(int j=0;j<d.N2();j++)
		{
			for(int k=0;k<d.N3();k++)
			{
				vector<double> phy=phis(d.x(i,j,k), d.y(i,j,k), d.z(i,j,k), nums, typesSpin);
				for(int l=0; l<d.Nval();l++)
				{
					g.set_Vijkl(phy[l],i,j,k,l);
				}
			}
		}
	}
	return g;
}
vector<double> Orbitals::phis(double x, double y, double z, const vector<int>& nums, const vector<int>& typesSpin)
{
	vector<double> v(_vcgtf.size());
        for(size_t k=0; k<_vcgtf.size(); k++)
	{
		v[k] = _vcgtf[k].func(x,y,z);
	}
	vector<double> values(nums.size(),0);
	for(size_t jj=0; jj<nums.size(); jj++)
	{
		int j=nums[jj];
		int i=typesSpin[jj];
		values[jj] = 0;
        	for(size_t k=0; k<_vcgtf.size(); k++)
		{
			values[jj] += _coefficients[i][j][k]*v[k];
		}
	}
	return values;
}
//epsilon=0 for Becke, epsilon =2.87e-5 for Savin. see Can. J. Chem. Vol. 74,1996 page 1088.
double Orbitals::ELF(const double& x, const double& y, const double& z, double epsilon)
{
	double rho=0.0;
	double sphi=0.0;
	double cf = 3.0/10.0*pow(3*M_PI*M_PI,2.0/3);
	int n;

	if(AlphaAndBeta())
		n=1;
	else
		n=2;

	vector<double> v(_vcgtf.size());
        for(size_t k=0; k<_vcgtf.size(); k++)
		v[k] = _vcgtf[k].func(x,y,z);

	vector<double> A(_vcgtf.size());
	vector< vector<double> > vg(3,A);
	for(size_t k=0;k<_vcgtf.size();k++)
	{
		for(int c=0;c<3;c++)
			vg[c][k]=_vcgtf[k].grad_CGTF(x,y,z,c);
	}

	double v1[3]={0,0,0};
	for(int j=0; j<NumberOfMo(); j++)
	for(int i=0; i<n; i++)
	{
		if(OccupationNumber()[i][j]>1e-10)
		{
			double phi=0;
        		for(size_t k=0; k<_vcgtf.size(); k++)
				phi += _coefficients[i][j][k]*v[k];
			rho+=OccupationNumber()[i][j] * phi*phi;

			for(int c=0;c<3;c++)
			{
				double dp=0;
				for(size_t k=0; k<_vcgtf.size(); k++)
					dp += _coefficients[i][j][k]*vg[c][k];
				sphi += OccupationNumber()[i][j] * dp*dp;
				v1[c] += OccupationNumber()[i][j] * phi*dp;
			}
		}
	}
	double  grho2=0;
        for(int c=0;c<3;c++)
		grho2 += v1[c]*v1[c]*4;
	double t = sphi/2 - grho2/8.0/rho;
	double th = cf*pow(rho,5.0/3.0);
	double XS2 = (t+epsilon)/th;
	XS2 = XS2*XS2;
	return 1.0/(1.0+XS2);
}

Grid Orbitals::makeELFgrid(const Domain& d,const double& epsilon)
{
	Grid g;
	g.set_str(_struct);
	g.set_dom(d);
	g.reset();
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
	for(int i=0;i<d.N1();i++)
	{
		for(int j=0;j<d.N2();j++)
		{
			for(int k=0;k<d.N3();k++)
			{
				double elf=ELF(d.x(i,j,k), d.y(i,j,k), d.z(i,j,k), epsilon);
				g.set_Vijkl(elf,i,j,k,0);
			}
		}
	}
	return g;
}
