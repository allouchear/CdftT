#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <Basis/GTF.h>
#include <Basis/CGTF.h>
#include <Common/PeriodicTable.h>
#include <Orbitals/Orbitals.h>
#include <Utils/Enums.hpp>
#include <Utils/FCHK.h>
#include <Utils/LM.h>
#include <Utils/LOG.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/Utils.h>
#include <Utils/WFX.h>


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Orbitals::Orbitals():
    _vcgtf(),
    _vcgtfUnnormalized(),
    _coefficients(),
    _numberOfAo(0),
    _numberOfGtf(0),
    _numberOfMo(0),
    _numberOfAlphaElectrons(0),
    _numberOfBetaElectrons(0),
    _numberOfAtoms(0),
    _primitiveCenters(),
    _struct(),
    _atomicNumbers(),
    _symbol(),
    _orbitalEnergy(),
    _all_f(),
    _homoLumoIndexes({ -1, -1 }),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(),
    _descriptors(),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{ }

Orbitals::Orbitals(WFX& wfxParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _vcgtfUnnormalized(),
    _coefficients(2),
    _numberOfAo(0),
    _numberOfGtf(0),
    _numberOfMo(0),
    _numberOfAlphaElectrons(0),
    _numberOfBetaElectrons(0),
    _numberOfAtoms(0),
    _primitiveCenters(),
    _struct(wfxParser, periodicTable),
    _atomicNumbers(),
    _symbol(),
    _orbitalEnergy(),
    _all_f(),
    _homoLumoIndexes({ -1, -1 }),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{
    _coordinates = wfxParser.Nuclear_Cartesian_Coordinates();

    _vcgtf = std::vector<CGTF> (wfxParser.Number_of_Primitives());

    std::vector<std::vector<double>> coord(wfxParser.Number_of_Nuclei(), std::vector<double>(0));
    
    for(int i = 0; i < wfxParser.Number_of_Nuclei(); ++i)
    {
        for(int j = i * 3; j < 3 * (1 + i); ++j)
        {
            coord[i].push_back(wfxParser.Nuclear_Cartesian_Coordinates()[j]);
        }
    }

    GTF gtf;
    for(int j = 0; j < wfxParser.Number_of_Primitives(); ++j)
    {
        std::array<double, 3> coords ({ coord[wfxParser.Primitive_Centers()[j] - 1][0],
                                        coord[wfxParser.Primitive_Centers()[j] - 1][1],
                                        coord[wfxParser.Primitive_Centers()[j] - 1][2] });
        
        gtf.push_back(wfxParser.Primitive_Exponents()[j], 1.0, coords, setLxyz(wfxParser.Primitive_Types()[j]), bino);
        _vcgtf[j].push_back(gtf);
        _vcgtf[j].setCoef(1.0);
        _vcgtf[j].setFactorCoef(1.0);
        _vcgtf[j].setNumCenter(wfxParser.Primitive_Centers()[j]);
        _vcgtf[j].setLtype(getLType(_vcgtf[j].gtf()[0].get_l()));
        _vcgtf[j].setFormat("Cart");
    }

    //_coefficients=std::vector<std::vector<std::vector<double>>> (2);
    _coefficients[0] = std::vector<std::vector<double>>(wfxParser.Molecular_Orbital_Primitive_Coefficients()[0].size());
    _coefficients[1] = std::vector<std::vector<double>>(wfxParser.Molecular_Orbital_Primitive_Coefficients()[1].size());

    for(size_t i = 0; i < wfxParser.Molecular_Orbital_Primitive_Coefficients()[0].size(); ++i)
    {
        _coefficients[0][i] = wfxParser.Molecular_Orbital_Primitive_Coefficients()[0][i].Coefficients();
    }
    for(size_t i = 0; i < wfxParser.Molecular_Orbital_Primitive_Coefficients()[1].size(); ++i)
    {
        _coefficients[1][i] = wfxParser.Molecular_Orbital_Primitive_Coefficients()[1][i].Coefficients();
    }

    _primitiveCenters = wfxParser.Primitive_Centers();
    _atomicNumbers = wfxParser.Atomic_Number();
    _numberOfMo = wfxParser.Number_of_Occupied_Molecular_Orbital();

    _numberOfGtf = wfxParser.Number_of_Primitives();
    if(!wfxParser.AlphaAndBeta())
    {
        _numberOfMo /= 2;
    }

    _numberOfAlphaElectrons = wfxParser.Number_of_Alpha_Electrons();
    _numberOfBetaElectrons = wfxParser.Number_of_Beta_Electrons();

    _numberOfAtoms = wfxParser.Number_of_Nuclei();
    _orbitalEnergy = wfxParser.Molecular_Orbital_Energies();
    _symbol = wfxParser.Nuclear_Names();
    _energy = wfxParser.Energy();
    _occupationNumber = wfxParser.Molecular_Orbital_Occupation_Numbers();
    _alphaAndBeta = wfxParser.AlphaAndBeta();
    _descriptors = Descriptors(wfxParser, periodicTable);
    _numberOfAo = _vcgtf.size();
    _vcgtfUnnormalized = _vcgtf;
}

Orbitals::Orbitals(FCHK& fchkParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _vcgtfUnnormalized(),
    _coefficients(2),
    _numberOfAo(0),
    _numberOfGtf(0),
    _numberOfMo(0),
    _numberOfAlphaElectrons(0),
    _numberOfBetaElectrons(0),
    _numberOfAtoms(0),
    _primitiveCenters(),
    _struct(fchkParser, periodicTable),
    _atomicNumbers(),
    _symbol(),
    _orbitalEnergy(),
    _all_f(),
    _homoLumoIndexes({ -1, -1 }),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{
    _numberOfMo = fchkParser.NumberOfBasisFunctions();
    _coordinates = fchkParser.CurrentCartesianCoordinates();
    _energy = fchkParser.ScfEnergy();
    _mixte = fchkParser.Mixte();

    int lmax = fchkParser.HighestAngularMomentum();
    int nShells = fchkParser.NumberOfContractedShells();
    int llmax = (lmax + 1) * (lmax + 2) / 2;

    std::vector<int> nCoefs(llmax);
    std::vector<std::vector<double>> coefs(llmax, std::vector<double>(llmax));
    std::vector<std::vector<std::vector<int>>> l(3, std::vector<std::vector<int>>(llmax, std::vector<int>(llmax)));

    std::vector<int> numAtoms = fchkParser.ShellToAtomMap();
    std::vector<int> nPrimitivesByShell = fchkParser.NumberOfPrimitivesPerShell();
    std::vector<int> shellTypes = fchkParser.ShellTypes();

    std::vector<double> contractionsCoefs = fchkParser.ContractionCoefficients();
    std::vector<double> contractionsCoefsSP = fchkParser.spContractionCoefficients();
    std::vector<double> coordinatesForShells = fchkParser.CoordinatesForShells();
    std::vector<double> primitiveExponents = fchkParser.PrimitiveExponents();

    int NOrb = 0;
    for(int nS = 0; nS < nShells; ++nS) 
    {
        if(shellTypes[nS]<-1)
        {
            NOrb += 2 * std::abs(shellTypes[nS]) + 1; /* Spherical D, F, G, ...*/
        }
        else if(shellTypes[nS] == -1)
        {
            NOrb += 4; /* This a SP.*/
        }
        else
        {
            NOrb += (shellTypes[nS] + 1) * (shellTypes[nS] + 2) / 2; /* Cartesian S,P,D,F,G,..*/
        }
    }

    _vcgtf = std::vector<CGTF>(NOrb);
    int kOrb = 0;
    int kPrimitive = 0;
    std::string format;

    for(int nS = 0; nS < nShells; ++nS)
    {
        int nM = 0;
        /* printf("begin primitive nS = %d\n",nS);*/
        if(shellTypes[nS] < -1)
        {
            nM = 2 * std::abs(shellTypes[nS]) + 1; /* Sperical D, F, G, ...*/
            format = "Sphe";
        }
        else if(shellTypes[nS] == -1)
        {
            nM = 1; /* This a SP. Make S before */
            format = "Cart";
        }
        else
        {
            nM = (shellTypes[nS] + 1) * (shellTypes[nS] + 2) / 2;
            format = "Cart";
        }

        /* printf("nM = %d\n",nM);*/
        if(shellTypes[nS] == -1)
        {
            getlTable(0, nCoefs, coefs, l, _bino); /* This a SP. Make S before */
        }
        else
        {
            getlTable(shellTypes[nS], nCoefs, coefs, l, _bino); 
        }

        /* printf("end getlTable\n");*/
        for(int m = 0; m < nM; ++m)
        {
            int ip,j,n;
            /* printf("P : m = %d nCoef = %d nPrim = %d\n",m,nCoefs[m],nPrimitivesByShell[nS]);*/

            _vcgtf[kOrb]= CGTF ();
            _vcgtf[kOrb].setNumCenter(numAtoms[nS]);
            _vcgtf[kOrb].setFactorCoef(1.0);
            _vcgtf[kOrb].setFormat(format);

            j = -1;
            for(ip = 0; ip < nPrimitivesByShell[nS]; ++ip)
            {
                for(n = 0; n < nCoefs[m]; ++n)
                {
                    j++;

                    std::array<double, 3> coord_;
                    std::vector<int> l_ (3);
                    for(int c = 0; c < 3; ++c)
                    {
                        coord_[c] = coordinatesForShells[c + nS * 3];
                        l_[c] = l[c][m][n];
                    }

                    GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
                    _vcgtf[kOrb].push_back(gtf);
                    _vcgtf[kOrb].setCoef(contractionsCoefs[kPrimitive+ip]*coefs[m][n]);
                    _vcgtf[kOrb].setLtype(getLType(l_));
                }
            }

            kOrb++;
        }

        if(shellTypes[nS] == -1) /* This a SP. Now make P*/
        {
            getlTable(-1, nCoefs, coefs, l, _bino);

            nM = 3;
            for(int m = 0; m < nM; ++m)
            {
                int ip,j,n;
                /* printf("P : m = %d nCoef = %d nPrim = %d\n",m,nCoefs[m],nPrimitivesByShell[nS]);*/

                _vcgtf[kOrb]= CGTF ();
                _vcgtf[kOrb].setNumCenter(numAtoms[nS]);
                _vcgtf[kOrb].setFactorCoef(1.0);
                _vcgtf[kOrb].setFormat(format);

                j = -1;
                for(ip = 0; ip < nPrimitivesByShell[nS]; ++ip)
                {
                    for(n = 0; n < nCoefs[m]; ++n)
                    {
                        j++;

                        std::array<double, 3> coord_ ;
                        std::vector<int> l_ (3);
                        for(int c = 0; c < 3; ++c)
                        {
                            coord_[c] = coordinatesForShells[c + nS * 3];
                            l_[c] = l[c][m][n];
                        }

                        GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
                        _vcgtf[kOrb].push_back(gtf);
                        _vcgtf[kOrb].setCoef(contractionsCoefsSP[kPrimitive+ip]*coefs[m][n]);
                        _vcgtf[kOrb].setLtype(getLType(l_));
                    }
                }

                kOrb++;
            }
        }
        /* printf("end primitive nS = %d\n",nS);*/
        kPrimitive += nPrimitivesByShell[nS];
    }

    _numberOfAo = _vcgtf.size();
    if(_numberOfAo != _numberOfMo)
    {
        std::stringstream errorMessage;
        errorMessage << "Error : There are " << _vcgtf.size() << " CGTFs for " << _numberOfMo << " basis in file." << std::endl;
        errorMessage << "Please check your file.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    _numberOfAlphaElectrons = fchkParser.NumberOfAlphaElectrons();
    _numberOfBetaElectrons = fchkParser.NumberOfBetaElectrons();
    _numberOfAtoms = fchkParser.NumberOfAtoms();
    _coefficients = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>());

    int nOrb_alpha = fchkParser.AlphaMOCoefficients().size() / fchkParser.AlphaOrbitalEnergies().size();
    int nOrb_beta = fchkParser.BetaMOCoefficients().size() / fchkParser.BetaOrbitalEnergies().size();

    _coefficients[0] = std::vector<std::vector<double>>(nOrb_alpha, std::vector<double>(_numberOfMo));
    _coefficients[1] = std::vector<std::vector<double>>(nOrb_beta, std::vector<double>(_numberOfMo));

    for(int i = 0; i < nOrb_alpha; ++i)
    {
        for(int j = 0; j < _numberOfMo; ++j)
        {
            _coefficients[0][i][j] = fchkParser.AlphaMOCoefficients()[i * _numberOfMo + j];
        }
    }

    for(int i = 0; i < nOrb_beta; ++i)
    {
        for(int j = 0; j < _numberOfMo; ++j)
        {
            _coefficients[1][i][j] = fchkParser.BetaMOCoefficients()[i * _numberOfMo + j];
        }
    }

    _orbitalEnergy = std::vector<std::vector<double>>(2, std::vector<double>());
    _orbitalEnergy[0] = fchkParser.AlphaOrbitalEnergies();
    _orbitalEnergy[1] = fchkParser.BetaOrbitalEnergies();

    _occupationNumber = std::vector<std::vector<double>>(2);

    _alphaAndBeta = fchkParser.AlphaAndBeta();
    if(_alphaAndBeta)
    {
        _occupationNumber[0] = std::vector<double>(_numberOfMo);
        for(int i = 0; i < _numberOfMo; ++i)
        {
            _occupationNumber[0][i] = fchkParser.AlphaOccupation()[i] + fchkParser.BetaOccupation()[i];
        }
        _occupationNumber[1] = _occupationNumber[0];
    }
    else
    {
        _occupationNumber[0] = fchkParser.AlphaOccupation();
        _occupationNumber[1] = fchkParser.BetaOccupation();
    }

    _atomicNumbers = fchkParser.AtomicNumbers();
    _symbol = std::vector<std::string>(_numberOfAtoms);
    
    for(int i = 0; i < _numberOfAtoms; ++i)
    {
        _symbol[i] = periodicTable.element(_atomicNumbers[i]).get_symbol();
    }

    _descriptors = Descriptors(fchkParser, periodicTable);

    for(size_t i = 0; i < _vcgtf.size(); i++)
    {
        _primitiveCenters.push_back(_vcgtf[i].NumCenter());
    }

    _vcgtfUnnormalized = _vcgtf;
    NormaliseAllBasis();

    for(size_t i = 0; i < _vcgtf.size(); i++)
    {
        _numberOfGtf += _vcgtf[i].numberOfFunctions();
    }
}

Orbitals::Orbitals(MOLDENGAB& moldengabParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _vcgtfUnnormalized(),
    _coefficients(2),
    _numberOfAo(0),
    _numberOfGtf(0),
    _numberOfMo(0),
    _numberOfAlphaElectrons(0),
    _numberOfBetaElectrons(0),
    _numberOfAtoms(0),
    _primitiveCenters(),
    _struct(moldengabParser, periodicTable),
    _atomicNumbers(),
    _symbol(),
    _orbitalEnergy(),
    _all_f(),
    _homoLumoIndexes({ -1, -1 }),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{
    _coefficients[0] = moldengabParser.AlphaMOCoefs();
    _coefficients[1] = moldengabParser.BetaMOCoefs();

    _numberOfAtoms = moldengabParser.NumberOfAtoms();

    _coordinates = std::vector<double>(_numberOfAtoms * 3);
    for (int i = 0; i < _numberOfAtoms; ++i)
    {
        int k = 0;
        for(int j = i * 3; j < (i + 1) * 3; ++j)
        {
            _coordinates[j] = moldengabParser.Coordinates()[i][k];
            k++;
        }
    }

    _numberOfMo = moldengabParser.NumberOfMO();
    _alphaAndBeta = moldengabParser.AlphaAndBeta();

    if(!_alphaAndBeta)
    {
        _numberOfMo /= 2;
    }

    if(_alphaAndBeta)
    {
        for(size_t i = 0; i < moldengabParser.AlphaOccupation().size(); i++)
        {
            if(moldengabParser.AlphaOccupation()[i] == 1)
            {
                _numberOfAlphaElectrons++;
            }
        }

        _numberOfBetaElectrons = _numberOfAlphaElectrons;
    }
    else
    {
        for(size_t i = 0; i < moldengabParser.AlphaOccupation().size(); i++)
        {
            if(moldengabParser.AlphaOccupation()[i] == 1)
            {
                _numberOfAlphaElectrons++;
            }
        }

        for(size_t i = 0; i < moldengabParser.BetaOccupation().size(); i++)
        {
            if(moldengabParser.BetaOccupation()[i] == 1)
            {
                _numberOfBetaElectrons++;
            }
        }
    }

    _atomicNumbers = moldengabParser.AtomicNumbers();
    _symbol = moldengabParser.Symbol();

    _orbitalEnergy = std::vector<std::vector<double>>(2);
    _orbitalEnergy[0] = moldengabParser.AlphaEnergies();
    _orbitalEnergy[1] = moldengabParser.BetaEnergies();

    _occupationNumber = std::vector<std::vector<double>>(2);
    _occupationNumber[0] = moldengabParser.AlphaOccupation();
    _occupationNumber[1] = moldengabParser.BetaOccupation();

    _descriptors = Descriptors(moldengabParser, periodicTable);

    std::vector<int> Ltypes = moldengabParser.Ltypes();

    int nShells = Ltypes.size();
    int lmax=0;
    for(int i = 0; i < nShells; ++i)
    {
        if(lmax < std::abs(Ltypes[i]))
        {
            lmax = std::abs(Ltypes[i]);
        }
    }

    int llmax = (lmax + 1) * (lmax + 2) / 2;
    std::vector<int> nPrimitivesByShell = moldengabParser.NumberOfGtf();
    std::vector<int> nCoefs (llmax);
    std::vector<std::vector<double>> coefs(llmax, std::vector<double>(llmax));
    std::vector<std::vector<std::vector<int>>> l(3, std::vector<std::vector<int>>(llmax, std::vector<int>(llmax)));

    std::vector<double> FactorCoefs = moldengabParser.FactorCoefficients();
    std::vector<double> CgtfCoefs = moldengabParser.CgtfCoefficients();
    std::vector<std::array<double, 3>> coordinatesForShells;
    std::vector<int> NatBasis = moldengabParser.NatBasis();

    for(int i = 0; i < _numberOfAtoms; ++i)
    {
        for(int j = 0; j < NatBasis[i]; ++j)
        {
            coordinatesForShells.push_back(moldengabParser.Coordinates()[i]);
        }
    }

    std::vector<double> primitiveExponents = moldengabParser.Exposants();

    int NOrb = moldengabParser.NumberOfMOCoefs();

    _vcgtf = std::vector<CGTF>(NOrb);

    std::string format;

    int kOrb = 0;
    int kPrimitive = 0;
    for(int nS = 0; nS < nShells; ++nS)
    {
        int nM = 0;

        if(Ltypes[nS] < -1)
        {
            nM = 2 * std::abs(Ltypes[nS]) + 1; /* Sperical D, F, G, ...*/
            format = "Sphe";
        }
        else if(Ltypes[nS] == -1)
        {
            nM = 1; /* This a SP. Make S before */
            format = "Cart";
        }
        else
        {
            nM = (Ltypes[nS] + 1) * (Ltypes[nS] + 2) / 2;
            format = "Cart";
        }

        if(Ltypes[nS] == -1)
        {
            getlTable(0, nCoefs, coefs, l, _bino); /* This a SP. Make S before */
        }
        else
        {
            getlTable(Ltypes[nS], nCoefs, coefs, l, _bino);
        }

        for(int m = 0; m < nM; ++m)
        {
            int ip,j,n;

            _vcgtf[kOrb]= CGTF ();
            _vcgtf[kOrb].setNumCenter(moldengabParser.Numcenter()[nS]);
            _vcgtf[kOrb].setLtype(moldengabParser.ShellTypes()[nS]);
            _vcgtf[kOrb].setFactorCoef(FactorCoefs[nS]);
            _vcgtf[kOrb].setFormat(format);

            j = -1;
            for(ip = 0; ip < nPrimitivesByShell[nS]; ++ip)
            {
                for(n = 0; n < nCoefs[m]; ++n)
                {
                    j++;

                    std::array<double, 3> coord_;
                    std::vector<int> l_(3);
                    for(int c = 0; c < 3; ++c)
                    {
                        coord_[c] = coordinatesForShells[kPrimitive + ip][c];
                        l_[c] = l[c][m][n];
                    }

                    GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
                    _vcgtf[kOrb].push_back(gtf);
                    _vcgtf[kOrb].setCoef(FactorCoefs[kPrimitive+ip]*CgtfCoefs[kPrimitive+ip]*coefs[m][n]);
                }
            }

            kOrb++;
        }

        if(Ltypes[nS] == -1) /* This a SP. Now make P*/
        {
            getlTable(-1, nCoefs, coefs, l, _bino);
            nM = 3;
            for(int m = 0; m < nM; ++m)
            {
                int ip,j,n;

                _vcgtf[kOrb]= CGTF ();
                _vcgtf[kOrb].setNumCenter(moldengabParser.Numcenter()[nS]);
                _vcgtf[kOrb].setLtype(moldengabParser.ShellTypes()[nS]);
                _vcgtf[kOrb].setFactorCoef(FactorCoefs[nS]);
                _vcgtf[kOrb].setFormat(format);

                j = -1;
                for(ip=0;ip<nPrimitivesByShell[nS];ip++)
                {
                    for(n=0;n<nCoefs[m];n++)
                    {
                        j++;
                        std::array<double, 3> coord_;
                        std::vector<int> l_(3);
                        for(int c = 0; c < 3; ++c)
                        {
                            coord_[c] = coordinatesForShells[kPrimitive + ip][c];
                            l_[c] = l[c][m][n];
                        }

                        GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
                        _vcgtf[kOrb].push_back(gtf);
                        _vcgtf[kOrb].setCoef(FactorCoefs[kPrimitive + ip] * CgtfCoefs[kPrimitive + ip] * coefs[m][n]);
                    }
                }

                kOrb++;
            }
        }

        kPrimitive += nPrimitivesByShell[nS];
    }

    _numberOfAo=_vcgtf.size();

    if(_numberOfAo != moldengabParser.NumberOfMOCoefs())
    {
        std::stringstream errorMessage;
        errorMessage << "Error : There are " << _vcgtf.size() << " CGTFs for " << moldengabParser.NumberOfMOCoefs() << " basis in file." << std::endl;
        errorMessage << "Please check your file.";

        print_error(errorMessage.str());

        std::exit(1);
    }

    _vcgtfUnnormalized = _vcgtf;

    for(size_t i = 0; i < _vcgtf.size(); ++i)
    {
        _primitiveCenters.push_back(_vcgtf[i].NumCenter());
    }

    NormaliseAllBasis();

    for(size_t i = 0; i < _vcgtf.size(); ++i)
    {
        _numberOfGtf += _vcgtf[i].numberOfFunctions();
    }

    _mixte = moldengabParser.Mixte();

    //Sorting(); Don't work
}

Orbitals::Orbitals(LOG& logParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _vcgtfUnnormalized(),
    _coefficients(2),
    _numberOfAo(0),
    _numberOfGtf(0),
    _numberOfMo(0),
    _numberOfAlphaElectrons(0),
    _numberOfBetaElectrons(0),
    _numberOfAtoms(0),
    _primitiveCenters(),
    _struct(logParser, periodicTable),
    _atomicNumbers(),
    _symbol(),
    _orbitalEnergy(),
    _all_f(),
    _homoLumoIndexes({ -1, -1 }),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{
    _struct = Structure(logParser, periodicTable);
    _vcgtfUnnormalized = std::vector<CGTF> ();
    _numberOfGtf = 0;
    _coefficients = std::vector<std::vector<std::vector<double>>> (2);
    _coefficients[0] = logParser.AlphaMOcoefs();
    _coefficients[1] = logParser.BetaMOcoefs();
    _energy = logParser.Energy();
    _numberOfAlphaElectrons = logParser.NumberOfAlphaElectrons();
    _numberOfBetaElectrons = logParser.NumberOfBetaElectrons();
    _numberOfAtoms = logParser.NumberOfAtoms();
    _coordinates = std::vector<double> (_numberOfAtoms*3);
    for(int i = 0; i < _numberOfAtoms; ++i)
    {
        int k = 0;
        for(int j = i * 3; j < (i +1) * 3; ++j)
        {
            _coordinates[j] = logParser.Coordinates()[i][k];
            k++;
        }
    }

    _numberOfMo = logParser.NumberOfMO();
    _alphaAndBeta = logParser.AlphaAndBeta();

    _bino = bino;
        
    _atomicNumbers = logParser.AtomicNumbers();
    _symbol = logParser.Symbol();
    _orbitalEnergy = std::vector<std::vector<double>> (2);
    _orbitalEnergy[0] = logParser.AlphaEnergy();
    _orbitalEnergy[1] = logParser.BetaEnergy();
    _occupationNumber = std::vector<std::vector<double>> (2);
    _occupationNumber[0] = logParser.AlphaOccupation();
    _occupationNumber[1] = logParser.BetaOccupation();
    _descriptors = Descriptors(logParser, periodicTable);

    std::vector<int> Ltypes = logParser.Ltypes();
    int nShells = Ltypes.size();

    int lmax = 0;
    for(int i = 0; i < nShells; ++i)
    {
        if(lmax < std::abs(Ltypes[i]))
        {
            lmax = std::abs(Ltypes[i]);
        }
    }

    int llmax = (lmax + 1) * (lmax + 2) / 2;
    std::vector<int> nPrimitivesByShell = logParser.NumberOfGtf();
    std::vector<int> nCoefs(llmax);
    
    std::vector<double> FactorCoefs = logParser.FactorCoefficients();
    std::vector<double> CgtfCoefs = logParser.CgtfCoefficients();
    std::vector<double> CgtfSpCoefs = logParser.CgtfSpCoefficients();
    std::vector<std::vector<double>> coordinatesForShells;
    std::vector<int> NatBasis = logParser.NatBasis();

    for(int i = 0; i < _numberOfAtoms; ++i)
    {
        for(int j = 0; j < NatBasis[i]; ++j)
        {
            coordinatesForShells.push_back(logParser.Coordinates()[i]);
        }
    }

    std::vector<double> primitiveExponents = logParser.Exposants();
    std::vector<std::vector<double>> coefs(llmax, std::vector<double> (llmax));
    std::vector<std::vector<std::vector<int>>> l(3, std::vector<std::vector<int>>(llmax, std::vector<int>(llmax)));

    int NOrb = _numberOfMo;
    _vcgtf = std::vector<CGTF> (NOrb);
    std::string format;

    int kOrb = 0;
    int kPrimitive = 0;
    for(int nS = 0; nS < nShells; ++nS)
    {
        int nM = 0;

        if(Ltypes[nS] < -1)
        {
            nM = 2 * std::abs(Ltypes[nS]) + 1; /* Spherical D, F, G, ...*/
            format = "Sphe";
        }
        else if(Ltypes[nS] == -1)
        {
            nM = 1; /* This a SP. Make S before */
            format = "Cart";
        }
        else
        {
            nM = (Ltypes[nS] + 1) * (Ltypes[nS] + 2) / 2;
            format = "Cart";
        }

        if(Ltypes[nS] == -1) // This a SP. Make S before
        {
            getlTable(0, nCoefs, coefs, l, _bino);
        } 
        else
        {
            getlTable(Ltypes[nS], nCoefs, coefs, l, _bino);
        } 

        for(int m = 0; m < nM; ++m)
        {
            int ip,j,n;
            _vcgtf[kOrb]= CGTF ();
            j = -1;
            
            for(ip = 0; ip < nPrimitivesByShell[nS]; ++ip)
            {
                for(n = 0; n < nCoefs[m]; ++n)
                {
                    j++;

                    std::array<double, 3> coord_;
                    std::vector<int> l_(3);
                    for(int c = 0; c < 3; ++c)
                    {
                        coord_[c] = coordinatesForShells[kPrimitive + ip][c];
                        l_[c] = l[c][m][n];
                    }

                    GTF gtf(primitiveExponents[kPrimitive + ip], 1.0, coord_, l_, _bino);
                    _vcgtf[kOrb].push_back(gtf);
                    _vcgtf[kOrb].setCoef(FactorCoefs[kPrimitive + ip] * CgtfCoefs[kPrimitive + ip] * coefs[m][n]);
                    _vcgtf[kOrb].setNumCenter(logParser.NumCenter()[kPrimitive + ip]);
                    _vcgtf[kOrb].setLtype(getLType(l_));
                    _vcgtf[kOrb].setFormat(format);
                    _vcgtf[kOrb].setFactorCoef(FactorCoefs[kPrimitive + ip]);
                }
            }

            kOrb++;
        }

        if(Ltypes[nS] == -1) /* This a SP. Now make P*/
        {
            getlTable(-1, nCoefs, coefs, l, _bino);
            nM = 3;

            for(int m = 0; m < nM; ++m)
            {
                int ip,j,n;

                _vcgtf[kOrb]= CGTF ();
                j = -1;
                for(ip = 0; ip < nPrimitivesByShell[nS]; ++ip)
                {
                    for(n = 0; n < nCoefs[m]; ++n)
                    {
                        j++;

                        std::array<double, 3> coord_;
                        std::vector<int> l_(3);
                        for(int c = 0; c < 3; ++c)
                        {
                            //coord_[c] = coordinatesForShells[nS][c];
                            coord_[c] = coordinatesForShells[kPrimitive + ip][c];
                            l_[c] = l[c][m][n];
                        }

                        GTF gtf(primitiveExponents[kPrimitive + ip], 1.0, coord_, l_, _bino);
                        _vcgtf[kOrb].push_back(gtf);

                        //_vcgtf[kOrb].setCoef(FactorCoefs[kOrb] * CgtfSpCoefs[kPrimitive + ip] * coefs[m][n]);
                        _vcgtf[kOrb].setCoef(FactorCoefs[kPrimitive + ip] * CgtfSpCoefs[kPrimitive + ip] * coefs[m][n]);
                        _vcgtf[kOrb].setNumCenter(logParser.NumCenter()[kPrimitive + ip]);
                        _vcgtf[kOrb].setLtype(getLType(l_));
                        _vcgtf[kOrb].setFormat(format);
                        _vcgtf[kOrb].setFactorCoef(FactorCoefs[kPrimitive + ip]);
                    }
                }

                kOrb++;
            }
        }

        kPrimitive += nPrimitivesByShell[nS];
    }

    _numberOfAo = _vcgtf.size();

    if(_numberOfAo != _numberOfMo)
    {
        std::cout << "Error : There are " << _vcgtf.size() << " CGTFs for " << _numberOfMo << " basis in file." << std::endl;
        std::cout << "Please check your file." << std::endl;

        std::exit(1);
    }

    if (logParser.NumberOfBasisFunctions() < _numberOfMo)
    {
        for(int i = 0; i < _numberOfMo - logParser.NumberOfBasisFunctions(); ++i)
        {
            std::vector<double> v(logParser.NumberOfBasisFunctions(), 0);
            _coefficients[0].push_back(v);
            _coefficients[1].push_back(v);
            _occupationNumber[0].push_back(0);
            _occupationNumber[1].push_back(0);
        }
    }

    _vcgtfUnnormalized = _vcgtf;

    for(size_t i = 0; i < _vcgtf.size(); ++i)
    {
        _primitiveCenters.push_back(_vcgtf[i].NumCenter());
    }

    NormaliseAllBasis();

    for(size_t i = 0; i < _vcgtf.size(); ++i)
    {
        _numberOfGtf += _vcgtf[i].numberOfFunctions();
    }

    _mixte = logParser.Mixte();
}


//----------------------------------------------------------------------------------------------------//
// PRIVATE METHODS
//----------------------------------------------------------------------------------------------------//

double Orbitals::density(int orbitalIndex, SpinType spinType, const std::vector<double>& evaluatedCgtfs) const
{
    double rho = 0.0;

    // Handle ALPHA_BETA case
    std::vector<SpinType> spins;
    if (spinType == SpinType::ALPHA || spinType == SpinType::ALPHA_BETA)
    {
        spins.push_back(SpinType::ALPHA);
    }
    if (spinType == SpinType::BETA || spinType == SpinType::ALPHA_BETA)
    {
        spins.push_back(SpinType::BETA);
    }

    for (SpinType spinType : spins)
    {
        int spin = static_cast<int>(spinType);

        if (_occupationNumber[spin][orbitalIndex] > 1e-10)
        {
            double phi = 0;

            for (size_t k = 0; k < evaluatedCgtfs.size(); ++k)
            {
                phi += _coefficients[spin][orbitalIndex][k] * evaluatedCgtfs[k];
            }

            rho += _occupationNumber[spin][orbitalIndex] * phi * phi;
        }
    }

    return rho;
}

void Orbitals::evaluateCgtfsAtPoint(std::vector<double>& evaluatedCgtfs, double x, double y, double z) const
{
    if (evaluatedCgtfs.size() != _vcgtf.size())
    {
        evaluatedCgtfs.resize(_vcgtf.size());
    }

    for (size_t k = 0; k < _vcgtf.size(); ++k)
    {
        evaluatedCgtfs[k] = _vcgtf[k].func(x, y, z);
    }
}

double Orbitals::phiSquared(int orbitalIndex, SpinType spinType, const std::vector<double>& evaluatedCgtfs) const
{
    double phiSquared = 0.0;

    // Handle ALPHA_BETA case
    std::vector<SpinType> spins;
    if (spinType == SpinType::ALPHA || spinType == SpinType::ALPHA_BETA)
    {
        spins.push_back(SpinType::ALPHA);
    }
    if (spinType == SpinType::BETA || spinType == SpinType::ALPHA_BETA)
    {
        spins.push_back(SpinType::BETA);
    }

    for (SpinType spinType : spins)
    {
        int spin = static_cast<int>(spinType);

        if (_occupationNumber[spin][orbitalIndex] > 1e-10)
        {
            double phi = 0;

            for (size_t k = 0; k < evaluatedCgtfs.size(); ++k)
            {
                phi += _coefficients[spin][orbitalIndex][k] * evaluatedCgtfs[k];
            }

            phiSquared += phi * phi;
        }
    }

    return phiSquared;
}

//----------------------------------------------------------------------------------------------------//
// GETTERS
//----------------------------------------------------------------------------------------------------//

const std::vector<std::vector<double>>& Orbitals::get_all_f() const
{
    return _all_f;
}

std::vector<CGTF> Orbitals::get_vcgtf() const
{
    return _vcgtf;
}

const std::vector<std::vector<std::vector<double>>>& Orbitals::get_coefficients() const
{
    return _coefficients;
}

int Orbitals::get_numberOfAo() const
{
    return _numberOfAo;
}

int Orbitals::get_numberOfMo() const
{
    return _numberOfMo;
}

int Orbitals::get_numberOfAtoms() const
{
    return _numberOfAtoms;
}

const std::vector<int>& Orbitals::get_primitiveCenters() const
{
    return _primitiveCenters;
}

const Structure& Orbitals::get_struct() const
{
    return _struct;
}

const std::vector<std::string>& Orbitals::get_symbol() const
{
    return _symbol;
}

const std::vector<std::vector<double>>& Orbitals::get_orbitalEnergy() const
{
    return _orbitalEnergy;
}

const std::vector<std::vector<double>>& Orbitals::get_occupationNumber() const
{
    return _occupationNumber;
}

bool Orbitals::get_alphaAndBeta() const
{
    return _alphaAndBeta;
}

const Descriptors& Orbitals::get_descriptors() const
{
    return _descriptors;
}

const double Orbitals::get_energy() const
{
    return _energy;
}


//----------------------------------------------------------------------------------------------------//
// SETTERS
//----------------------------------------------------------------------------------------------------//

void Orbitals::set_coefficients(const std::vector<std::vector<std::vector<double>>>& coefficients)
{
    _coefficients = coefficients;
}

void Orbitals::set_energy(double energy)
{
    _energy = energy;
}

void Orbitals::set_homoLumoIndexes(int homoIndex, int lumoIndex)
{
    _homoLumoIndexes[0] = homoIndex;
    _homoLumoIndexes[1] = lumoIndex;
}

void Orbitals::set_numberOfAlphaElectrons(int numberOfAlphaElectrons)
{
    _numberOfAlphaElectrons = numberOfAlphaElectrons;
}

void Orbitals::set_numberOfBetaElectrons(int numberOfBetaElectrons)
{
    _numberOfBetaElectrons = numberOfBetaElectrons;
}

void Orbitals::set_numberOfMo(int numberOfMo)
{
    _numberOfMo = numberOfMo;
}

void Orbitals::set_occupationNumber(const std::vector<std::vector<double>>& occupationNumber)
{
    _occupationNumber = occupationNumber;
}

void Orbitals::set_orbitalEnergy(const std::vector<std::vector<double>>& orbitalEnergy)
{
    _orbitalEnergy = orbitalEnergy;
}


//----------------------------------------------------------------------------------------------------//
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

std::vector<std::vector<std::vector<double>>> Orbitals::getHomoCoefficients()
{
    std::vector<std::vector<std::vector<double>>> homoCoefficients(2, std::vector<std::vector<double>>(1));

    if (_homoLumoIndexes[0] == -1)
    {
        HOMO();
    }

    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);

    homoCoefficients[ALPHA][0] = _coefficients[ALPHA][_homoLumoIndexes[0]];
    homoCoefficients[BETA][0] = _coefficients[BETA][_homoLumoIndexes[0]];
    
    return homoCoefficients;
}

std::vector<double> Orbitals::getHomoEnergy()
{
    int ALPHA = static_cast<int>(SpinType::ALPHA);
    int BETA = static_cast<int>(SpinType::BETA);

    if (_homoLumoIndexes[0] == -1)
    {
        HOMO();
    }

    return std::vector<double>({ _orbitalEnergy[ALPHA][_homoLumoIndexes[0]], _orbitalEnergy[BETA][_homoLumoIndexes[0]] });
}

int Orbitals::getHomoIndex()
{
    if (_homoLumoIndexes[0] == -1)
    {
        HOMO();
    }

    return _homoLumoIndexes[0];
}

std::vector<std::vector<std::vector<double>>> Orbitals::getLumoCoefficients()
{
    std::vector<std::vector<std::vector<double>>> lumoCoefficients(2, std::vector<std::vector<double>>(1));

    if (_homoLumoIndexes[1] == -1)
    {
        LUMO();
    }

    const int ALPHA = static_cast<int>(SpinType::ALPHA);
    const int BETA = static_cast<int>(SpinType::BETA);

    lumoCoefficients[ALPHA][0] = _coefficients[ALPHA][_homoLumoIndexes[1]];
    lumoCoefficients[BETA][0] = _coefficients[BETA][_homoLumoIndexes[1]];

    return lumoCoefficients;
}

std::vector<double> Orbitals::getLumoEnergy()
{
    int ALPHA = static_cast<int>(SpinType::ALPHA);
    int BETA = static_cast<int>(SpinType::BETA);

    if (_homoLumoIndexes[1] == -1)
    {
        LUMO();
    }

    return std::vector<double>({ _orbitalEnergy[ALPHA][_homoLumoIndexes[1]], _orbitalEnergy[BETA][_homoLumoIndexes[1]] });
}

int Orbitals::getLumoIndex()
{
    if (_homoLumoIndexes[1] == -1)
    {
        LUMO();
    }

    return _homoLumoIndexes[1];
}

std::vector<std::vector<int>> Orbitals::getOccupiedOrbitalNumbers() const
{
    std::vector<std::vector<int>> occupiedOrbitalNumbers(2);

    for (int spin = 0; spin < 2; ++spin)
    {
        for (size_t i = 0; i < _occupationNumber[spin].size(); ++i)
        {
            if (_occupationNumber[spin][i] != 0.0)
            {
                occupiedOrbitalNumbers[spin].push_back(i + 1); // +1 because orbital numbers are 1-based indexing
            }
        }
    }

    return occupiedOrbitalNumbers;
}

std::vector<std::vector<int>> Orbitals::getVirtualOrbitalNumbers() const
{
    std::vector<std::vector<int>> virtualOrbitalNumbers(2);

    for (int spin = 0; spin < 2; ++spin)
    {
        for (size_t i = 0; i < _occupationNumber[spin].size(); ++i)
        {
            if (_occupationNumber[spin][i] == 0.0)
            {
                virtualOrbitalNumbers[spin].push_back(i + 1); // +1 because orbital numbers are 1-based indexing
            }
        }
    }
    
    return virtualOrbitalNumbers;
}

void Orbitals::getOccupiedAndVirtualOrbitalNumbers(std::vector<std::vector<int>>& occupiedOrbitalNumbers, std::vector<std::vector<int>>& virtualOrbitalNumbers) const
{
    occupiedOrbitalNumbers.resize(2);
    virtualOrbitalNumbers.resize(2);
    
    for (int spin = 0; spin < 2; ++spin)
    {
        for (size_t i = 0; i < _occupationNumber[spin].size(); ++i)
        {
            if (_occupationNumber[spin][i] == 0.0)
            {
                virtualOrbitalNumbers[spin].push_back(i + 1); // +1 because orbital numbers are 1-based indexing
            }
            else
            {
                occupiedOrbitalNumbers[spin].push_back(i + 1); // +1 because orbital numbers are 1-based indexing
            }
        }
    }
}

int Orbitals::getPrimitiveCenter(int i) const
{
    return _primitiveCenters[i];
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
                    sum += _vcgtf[np].ERICGTF(q.get_vcgtf()[nq],r.get_vcgtf()[nr],s.get_vcgtf()[ns]); 

    return sum;
}

double Orbitals::overlap(const int i, const int j, const SpinType spinType)
{
    if (spinType == SpinType::ALPHA_BETA)
    {
        print_error("Error in Orbitals::overlap(): spinType must be either ALPHA or BETA but not ALPHA_BETA.");
        std::exit(1);
    }

    int alpha = static_cast<int>(spinType);

    double sum = 0.0;

    #ifdef ENABLE_OMP
    #pragma omp parallel for reduction(+:sum)
    #endif
    for(size_t m = 0; m < _coefficients[alpha][i].size(); ++m)
    {
        for(size_t n = 0; n < _coefficients[alpha][j].size(); ++n)
        {
            sum += _coefficients[alpha][i][m] * _coefficients[alpha][j][n] * _vcgtf[m].overlapCGTF(_vcgtf[n]);
        }
    }

    return sum;
}

void Orbitals::printOverlap(const int i, const int j, const SpinType spinType)
{
    std::cout << "Overlap <" << i << "|" << j << "> = " << overlap(i, j, spinType) << std::endl;
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
    double sum=0.0;

    for (int n = 0; n < _numberOfAo; ++n)
    {
        for (int np = 0; np < _numberOfAo; ++np)
        {
            sum += _vcgtf[n].kineticCGTF(_vcgtf[np]);
        }
    }

    return sum;
}

std::vector<std::vector<std::vector<double>>> Orbitals::getIonicPotentialMatrix(const std::array<double, 3>& chargePosition, double charge, bool debug, bool printAOMatrix, bool printMOMatrix)
{
    // debug
    if (__debug_AOMatrix.size() == 0)
    {
        __debug_AOMatrix = std::vector<std::vector<double>>(_numberOfAo, std::vector<double>());
    }

    // Build ionic potential matrix in AO basis    
    std::vector<std::vector<double>> ionicMatrixAO = std::vector<std::vector<double>>(_numberOfAo, std::vector<double>());
    for (int i = 0; i < _numberOfAo; ++i)
    {
        ionicMatrixAO[i].resize(i + 1, 0.0);

        // debug
        __debug_AOMatrix[i].resize(i + 1, 0.0);
    }

    // Compute ionic potential matrix elements in AO basis
    for (int i = 0; i < _numberOfAo; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            ionicMatrixAO[i][j] = _vcgtf[i].ionicPotentialCGTF(_vcgtf[j], chargePosition, charge);
        }
    }

    // debug
    if (debug)
    {
        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        if (printAOMatrix) std::cout << "Ionic potential matrix in AO basis:" << std::endl;
        for (int ii = 0; ii < _numberOfAo; ++ii)
        {
            for (int jj = 0; jj <= ii; ++jj)
            {
                __debug_AOMatrix[ii][jj] += ionicMatrixAO[ii][jj];
                __debug_totalSumAO += (ii == jj ? ionicMatrixAO[ii][jj] : 2.0 * ionicMatrixAO[ii][jj]);
                if (printAOMatrix) std::cout << std::right << std::setw(17) << __debug_AOMatrix[ii][jj] << '\t';
            }
            if (printAOMatrix) std::cout << std::endl;
        }
        if (printAOMatrix) std::cout << std::defaultfloat << "Total sum of AO matrix elements: " << std::setprecision(10) << __debug_totalSumAO << std::endl << std::endl;
    }

    // debug
    if (__debug_MOMatrix.size() == 0)
    {
        __debug_MOMatrix = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(_numberOfMo, std::vector<double>()));
    }

    // Build ionic potential matrix in MO basis
    std::vector<std::vector<std::vector<double>>> ionicMatrixMO = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(_numberOfMo, std::vector<double>()));
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int i = 0; i < _numberOfMo; ++i)
        {
            ionicMatrixMO[spin][i].resize(i + 1, 0.0);

            // debug
            __debug_MOMatrix[spin][i].resize(i + 1, 0.0);
        }
    }

    // Compute ionic potential matrix elements in MO basis
    int spin, i, j;
    #ifdef ENABLE_OMP
    #pragma omp parallel for private(spin, i, j)
    #endif
    for (spin = 0; spin < 2; ++spin)
    {
        for (i = 0; i < _numberOfMo; ++i)
        {
            for (j = 0; j <= i; ++j)
            {
                double sum = 0.0;

                for (size_t m = 0; m < _coefficients[spin][i].size(); ++m)
                {
                    for (size_t n = 0; n < _coefficients[spin][j].size(); ++n)
                    {
                        sum += _coefficients[spin][i][m] * _coefficients[spin][j][n] * (n <= m ? ionicMatrixAO[m][n] : ionicMatrixAO[n][m]);
                    }
                }

                ionicMatrixMO[spin][i][j] = sum;
            }
        }
    }

    // debug
    if (debug)
    {
        for (spin = 0; spin < 2; ++spin)
        {
            std::cout << std::scientific;
            std::cout << std::setprecision(10);
            
            if (printMOMatrix) std::cout << "Ionic potential matrix in MO basis for " << to_string(static_cast<SpinType>(spin)) << " spin:" << std::endl;
            for (int ii = 0; ii < _numberOfMo; ++ii)
            {
                for (int jj = 0; jj <= ii; ++jj)
                {
                    __debug_MOMatrix[spin][ii][jj] += ionicMatrixMO[spin][ii][jj];
                    __debug_totalSumMO[spin] += (ii == jj ? ionicMatrixMO[spin][ii][jj] : 2.0 * ionicMatrixMO[spin][ii][jj]);

                    if (printMOMatrix) std::cout << std::right << std::setw(17) << __debug_MOMatrix[spin][ii][jj] << '\t';
                }
                if (printMOMatrix) std::cout << std::endl;
            }
            
            if (printMOMatrix) std::cout << std::defaultfloat <<  "Total sum of MO matrix elements for " << to_string(static_cast<SpinType>(spin)) << " spin: " << std::setprecision(10) << __debug_totalSumMO[spin] << std::endl;
        }
    }

    return ionicMatrixMO;
}

std::vector<std::vector<std::vector<std::vector<double>>>> Orbitals::getTripleOrbitalIntegralMatrix()
{
    // Build LRF matrix in AO basis
    std::vector<std::vector<std::vector<double>>> lrfMatrixAO = std::vector<std::vector<std::vector<double>>>(_numberOfAo, std::vector<std::vector<double>>());
    for (int i = 0; i < _numberOfAo; ++i)
    {
        lrfMatrixAO[i] = std::vector<std::vector<double>>(i + 1, std::vector<double>());

        for (int j = 0; j <= i; ++j)
        {
            lrfMatrixAO[i][j] = std::vector<double>(j + 1, 0.0);
        }
    }

    // Compute LRF matrix elements in AO basis
    for (int i = 0; i < _numberOfAo; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            for (int k = 0; k <= j; ++k)
            {
                lrfMatrixAO[i][j][k] = _vcgtf[i].overlap3CGTF(_vcgtf[j], _vcgtf[k]);
            }
        }
    }

    // debug
    /*
    std::cout << std::setprecision(10);
    for (int ii = 0; ii < _numberOfAo; ++ii)
    {
        std::cout << "ii = " << ii << std::endl;

        for (int jj = 0; jj <= ii; ++jj)
        {
            for (int kk = 0; kk <= jj; ++kk)
            {
                std::cout << std::right << std::setw(16) << lrfMatrixAO[ii][jj][kk] << ' ';
            }

            std::cout << std::endl;
        }
    }
    */

    // Build LRF matrix in MO basis
    std::vector<std::vector<std::vector<std::vector<double>>>> lrfMatrixMO = std::vector<std::vector<std::vector<std::vector<double>>>>(2, std::vector<std::vector<std::vector<double>>>(_numberOfMo, std::vector<std::vector<double>>()));
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int i = 0; i < _numberOfMo; ++i)
        {
            lrfMatrixMO[spin][i].resize(i + 1);

            for (int j = 0; j <= i; ++j)
            {
                lrfMatrixMO[spin][i][j].resize(j + 1);
            }
        }
    }

    // Compute LRF matrix elements in MO basis
    int spin, i, j, k;
    #ifdef ENABLE_OMP
    #pragma omp parallel for private(spin, i, j, k)
    #endif
    for (spin = 0; spin < 2; ++spin)
    {
        for (i = 0; i < _numberOfMo; ++i)
        {
            for (j = 0; j <= i; ++j)
            {
                for (k = 0; k <= j; ++k)
                {
                    double sum = 0.0;

                    for (size_t p = 0; p < _coefficients[spin][i].size(); ++p)
                    {
                        for (size_t q = 0; q < _coefficients[spin][j].size(); ++q)
                        {
                            for (size_t r = 0; r < _coefficients[spin][k].size(); ++r)
                            {
                                // Get the correct AO matrix element (considering symmetry)
                                std::array <size_t, 3> indices = { p, q, r };
                                std::sort(indices.begin(), indices.end(), std::greater<size_t>());
                                double lrfAOElement = lrfMatrixAO[indices[0]][indices[1]][indices[2]];

                                sum += _coefficients[spin][i][p] * _coefficients[spin][j][q] * _coefficients[spin][k][r] * lrfAOElement;
                            }
                        }
                    }

                    lrfMatrixMO[spin][i][j][k] = sum;
                }
            }
        }
    }

    return lrfMatrixMO;
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
    std::array<double, 3> C({ 0.0, 0.0, 0.0 });
    std::vector<int> l {ix, iy, iz};
    GTF m1(0.0, 1.0, C, l, _bino);
    std::vector<GTF> mbis (1,m1);
    CGTF m2(mbis);

    for(n=0;n<_numberOfAo;n++)
        for(ns=0;ns<_numberOfAo;ns++)
                sum += _vcgtf[n].gtf()[ns].overlap3GTF(m2.gtf()[0],_vcgtf[n].gtf()[ns]);

    return sum;
}

void Orbitals::NormaliseAllBasis()
{
    for(int k = 0; k < _numberOfAo; ++k)
    {
        _vcgtf[k].normaliseCGTF();
    }
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

    std::vector<SpinType> spins;
    if (_alphaAndBeta)
    {
        spins.push_back(SpinType::ALPHA);
    }
    else
    {
        spins.push_back(SpinType::ALPHA);
        spins.push_back(SpinType::BETA);
    }

    std::vector<double> evaluatedCgtfs(_vcgtf.size());
    evaluateCgtfsAtPoint(evaluatedCgtfs, x, y, z);

    for(SpinType spinType : spins)
    {
        int spin = static_cast<int>(spinType);

        for(int i = 0; i < _numberOfMo; ++i)
        {
            if(_coefficients[spin][i].size() != _vcgtf.size())
            {
                std::cout<<"Error in Orbitals::func(): there are " << _coefficients[spin][i].size() << " coefficients for " << _vcgtf.size() << " CGTFs." << std::endl;
                std::exit(1);
            }

            for(int k =0 ; k < _numberOfMo; ++k)
            {
                if(std::abs(_coefficients[spin][i][k]) > 1e-10)
                {
                    r += _coefficients[spin][i][k] * evaluatedCgtfs[k];
                }
            }
        }
    }

    return r;
}

void Orbitals::HOMO()
{
    if(_numberOfAlphaElectrons >= _numberOfBetaElectrons)
    {
        _homoLumoIndexes[0] = _numberOfAlphaElectrons - 1; // 0-based index
    }
    else
    {
        _homoLumoIndexes[0] = _numberOfBetaElectrons - 1; // 0-based index
    }
}

void Orbitals::LUMO()
{
    if (_homoLumoIndexes[0] == -1)
    {
        HOMO();
    }

    _homoLumoIndexes[1] = _homoLumoIndexes[0] + 1;

    if(_homoLumoIndexes[1] >= _numberOfMo)
    {
        print_error("Error in Orbitals::LUMO(): LUMO is not available because its orbital index is out of bounds.");

        std::exit(1);
    }
}

std::vector<std::vector<double>> Orbitals::get_S()
{
    int i,j;
    std::vector<std::vector<double>> S (_numberOfAo, std::vector<double> (_numberOfAo,0.0));

    #ifdef ENABLE_OMP
    #pragma omp parallel for private(i, j)
    #endif
    for(i=0; i<_numberOfAo; i++)
    {
        for(j=i; j<_numberOfAo; j++)
        {
            S[i][j]=S[j][i]=_vcgtf[i].overlapCGTF(_vcgtf[j]);
        }
    }

    return S;
}

std::vector<double> Orbitals::get_f(int orb, int alpha)
{
    int i;
    std::vector<std::vector<double>> S =get_S();
    size_t nu,xi;
    std::vector<double> f(_numberOfAtoms,0.0);

    #ifdef ENABLE_OMP
    #pragma omp parallel for private(i,nu,xi)
    #endif
    for(i=0; i<_numberOfAtoms; i++)
        for(nu=0; nu<_coefficients[alpha][orb].size(); nu++)
        {
            if(i+1 == _primitiveCenters[nu])
                f[i]+=_coefficients[alpha][orb][nu]*_coefficients[alpha][orb][nu];

            for(xi=0; xi<_coefficients[alpha][orb].size(); xi++)
                if(xi!=nu && i+1 == _primitiveCenters[nu])
                    f[i]+=_coefficients[alpha][orb][xi]*_coefficients[alpha][orb][nu];
        }

    return f;
}

void Orbitals::get_f(int alpha)
{
    if (_homoLumoIndexes[0] == -1)
    {
        HOMO();
    }

    if (_homoLumoIndexes[1] == -1)
    {
        LUMO();
    }

    std::vector<std::vector<double>> S = get_S();
    std::vector<std::vector<double>> f(2, std::vector<double>(_numberOfAtoms, 0.0));

    for(int i = 0; i < _numberOfAtoms; ++i)
    {
        for(size_t j = 0; j < 2; ++j)
        {
            for(size_t nu = 0; nu < _coefficients[alpha][_homoLumoIndexes[j]].size(); ++nu)
            {
                if(i + 1 == _primitiveCenters[nu])
                {
                    f[j][i] += _coefficients[alpha][_homoLumoIndexes[j]][nu] * _coefficients[alpha][_homoLumoIndexes[j]][nu];
                }

                for(size_t xi = 0; xi < _coefficients[alpha][_homoLumoIndexes[j]].size(); ++xi)
                {
                    if(xi != nu && i + 1 == _primitiveCenters[nu])
                    {
                        f[j][i] += _coefficients[alpha][_homoLumoIndexes[j]][xi] * _coefficients[alpha][_homoLumoIndexes[j]][nu] * S[xi][nu];
                    }
                }
            }
        }
    }

    _all_f = f;
}

void Orbitals::init_homoLumoIndexes()
{
    HOMO();
    LUMO();
}

void Orbitals::printDescriptors()
{
    init_homoLumoIndexes();
    std::cout<<"end HOMOLUMO"<<std::endl;
    get_f();
    std::cout<<"end get_f"<<std::endl;
    _descriptors.set_mu_fk_data(_all_f, getHomoEnergy()[static_cast<int>(SpinType::ALPHA)], getLumoEnergy()[static_cast<int>(SpinType::ALPHA)]);
    _descriptors.compute_all();
    std::cout<<_descriptors<<std::endl;
}

void Orbitals::printDescriptors(int homoIndex, int lumoIndex)
{
    set_homoLumoIndexes(homoIndex, lumoIndex);
    get_f();
    _descriptors.set_mu_fk_data(_all_f, getHomoEnergy()[static_cast<int>(SpinType::ALPHA)], getLumoEnergy()[static_cast<int>(SpinType::ALPHA)]);
    _descriptors.compute_all();
    std::cout<<_descriptors<<std::endl;
}

double operator*(const Orbitals& a, const std::vector<double>& coord)
{
    double r=1.0;
    for(size_t i=1; i<coord.size(); i++)
        r*=a.func(coord[0],coord[1],coord[2]);
    
    return r;
}

std::ostream& operator<<(std::ostream& stream, const Orbitals& orbitals)
{
    stream << std::scientific;
    stream << std::setprecision(10);
    stream << std::setw(20);
    stream << std::left << std::setw(20) << "Coef CGTF"
                   << std::setw(20) << "Coef GTF"
                   << std::setw(20) << "Exp"
                   << std::setw(5)  << "Lx"
                   << std::setw(5)  << "Ly"
                   << std::setw(5)  << "Lz"
                   << std::setw(20) << "x"
                   << std::setw(20) << "y"
                   << std::setw(20) << "z" << std::endl;
    
    for(int i = 0; i < orbitals.get_numberOfAo(); ++i)
    {
        stream << std::left << orbitals.get_vcgtf()[i] << std::endl;
    }

    int n = 2;
    if(orbitals._alphaAndBeta)
    {
        n = 1;
    }

    for(int j = 0; j < orbitals.get_numberOfMo(); ++j)
    {
        for(int i = 0; i < n; ++i)
        {
            if(i == 0)
            {
                stream << "Alpha, Occ= " << orbitals._occupationNumber[i][j] << " num= " << j << std::endl;
            }
            else
            {
                stream << "Beta, Occ= " << orbitals._occupationNumber[i][j] << " num= " << j << std::endl;
            }

            for(size_t k = 0; k < orbitals._vcgtf.size(); ++k)
            {
                stream << ' ' << std::left << k << ' ' << orbitals.get_coefficients()[i][j][k] << std::endl;
            }
        }
    }

    stream << std::defaultfloat;

    return stream;
}

Grid Orbitals::makeGrid(const Domain& d)
{
    Grid g;
    g.set_structure(_struct);
    g.set_domain(d);
    g.reset();
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
    for(int i=0;i<d.get_N1();i++)
    {
        for(int j=0;j<d.get_N2();j++)
        {
            for(int k=0;k<d.get_N3();k++)
            {
                double rho=density(d.x(i,j,k), d.y(i,j,k), d.z(i,j,k));
                g.set_Vijkl(rho,i,j,k,0);
            }
        }
    }
    return g;
}

double Orbitals::density(double x, double y, double z) const
{
    double rho = 0.0;
    SpinType spinType = _alphaAndBeta ? SpinType::ALPHA : SpinType::ALPHA_BETA;

    // Evaluate all CGTFs at the given point (x, y, z) once and reuse the values for each orbital
    std::vector<double> evaluatedCgtfs(_vcgtf.size());
    evaluateCgtfsAtPoint(evaluatedCgtfs, x, y, z);

    for(int orbitalIndex = 0; orbitalIndex < get_numberOfMo(); ++orbitalIndex)
    {
        rho += density(orbitalIndex, spinType, evaluatedCgtfs);
    }

    return rho;
}

double Orbitals::density(int orbitalIndex, SpinType spinType, double x, double y, double z) const
{
    // Evaluate all CGTFs at the given point (x, y, z) once and reuse the values for each orbital
    std::vector<double> evaluatedCgtfs(_vcgtf.size());
    evaluateCgtfsAtPoint(evaluatedCgtfs, x, y, z);

    return density(orbitalIndex, spinType, evaluatedCgtfs);
}

double Orbitals::density(const std::vector<int>& orbitalIndexes, const std::vector<SpinType>& orbitalSpins, double x, double y, double z) const
{
    if (orbitalIndexes.size() != orbitalSpins.size())
    {
        print_error("Error in Orbitals::density(): orbitalIndexes and orbitalSpins vectors must have the same size.");
        std::exit(1);
    }

    double rho = 0.0;

    // Evaluate all CGTFs at the given point (x, y, z) once and reuse the values for each orbital
    std::vector<double> evaluatedCgtfs(_vcgtf.size());
    evaluateCgtfsAtPoint(evaluatedCgtfs, x, y, z);

    for (size_t i = 0; i < orbitalIndexes.size(); ++i)
    {
        rho += density(orbitalIndexes[i], orbitalSpins[i], evaluatedCgtfs);
    }

    return rho;
}

Grid Orbitals::makeOrbGrid(const Domain& domain, const std::vector<int>& orbitalsNumbers, const std::vector<SpinType>& orbitalsSpins, bool showProgress)
{
    Grid g;
    g.set_structure(_struct);
    g.set_domain(domain);
    g.reset();

    int N1 = domain.get_N1();
    int N2 = domain.get_N2();
    int N3 = domain.get_N3();

    const int nbStepsTotal = N1 * N2 * N3;
    std::atomic<int> progress(0);
    int lastProgress = -1;

    // Afficher la barre de progression à 0% dès le début
    if(showProgress)
    {
        print_progressBar(0, nbStepsTotal, lastProgress);
    }

    #ifdef ENABLE_OMP
    #pragma omp parallel
    #endif
    {
        #ifdef ENABLE_OMP
        #pragma omp for collapse(2)
        #endif
        for (int i = 0; i < N1; ++i)
        {
            for(int j = 0; j < N2; ++j)
            {
                for(int k = 0; k < N3; ++k)
                {
                    std::vector<double> phy = phis(domain.x(i, j, k), domain.y(i, j, k), domain.z(i, j, k), orbitalsNumbers, orbitalsSpins);

                    for(int l = 0; l < domain.get_Nval(); ++l)
                    {
                        g.set_Vijkl(phy[l], i, j, k, l);
                    }
                }
                
                if(showProgress)
                {
                    // Mise à jour à chaque itération de N2 pour un affichage plus fluide
                    int currentStep = progress.fetch_add(N3) + N3;
                    
                    #ifdef ENABLE_OMP
                    #pragma omp critical
                    #endif
                    {
                        print_progressBar(currentStep, nbStepsTotal, lastProgress);
                    }
                }
            }
        }
    }
    
    if(showProgress)
    {
        std::cout << std::endl;
    }

    return g;
}

std::vector<double> Orbitals::phis(double x, double y, double z, const std::vector<int>& orbitalIndexes, const std::vector<SpinType>& orbitalSpins)
{
    std::vector<double> values(orbitalIndexes.size(), 0.0);

    std::vector<double> evaluatedCgtfs(_vcgtf.size());
    evaluateCgtfsAtPoint(evaluatedCgtfs, x, y, z);

    for (size_t i = 0; i < orbitalIndexes.size(); ++i)
    {
        int orbitalIndex = orbitalIndexes[i];
        int spin = static_cast<int>(orbitalSpins[i]);

        for(size_t k = 0; k < _vcgtf.size(); ++k)
        {
            values[i] += _coefficients[spin][orbitalIndex][k] * evaluatedCgtfs[k];
        }
    }

    return values;
}
//epsilon=0 for Becke, epsilon =2.87e-5 for Savin. see Can. J. Chem. Vol. 74,1996 page 1088.
double Orbitals::ELF(double x, double y, double z, double epsilon)
{
    double rho = 0.0;
    double sphi = 0.0;
    double cf = 3.0 / 10.0 * std::pow(3 * M_PI * M_PI, 2.0 / 3);
    int n = _alphaAndBeta ? 1 : 2;

    std::vector<double> v(_vcgtf.size());
    evaluateCgtfsAtPoint(v, x, y, z);

    std::vector<double> A(_vcgtf.size());
    std::vector< std::vector<double> > vg(3,A);
    for(size_t k=0;k<_vcgtf.size();k++)
    {
        for(int c=0;c<3;c++)
            vg[c][k]=_vcgtf[k].grad_CGTF(x,y,z,c);
    }

    double v1[3]={0,0,0};
    for(int j=0; j<get_numberOfMo(); j++)
    for(int i=0; i<n; i++)
    {
        if(get_occupationNumber()[i][j]>1e-10)
        {
            double phi=0;
                for(size_t k=0; k<_vcgtf.size(); k++)
                phi += _coefficients[i][j][k]*v[k];
            rho+=get_occupationNumber()[i][j] * phi*phi;

            for(int c=0;c<3;c++)
            {
                double dp=0;
                for(size_t k=0; k<_vcgtf.size(); k++)
                    dp += _coefficients[i][j][k]*vg[c][k];
                sphi += get_occupationNumber()[i][j] * dp*dp;
                v1[c] += get_occupationNumber()[i][j] * phi*dp;
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
    g.set_structure(_struct);
    g.set_domain(d);
    g.reset();
#ifdef ENABLE_OMP
#pragma omp parallel
#endif
    for(int i=0;i<d.get_N1();i++)
    {
        for(int j=0;j<d.get_N2();j++)
        {
            for(int k=0;k<d.get_N3();k++)
            {
                double elf=ELF(d.x(i,j,k), d.y(i,j,k), d.z(i,j,k), epsilon);
                g.set_Vijkl(elf,i,j,k,0);
            }
        }
    }
    return g;
}
