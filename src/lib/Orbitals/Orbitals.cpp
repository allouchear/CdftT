#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "../Basis/GTF.h"
#include "../Basis/CGTF.h"
#include "../Common/PeriodicTable.h"
#include "../Orbitals/Orbitals.h"
#include <Utils/Enums.hpp>
#include "../Utils/FCHK.h"
#include "../Utils/LM.h"
#include "../Utils/LOG.h"
#include "../Utils/MOLDENGAB.h"
#include "../Utils/Utils.h"
#include "../Utils/WFX.h"


//----------------------------------------------------------------------------------------------------//
// CONSTRUCTORS
//----------------------------------------------------------------------------------------------------//

Orbitals::Orbitals():
    _vcgtf(),
    _coefficients(),
    _numberOfAo(0),
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
    _numOrb(),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(),
    _descriptors(),
    _vcgtfNonNormalise(),
    _numberOfGtf(0),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{ }

Orbitals::Orbitals(WFX& wfxParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _coefficients(2),
    _numberOfAo(0),
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
    _numOrb(),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _vcgtfNonNormalise(),
    _numberOfGtf(0),
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
    _numOrb = std::vector<int> (2,0);
    _energy = wfxParser.Energy();
    _occupationNumber = wfxParser.Molecular_Orbital_Occupation_Numbers();
    _alphaAndBeta = wfxParser.AlphaAndBeta();
    _descriptors = Descriptors(wfxParser, periodicTable);
    _numberOfAo = _vcgtf.size();
    _vcgtfNonNormalise = _vcgtf;
}

Orbitals::Orbitals(FCHK& fchkParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _coefficients(2),
    _numberOfAo(0),
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
    _numOrb(),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _vcgtfNonNormalise(),
    _numberOfGtf(0),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{
    _numberOfMo = fchkParser.NumberOfBasisFunctions();
    _coordinates = fchkParser.CurrentCartesianCoordinates();
    _energy = fchkParser.TotalEnergy();
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

        exit(1);
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

    _numOrb = std::vector<int>(2, 0);
    _descriptors = Descriptors(fchkParser, periodicTable);

    for(size_t i = 0; i < _vcgtf.size(); i++)
    {
        _primitiveCenters.push_back(_vcgtf[i].NumCenter());
    }

    _vcgtfNonNormalise = _vcgtf;
    NormaliseAllBasis();

    for(size_t i = 0; i < _vcgtf.size(); i++)
    {
        _numberOfGtf += _vcgtf[i].numberOfFunctions();
    }
}

Orbitals::Orbitals(MOLDENGAB& moldengabParser, Binomial& bino, const PeriodicTable& periodicTable):
    _vcgtf(),
    _coefficients(2),
    _numberOfAo(0),
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
    _numOrb(),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _vcgtfNonNormalise(),
    _numberOfGtf(0),
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

    _numOrb = std::vector<int>(2, 0);

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
    std::vector<std::vector<double>> coordinatesForShells;
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
                            coord_[c] = coordinatesForShells[nS][c];
                            l_[c] = l[c][m][n];
                        }

                        GTF gtf (primitiveExponents[kPrimitive+ip], 1, coord_, l_, _bino);
                        _vcgtf[kOrb].push_back(gtf);
                        _vcgtf[kOrb].setCoef(FactorCoefs[kOrb]*CgtfCoefs[kPrimitive+ip]*coefs[m][n]);
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

        exit(1);
    }

    _vcgtfNonNormalise = _vcgtf;

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
    _coefficients(2),
    _numberOfAo(0),
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
    _numOrb(),
    _occupationNumber(),
    _alphaAndBeta(false),
    _bino(bino),
    _descriptors(),
    _vcgtfNonNormalise(),
    _numberOfGtf(0),
    _energy(0.0),
    _coordinates(),
    _mixte(false)
{
    _struct = Structure(logParser, periodicTable);
    _vcgtfNonNormalise = std::vector<CGTF> ();
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
    _numOrb = std::vector<int> (2,0);
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
                            coord_[c] = coordinatesForShells[nS][c];
                            l_[c] = l[c][m][n];
                        }

                        GTF gtf(primitiveExponents[kPrimitive + ip], 1.0, coord_, l_, _bino);
                        _vcgtf[kOrb].push_back(gtf);
                        _vcgtf[kOrb].setCoef(FactorCoefs[kOrb] * CgtfSpCoefs[kPrimitive + ip] * coefs[m][n]);
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
        cout << "Error : There are " << _vcgtf.size() << " CGTFs for " << _numberOfMo << " basis in file." << endl;
        cout << "Please check your file." << endl;

        exit(1);
    }

    if(logParser.NumberOfBasisFunctions() < _numberOfMo)
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

    _vcgtfNonNormalise = _vcgtf;

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
// GETTERS
//----------------------------------------------------------------------------------------------------//

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
// OTHER PUBLIC METHODS
//----------------------------------------------------------------------------------------------------//

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
        exit(1);
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

std::vector<std::vector<double>> displayM;
double totalSumM = 0.0;

std::vector<std::vector<std::vector<double>>> Orbitals::getIonicPotentialMatrix(const std::array<double, 3>& chargePosition, double charge, bool debug)
{
    // debug
    if (displayM.size() == 0)
    {
        displayM = std::vector<std::vector<double>>(_numberOfAo, std::vector<double>());
    }

    // Build ionic potential matrix in AO basis    
    std::vector<std::vector<double>> ionicMatrixAO = std::vector<std::vector<double>>(_numberOfAo, std::vector<double>());
    for (int i = 0; i < _numberOfAo; ++i)
    {
        ionicMatrixAO[i].resize(i + 1, 0.0);

        // debug
        displayM[i].resize(i + 1, 0.0);
    }

    // Compute ionic potential matrix elements in AO basis
    for (int i = 0; i < _numberOfAo; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            ionicMatrixAO[i][j] = _vcgtf[i].ionicPotentialCGTF(_vcgtf[j], chargePosition, charge);

            // debug
            // ionicMatrixAO[i][j] = _vcgtf[i].ionicPotentialCGTF(_vcgtf[j], chargePosition, charge, (i == 9 && j == 0));
        }
    }

    // debug
    if (debug)
    {
        for (int ii = 0; ii < _numberOfAo; ++ii)
        {
            for (int jj = 0; jj <= ii; ++jj)
            {
                displayM[ii][jj] += ionicMatrixAO[ii][jj];
                totalSumM += (ii == jj ? ionicMatrixAO[ii][jj] : 2.0 * ionicMatrixAO[ii][jj]);
                std::cout << std::setprecision(6) << std::setw(10) << displayM[ii][jj] << '\t';
            }
            std::cout << std::endl;
        }
        std::cout << "Total sum: " << totalSumM << std::endl << std::endl;
    }

    // Build ionic potential matrix in MO basis
    std::vector<std::vector<std::vector<double>>> ionicMatrixMO = std::vector<std::vector<std::vector<double>>>(2, std::vector<std::vector<double>>(_numberOfMo, std::vector<double>()));
    for (int spin = 0; spin < 2; ++spin)
    {
        for (int i = 0; i < _numberOfMo; ++i)
        {
            ionicMatrixMO[spin][i].resize(i + 1, 0.0);
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

    return ionicMatrixMO;
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

    if(_alphaAndBeta)
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
                if(std::abs(_coefficients[i][j][k])>1e-10)
                    r+=_coefficients[i][j][k] * _vcgtf[k].func(x,y,z);
            }
        }
    }

    return r;
}

void Orbitals::HOMO()
{
    if(_numberOfAlphaElectrons>=_numberOfBetaElectrons)
        _numOrb[0] = _numberOfAlphaElectrons-1;
    else
        _numOrb[0] = _numberOfBetaElectrons-1;
}

void Orbitals::LUMO()
{
    _numOrb[1] = _numOrb[0]+1;
    if(_numOrb[1] +1 >_numberOfMo)
    {
        cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        cerr<<" Lumo is not available in your file orbitals file"<<endl;
        cerr<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
        exit(1);
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
    std::vector<std::vector<double>> S=get_S();
    int i;
    size_t j,nu,xi;
    std::vector<double> V(_numberOfAtoms,0.0);
    std::vector<std::vector<double>> f(_numOrb.size(), V);

    for(i=0; i<_numberOfAtoms; i++)
        for(j=0; j<_numOrb.size(); j++)
            for(nu=0; nu<_coefficients[alpha][_numOrb[j]].size(); nu++)
            {
                if(i+1 == _primitiveCenters[nu])
                {
                    f[j][i]+=_coefficients[alpha][_numOrb[j]][nu]*_coefficients[alpha][_numOrb[j]][nu];
                }

                for(xi=0; xi<_coefficients[alpha][_numOrb[j]].size(); xi++)
                    if(xi!=nu && i+1 == _primitiveCenters[nu])
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
    cout<<"end HOMOLUMO"<<endl;
    get_f();
    cout<<"end get_f"<<endl;
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

double operator*(const Orbitals& a, const std::vector<double>& coord)
{
    double r=1.0;
    for(size_t i=1; i<coord.size(); i++)
        r*=a.func(coord[0],coord[1],coord[2]);
    
    return r;
}

std::ostream& operator<<(std::ostream& stream, Orbitals& orbitals)
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
double Orbitals::density(double x, double y, double z)
{
    double rho = 0.0;
    int n = _alphaAndBeta ? 1 : 2;

    std::vector<double> v(_vcgtf.size());
    for(size_t k = 0; k < _vcgtf.size(); ++k)
    {
        v[k] = _vcgtf[k].func(x,y,z);
    }

    for(int j = 0; j < get_numberOfMo(); ++j)
    {
        for(int i = 0; i < n; ++i)
        {
            if(get_occupationNumber()[i][j] > 1e-10)
            {
                double phi = 0;
                
                for(size_t k = 0; k < _vcgtf.size(); ++k)
                {
                    phi += _coefficients[i][j][k] * v[k];
                }

                rho += get_occupationNumber()[i][j] * phi * phi;
            }
        }
    }

    return rho;
}

Grid Orbitals::makeOrbGrid(const Domain& d, const std::vector<int>& nums, const std::vector<int>& typesSpin)
{
    Grid g;
    g.set_structure(_struct);
    g.set_domain(d);
    g.reset();

    #ifdef ENABLE_OMP
    #pragma omp parallel
    #endif
    for(int i = 0; i < d.get_N1() ; ++i)
    {
        for(int j = 0; j < d.get_N2(); ++j)
        {
            for(int k = 0; k < d.get_N3(); ++k)
            {
                std::vector<double> phy = phis(d.x(i,j,k), d.y(i,j,k), d.z(i,j,k), nums, typesSpin);

                for(int l = 0; l < d.get_Nval(); ++l)
                {
                    g.set_Vijkl(phy[l],i,j,k,l);
                }
            }
        }
    }

    return g;
}

std::vector<double> Orbitals::phis(double x, double y, double z, const std::vector<int>& nums, const std::vector<int>& typesSpin)
{
    std::vector<double> v(_vcgtf.size());
    for(size_t k = 0; k < _vcgtf.size(); ++k)
    {
        v[k] = _vcgtf[k].func(x, y, z);
    }

    std::vector<double> values(nums.size(), 0);
    for(size_t jj = 0; jj < nums.size(); ++jj)
    {
        int j = nums[jj];
        int i = typesSpin[jj];
        values[jj] = 0;

        for(size_t k = 0; k < _vcgtf.size(); ++k)
        {
            values[jj] += _coefficients[i][j][k] * v[k];
        }
    }

    return values;
}
//epsilon=0 for Becke, epsilon =2.87e-5 for Savin. see Can. J. Chem. Vol. 74,1996 page 1088.
double Orbitals::ELF(const double& x, const double& y, const double& z, double epsilon)
{
    double rho = 0.0;
    double sphi = 0.0;
    double cf = 3.0 / 10.0 * std::pow(3 * M_PI * M_PI, 2.0 / 3);
    int n = _alphaAndBeta ? 1 : 2;

    std::vector<double> v(_vcgtf.size());
        for(size_t k=0; k<_vcgtf.size(); k++)
        v[k] = _vcgtf[k].func(x,y,z);

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
