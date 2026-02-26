#include <array>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <Common/Constants.h>
#include <Utils/MOLDENGAB.h>
#include <Utils/Utils.h>


MOLDENGAB::MOLDENGAB() :
    _symbols(),
    _atomic_numbers(),
    _coordinates(),
    _shell_types(),
    _L_types(),
    _exposants(),
    _number_of_gtf(),
    _num_center(),
    _cgtf_coefs(),
    _factor_coefs(),
    _MO_energy(),
    _MO_coefs(),
    _occupation(),
    _spin_types(),
    _basis_or_gto("None"),
    _number_of_atoms(0),
    _number_of_MO_coefs(0),
    _number_of_MO(0),
    _alpha_and_beta(true),
    _n_at_basis(),
    _alpha_occupation(),
    _beta_occupation(),
    _alpha_MO_coefs(),
    _beta_MO_coefs(),
    _alpha_energies(),
    _beta_energies(),
    _format("None"),
    _cart_sphe("none"),
    _mixte(false)
{ }

MOLDENGAB::MOLDENGAB(std::ifstream& moldenGabFile) :
    _symbols(),
    _atomic_numbers(),
    _coordinates(),
    _shell_types(),
    _L_types(),
    _exposants(),
    _number_of_gtf(),
    _num_center(),
    _cgtf_coefs(),
    _factor_coefs(),
    _MO_energy(),
    _MO_coefs(),
    _occupation(),
    _spin_types(),
    _basis_or_gto("None"),
    _number_of_atoms(0),
    _number_of_MO_coefs(0),
    _number_of_MO(0),
    _alpha_and_beta(true),
    _n_at_basis(),
    _alpha_occupation(),
    _beta_occupation(),
    _alpha_MO_coefs(),
    _beta_MO_coefs(),
    _alpha_energies(),
    _beta_energies(),
    _format("None"),
    _cart_sphe("none"),
    _mixte(false)
{

    moldenGabFile.clear();
    moldenGabFile.seekg(0,moldenGabFile.beg);
    std::string p;
    getline(moldenGabFile,p);

    if(p.find("[Molden Format]") != std::string::npos)
    {
        _format = "molden";
        _basis_or_gto = "[GTO]";
    }
    else if(p.find("[Gabedit Format]") != std::string::npos)
    {
        _format = "gabedit";
        _basis_or_gto = "[Basis]";

        if(p.find("Sphe") != std::string::npos)
        {
            _cart_sphe = "sphe";
        }
        else if(p.find("Cart") != std::string::npos)
        {
            _cart_sphe = "cart";
        }
        else
        {
            std::stringstream errorMessage;
            errorMessage << "Error in MOLDENGAB::MOLDENGAB(): can't recognize data format (sphe/cart)." << std::endl;
            errorMessage << "Please check your file.";
            print_error(errorMessage.str());

            std::exit(1);
        }
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Error in MOLDENGAB::MOLDENGAB(): can't recognize file format." << std::endl;
        errorMessage << "Please check your file.";
        print_error(errorMessage.str());

        std::exit(1);
    }

    read_atom_data(moldenGabFile);

    _number_of_atoms = _atomic_numbers.size();
    _n_at_basis = std::vector<int>(_number_of_atoms);
    
    read_basis_data(moldenGabFile);
    read_MO_data(moldenGabFile);

    if(_alpha_and_beta)
    {
        _alpha_energies = _beta_energies = _MO_energy;
        _alpha_occupation = _beta_occupation = _occupation;
        _alpha_MO_coefs = _beta_MO_coefs = _MO_coefs;
    }
    else
    {
        for(int i = 0; i < _number_of_MO; ++i)
        {
            if(_spin_types[i] == "Alpha")
            {
                _alpha_energies.push_back(_MO_energy[i]);
                _alpha_occupation.push_back(_occupation[i]);
                _alpha_MO_coefs.push_back(_MO_coefs[i]);
            }
            
            else if(_spin_types[i] == "Beta")
            {
                _beta_energies.push_back(_MO_energy[i]);
                _beta_occupation.push_back(_occupation[i]);
                _beta_MO_coefs.push_back(_MO_coefs[i]);
            }
        }
    }
}

void MOLDENGAB::read_atom_data(std::ifstream& moldenGabFile)
{
    std::string line;
    std::string word;

    std::string coord_type;

    int atomicNumber;
    std::array<double, 3> coords;

    long int position = LocaliseDataMolGabBefore(moldenGabFile,"[Atoms]");

    if(position == -1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in MOLDENGAB::read_atom_data(): Atoms data not found." << std::endl;
        errorMessage << "Data required, please check your file.";
        print_error(errorMessage.str());

        std::exit(1);
    }

    moldenGabFile.seekg(position);
    getline(moldenGabFile, line);
    std::stringstream sstream(line);
    sstream >> word;
    sstream >> coord_type;
    
    getline(moldenGabFile, line);
    do
    {
        std::stringstream sstream2(line);

        // Atom symbol
        sstream2 >> word;    
        _symbols.push_back(word);

        // Atom number in the system (not stored)
        sstream2 >> word;

        // Atomic number
        sstream2 >> atomicNumber;
        _atomic_numbers.push_back(atomicNumber);

        // Coordinates
        sstream2 >> coords[0];
        sstream2 >> coords[1];
        sstream2 >> coords[2];

        if (coord_type == "Angs")
        {
            for(int i = 0; i < 3; ++i)
            {
                coords[i] *= Constants::ANGSTROM_TO_BOHR_RADIUS;
            }
        }

        _coordinates.push_back(coords);

        getline(moldenGabFile, line);
    } while(line.find("[") == std::string::npos);
}

void MOLDENGAB::read_basis_data(std::ifstream& moldenGabFile)
{
    long int position = LocaliseDataMolGab(moldenGabFile, _basis_or_gto);

    if(position == -1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in MOLDENGAB::read_basis_data(): Basis (GTO) data not found." << std::endl;
        errorMessage << "Data required, please check your file.";
        print_error(errorMessage.str());

        std::exit(1);
    }

    moldenGabFile.seekg(position);
    int nbAtoms = _atomic_numbers.size();

    for(int i = 0; i < nbAtoms; ++i)
    {
        read_one_basis_data(moldenGabFile);
    }

    int m = _shell_types.size();
    
    long int posd, posf, posg;

    if(_format == "molden")
    {
        posd = LocaliseDataMolGab(moldenGabFile, "[5D]");
        if (posd == -1)
        {
            posd = LocaliseDataMolGab(moldenGabFile, "[5d]");
        }

        posf = LocaliseDataMolGab(moldenGabFile, "[7F]");
        if (posf == -1)
        {
            posf = LocaliseDataMolGab(moldenGabFile, "[7f]");
        }

        posg = LocaliseDataMolGab(moldenGabFile, "[9G]");
        if (posg == -1)
        {
            posg = LocaliseDataMolGab(moldenGabFile, "[9g]");
        }

        if((posd == -1 && (posf !=- 1 || posg != 1))
           || (posf == -1 && (posd !=- 1 || posg != 1 ))
           || (posg == -1 && (posf !=- 1 || posd != 1)))
        {
            _mixte = true;
        }
    }
    else if(_format == "gabedit" && _cart_sphe == "sphe")
    {
        posd = posf = posg = 0;
    }
    else
    {
        posd = posf = posg = -1;
    }

    for(int i = 0; i < m; ++i)
    {
        if(_shell_types[i] == "s" || _shell_types[i] == "S")
        {
            _L_types.push_back(0);
            _number_of_MO_coefs += 1;
        }
        else if(_shell_types[i] == "p" || _shell_types[i] == "P")
        {
            _L_types.push_back(1);
            _number_of_MO_coefs += 3;
        }
        else if((_shell_types[i] == "d" || _shell_types[i] == "D") && posd == -1)
        {
            _L_types.push_back(2);
            _number_of_MO_coefs += 6;
        }
        else if((_shell_types[i] == "d" || _shell_types[i] == "D") && posd != -1)
        {
            _L_types.push_back(-2);
            _number_of_MO_coefs += 5;
        }
        else if((_shell_types[i] == "f" || _shell_types[i] == "F") && posf == -1)
        {
            _L_types.push_back(3);
            _number_of_MO_coefs += 10;
        }
        else if((_shell_types[i] == "f" || _shell_types[i] == "F") && posf != -1)
        {
            _L_types.push_back(-3);
            _number_of_MO_coefs += 7;
        }
        else if((_shell_types[i] == "g" || _shell_types[i] == "G") && posg == -1)
        {
            _L_types.push_back(4);
            _number_of_MO_coefs += 15;
        }
        else if((_shell_types[i] == "g" || _shell_types[i] == "G") && posg != -1)
        {
            _L_types.push_back(-4);
            _number_of_MO_coefs += 9;
        }
        else if(_format == "molden")
        {
            int N = int(toupper(_shell_types[i][0])) - int('F') + 3;
            _L_types.push_back(-N);
            _number_of_MO_coefs += 2 * std::abs(N) + 1;
        }
        else if(_cart_sphe == "sphe")
        {
            int N = int(toupper(_shell_types[i][0])) - int('F') + 3;
            _L_types.push_back(-N);
            _number_of_MO_coefs += 2 * std::abs(N) + 1;
        }
        else if(_cart_sphe == "cart")
        {
            int N = int(toupper(_shell_types[i][0])) - int('F') + 3;
            _L_types.push_back(N);
            _number_of_MO_coefs += 2 * N + 1;
        }
        else
        {
            std::stringstream errorMessage;
            errorMessage << "Error in MOLDENGAB::read_basis_data(): shell type not recognized." << std::endl;
            errorMessage << "Please check your file.";
            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    _number_of_MO=_number_of_MO_coefs;
}

void MOLDENGAB::read_one_basis_data(std::istream& moldenGabFile)
{
    std::string line;
    std::string shellType;

    int centerNumber;
    int nbGTF;
    int nbGTFCenter = 0;

    double factorCoeff;
    double exposant, coeff;

    getline(moldenGabFile, line);
    std::stringstream sstream(line);
    sstream >> centerNumber;

    getline(moldenGabFile, line);
    do
    {
        _num_center.push_back(centerNumber);

        std::stringstream sstream2(line);
        sstream2 >> shellType;
        _shell_types.push_back(shellType);

        sstream2 >> nbGTF;
        _number_of_gtf.push_back(nbGTF);

        sstream2 >> factorCoeff;

        for(int i = 0; i < nbGTF; ++i)
        {
            ++nbGTFCenter;
            _factor_coefs.push_back(factorCoeff);

            getline(moldenGabFile, line);
            line = std::regex_replace(line, std::regex("D|d"), "e");

            std::stringstream sstream3(line);
            sstream3 >> exposant;
            _exposants.push_back(exposant);

            sstream3 >> coeff;
            _cgtf_coefs.push_back(coeff);
        }

        getline(moldenGabFile,line);
    } while(!trim_whitespaces(line, true, true).empty());

    _n_at_basis[centerNumber - 1] = nbGTFCenter;
}

void MOLDENGAB::read_MO_data(std::ifstream& moldenGabFile)
{
    std::string p,t;
    double a;
    std::vector<double> aa;

    if(LocaliseDataMolGab(moldenGabFile,"Spin= Beta")!=-1)
    {
        _alpha_and_beta=false;
    }

    long int pos=LocaliseDataMolGab(moldenGabFile, "[MO]");

    if(pos==-1)
    {
        std::stringstream errorMessage;
        errorMessage << "Error in MOLDENGAB::read_MO_data(): Basis (GTO) data not found." << std::endl;
        errorMessage << "Data required, please check your file.";
        print_error(errorMessage.str());
        
        std::exit(1);
    }

    moldenGabFile.seekg(pos);
    getline(moldenGabFile,p);

    do{
        for(int j=0; j<4; j++)
        {
            if(j!=0)
                getline(moldenGabFile,p);
            std::stringstream s(p);
            s>>p;

            if(p=="Ene=")
            {
                s>>a;
                _MO_energy.push_back(a);
            }
            else if(p=="Spin=")
            {
                s>>p;
                _spin_types.push_back(p);
            }
            else if(p=="Occup=")
            {
                s>>a;
                _occupation.push_back(a);
            }
        }

        for(int k=0; k<_number_of_MO_coefs; k++)
        {
            getline(moldenGabFile,p);
            std::stringstream ss(p);
            ss>>p;
            ss>>a;
            aa.push_back(a);
        }
        _MO_coefs.push_back(aa);
        aa=std::vector<double> ();
        getline(moldenGabFile,p);
        t=p;
        while(t.find(" ")!=std::string::npos)
            t.erase(t.find(" "),1);
    }while(!t.empty());

    _number_of_MO=_MO_energy.size();
}

void MOLDENGAB::PrintData()
{
    std::cout<<"Number of atoms = "<<_number_of_atoms<<std::endl;
    for(size_t i=0; i<_symbols.size(); i++)
        std::cout<<"Symbol "<<i<<" = "<<_symbols[i]<<std::endl;
    for(size_t i=0; i<_atomic_numbers.size(); i++)
        std::cout<<"Atomic number "<<i<<" = "<<_atomic_numbers[i]<<std::endl;
    for(size_t i=0; i<_coordinates.size(); i++)
        for(size_t j=0; j<_coordinates[i].size(); j++)
            std::cout<<"Coordinates "<<j<<" for atom "<<i<<" = "<<_coordinates[i][j]<<std::endl;
    for(size_t i=0; i<_num_center.size(); i++)
        std::cout<<"Num center "<<i<<" = "<<_num_center[i]<<std::endl;
    for(size_t i=0; i<_n_at_basis.size(); i++)
        std::cout<<"N at basis "<<i<<" = "<<_n_at_basis[i]<<std::endl;
    for(size_t i=0; i<_shell_types.size(); i++)
        std::cout<<"Shell type "<<i<<" = "<<_shell_types[i]<<std::endl;
    for(size_t i=0; i<_L_types.size(); i++)
        std::cout<<"L type "<<i<<" = "<<_L_types[i]<<std::endl;
    for(size_t i=0; i<_exposants.size(); i++)
        std::cout<<"Exposant "<<i<<" = "<<_exposants[i]<<std::endl;
    for(size_t i=0; i<_number_of_gtf.size(); i++)
        std::cout<<"Number of GTF "<<i<<" = "<<_number_of_gtf[i]<<std::endl;
    for(size_t i=0; i<_cgtf_coefs.size(); i++)
        std::cout<<"GTF coefficient "<<i<<" = "<<_cgtf_coefs[i]<<std::endl;
    for(size_t i=0; i<_factor_coefs.size(); i++)
        std::cout<<"CGTF coefficient "<<i<<" = "<<_factor_coefs[i]<<std::endl;
    for(size_t i=0; i<_alpha_energies.size(); i++)
        std::cout<<"Alpha MO energy "<<i<<" = "<<_MO_energy[i]<<std::endl;
    for(size_t i=0; i<_beta_energies.size(); i++)
        std::cout<<"Beta MO energy "<<i<<" = "<<_MO_energy[i]<<std::endl;
    for(size_t i=0; i<_alpha_MO_coefs.size(); i++)
        for(size_t j=0; j<_alpha_MO_coefs[i].size(); j++)
            std::cout<<"Alpha MO coefficient ["<<i<<"]["<<j<<"] = "<<_alpha_MO_coefs[i][j]<<std::endl;
    for(size_t i=0; i<_beta_MO_coefs.size(); i++)
        for(size_t j=0; j<_beta_MO_coefs[i].size(); j++)
            std::cout<<"Beta MO coefficient ["<<i<<"]["<<j<<"] = "<<_beta_MO_coefs[i][j]<<std::endl;
    for(size_t i=0; i<_alpha_occupation.size(); i++)
        std::cout<<"Alpha Occupation "<<i<<" = "<<_alpha_occupation[i]<<std::endl;
    for(size_t i=0; i<_beta_occupation.size(); i++)
        std::cout<<"Beta Occupation "<<i<<" = "<<_beta_occupation[i]<<std::endl;
    for(size_t i=0; i<_spin_types.size(); i++)
        std::cout<<"Spin types "<<i<<" = "<<_spin_types[i]<<std::endl;
    std::cout<<"Number of MO = "<<_number_of_MO<<std::endl;
    std::cout<<"Number of MO coefficients = "<<_number_of_MO_coefs<<std::endl;

    if(_alpha_and_beta)
        std::cout<<"Alpha == Beta"<<std::endl;
    else
        std::cout<<"Alpha != Beta"<<std::endl;
}

long int LocaliseDataMolGab(std::ifstream& f, std::string b)
{
    long int position;
    f.clear();
    f.seekg(0,f.beg);
    std::string test;
    bool ok=false;
    while(!f.eof())
    {    
        getline(f, test);
        if(test.find(b)!=std::string::npos)
        {
            ok=true;
            position=f.tellg();
            break;
        }
    }

    if(!ok) 
        return -1;    

    return position;
}

long int LocaliseDataMolGabBefore(std::ifstream& f, std::string b)
{
    long int position;
    f.clear();
    f.seekg(0,f.beg);
    std::string test;
    bool ok=false;
    while(!f.eof())
    {    
        position=f.tellg();
        getline(f, test);
        if(test.find(b)!=std::string::npos)
        {
            ok=true;
            break;
        }
    }

    if(!ok) 
        return -1;    

    return position;
}
