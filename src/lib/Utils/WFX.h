#ifndef CDFTT_WFX_H_INCLUDED
#define CDFTT_WFX_H_INCLUDED

#include <fstream>
#include <string>
#include <vector>

#include <Utils/Utils.h>


    //! A MOPC class.
    /*! This class will be used in the WFX class. */
class MOPC
{
    private:
        int _MO_Number;
        std::vector<double> _Coefficients;
    public:

            //! A default constructor.
            /*! This constructor is used to set all of the parameters on 0 or "None" value. */
        MOPC();

            //! A default desctructor.
            /*! We don't use it. */
        ~MOPC(){}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The MO number. */
        int MO_Number() {return _MO_Number;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The MO coefficients. */
        std::vector<double> Coefficients() {return _Coefficients;}

            //! A normal member taking two arguments and returning a void value.
            /*! Push back one MO number and their coefficients */
        void push_back(int, std::vector<double>);

};

    //! A NCEG class.
    /*! This class will be used in the WFX class. */
class NCEG
{
    private:
        std::string _symbol;
        std::vector<double> _gradient;
    public:

            //! A default constructor.
            /*! This constructor is used to set all of the parameters on 0 or "None" value. */
        NCEG();

            //! A default desctructor.
            /*! We don't use it. */
        ~NCEG(){}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The symbol of the atom. */
        std::string symbol() {return _symbol;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The gradient of one atom. */
        std::vector<double> gradient() {return _gradient;}

            //! A normal member taking two arguments and returning a void value.
            /*! Push back the symbol of an atom and their gradient */
        void push_back(std::string, std::vector<double>);
};

    //! A AEDF class.
    /*! This class will be used in the WFX class. */
class AEDF
{
    private:
        int _Number_of_EDF_Primitives;
        std::vector<int> _EDF_Primitives_Centers;
        std::vector<int> _EDF_Primitives_Types;
        std::vector<double> _EDF_Primitives_Exponents;
        std::vector<double> _EDF_Primitives_Coefficients;
    public:

            //! A default constructor.
            /*! This constructor is used to set all of the parameters on 0 or "None" value. */
        AEDF();

            //! A default desctructor.
            /*! We don't use it. */
        ~AEDF(){}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of EDF primitives. */
        int Number_of_EDF_Primitives() {return _Number_of_EDF_Primitives;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The EDF primitives centers. */
        std::vector<int> EDF_Primitives_Centers() {return _EDF_Primitives_Centers;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The EDF primitives types. */
        std::vector<int> EDF_Primitives_Types() {return _EDF_Primitives_Types;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The EDF primitives exponents. */
        std::vector<double> EDF_Primitives_Exponents() {return _EDF_Primitives_Exponents;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The EDF primitives coefficients. */
        std::vector<double> EDF_Primitives_Coefficients() {return _EDF_Primitives_Coefficients;}

            //! A normal member taking five arguments and returning a void value.
            /*! Push back the number, the centers, the types, the expoenents and the coefficients of EDF primitives. */
        void push_back(int, std::vector<int>, std::vector<int>, std::vector<double>, std::vector<double>);
};

    //! A WFX class.
    /*! This class will be used to read and write in the wfx format. */
class WFX
{
    private:
        std::string _Title;                                                        
        std::string _Keywords;                                                    
        int _Number_of_Nuclei;                                                
        int _Number_of_Primitives;    
        int _Number_of_Occupied_Molecular_Orbital;    
        int _Number_of_Perturbations;                                    
        std::vector<std::string> _Nuclear_Names;    
        std::vector<int> _Atomic_Number;
        std::vector<double> _Nuclear_Charges;    
        std::vector<double> _Nuclear_Cartesian_Coordinates;
        double _Net_Charge;                                                    
        int _Number_of_Electrons;                                            
        int _Number_of_Alpha_Electrons;                                        
        int _Number_of_Beta_Electrons;                                        
        int _Electronic_Spin_Multiplicity;
        std::string _Model;                            
        std::vector<int> _Primitive_Centers;                                        
        std::vector<int> _Primitive_Types;
        std::vector<std::vector<int>> _Lxyz;                                            
        std::vector<double> _Primitive_Exponents;                                    
        std::vector<std::vector<double>> _Molecular_Orbital_Occupation_Numbers;                
        std::vector<std::vector<double>> _Molecular_Orbital_Energies;
        std::vector<std::string> _Molecular_Orbital_Spin_Types;                    
        std::vector<std::vector<MOPC>> _Molecular_Orbital_Primitive_Coefficients;
        double _Energy; // Energy = T + Vne + Vee + Vnn
        double _Virial_Ratio; // (-V/T)
        std::vector<NCEG> _Nuclear_Cartesian_Energy_Gradients;
        double _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei; // ,W
        double _Full_Virial_Ratio; // ,-(V-W)/T
        int _Number_of_Core_Electrons;
        AEDF _Additionnal_Electron_Density_Function; //(EDF)

        bool _alpha_and_beta;
    public:

            //! A default constructor.
            /*! This constructor is used to set all of the parameters on 0 or "None" value. */
        WFX();

            //! A constructor.
            /*! This constructor is used to set all of the parameters with the data in the file. */
        WFX(std::ifstream&);

            //! A default desctructor.
            /*! We don't use it. */
        ~WFX(){};

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The title. */
        std::string Title() {return _Title;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The keywords. */
        std::string Keywords() {return _Keywords;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of nuclei. */
        int Number_of_Nuclei() {return _Number_of_Nuclei;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of primitives. */
        int Number_of_Primitives() {return _Number_of_Primitives;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of occupied molecular orbital. */
        int Number_of_Occupied_Molecular_Orbital() {return _Number_of_Occupied_Molecular_Orbital;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of perturbations. */
        int Number_of_Perturbations() {return _Number_of_Perturbations;}    

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of nuclear names. */
        std::vector<std::string> Nuclear_Names() {return _Nuclear_Names;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of atomic number. */
        std::vector<int> Atomic_Number() {return _Atomic_Number;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of nuclear charges. */
        std::vector<double> Nuclear_Charges() {return _Nuclear_Charges;}    

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of nuclear cartesian coordinates. */
        std::vector<double> Nuclear_Cartesian_Coordinates() {return _Nuclear_Cartesian_Coordinates;}

            //! A normal member taking no arguments and returning a double value.
            /*! \return The net charge. */
        double Net_Charge() {return _Net_Charge;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of electrons. */
        int Number_of_Electrons() {return _Number_of_Electrons;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of alpha electrons. */
        int Number_of_Alpha_Electrons() {return _Number_of_Alpha_Electrons;}                            
        
            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of beta electrons. */
        int Number_of_Beta_Electrons() {return _Number_of_Beta_Electrons;}

            //! A normal member taking no arguments and returning an int value.
            /*! \return The electronic spin multiplicity. */
        int Electronic_Spin_Multiplicity() {return _Electronic_Spin_Multiplicity;}

            //! A normal member taking no arguments and returning a std::string value.
            /*! \return The model. */
        std::string Model() {return _Model;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of primitives centers. */
        std::vector<int> Primitive_Centers() {return _Primitive_Centers;}

            //! A normal member taking no arguments and returning a std::vector<int> value.
            /*! \return The table of primitives types. */
        std::vector<int> Primitive_Types() {return _Primitive_Types;}

            //! A normal member taking no arguments and returning a std::vector<std::vector<int>> value.
            /*! \return The table of Lx, Ly, and Lz values for each primitives. */
        std::vector<std::vector<int>> Lxyz() {return _Lxyz;}    

            //! A normal member taking one argument and returning a std::vector<int> value.
            /*! \return The table of Lx, Ly, and Lz values for one primitive. */
        std::vector<int> Lxyz(int i) {return _Lxyz[i];}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of primitive exponents. */
        std::vector<double> Primitive_Exponents() {return _Primitive_Exponents;}            

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of molecular orbital occupation numbers. */
        std::vector<std::vector<double>> Molecular_Orbital_Occupation_Numbers() {return _Molecular_Orbital_Occupation_Numbers;}

            //! A normal member taking no arguments and returning a std::vector<double> value.
            /*! \return The table of molecular orbital energies. */
        std::vector<std::vector<double>> Molecular_Orbital_Energies() {return _Molecular_Orbital_Energies;}

            //! A normal member taking no arguments and returning a std::vector<std::string> value.
            /*! \return The table of molecular orbital spin types. */
        std::vector<std::string> Molecular_Orbital_Spin_Types() {return _Molecular_Orbital_Spin_Types;}

            //! A normal member taking no arguments and returning a std::vector<MOPC> value.
            /*! \return The table of molecular orbital primitive coefficients. */
        std::vector<std::vector<MOPC>> Molecular_Orbital_Primitive_Coefficients() {return _Molecular_Orbital_Primitive_Coefficients;}

            //! A normal member taking no arguments and returning a double value.
            /*! \return The energy. */
        double Energy() {return _Energy;} // Energy = T + Vne + Vee + Vnn

            //! A normal member taking no arguments and returning a double value.
            /*! \return The viriral ratio. */
        double Virial_Ratio() {return _Virial_Ratio;} // (-V/T)

            //! A normal member taking no arguments and returning a std::vector<NCEG> value.
            /*! \return The table of nuclear cartesian energy gradients. */
        std::vector<NCEG> Nuclear_Cartesian_Energy_Gradients() {return _Nuclear_Cartesian_Energy_Gradients;}

            //! A normal member taking no arguments and returning a double value.
            /*! \return The nuclear virial of energy gradient based forces on nuclei. */
        double Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei() {return _Nuclear_Virial_of_Energy_Gradient_Based_Forces_on_Nuclei;} // ,W

            //! A normal member taking no arguments and returning a double value.
            /*! \return The full virial ratio. */
        double Full_Virial_Ratio() {return _Full_Virial_Ratio;} // ,-(V-W)/T

            //! A normal member taking no arguments and returning an int value.
            /*! \return The number of core electrons. */
        int Number_of_Core_Electrons() {return _Number_of_Core_Electrons;}

            //! A normal member taking no arguments and returning an AEDF value.
            /*! \return The additionnal electron density function. */
        AEDF Additionnal_Electron_Density_Function() {return _Additionnal_Electron_Density_Function;} //(EDF)

            //! A normal member taking no arguments and returning a boolean value.
            /*! \return The boolean of alpha and beta (true if alpha and beta have the same molecular orbitals coefficients). */
        bool AlphaAndBeta() {return _alpha_and_beta;}


            //! A normal member taking four arguments and returning a std::vector<int> value.
            /*! \return The one block of int read. */
        std::vector<int> read_one_block_int(std::ifstream&, std::string, bool, int);

            //! A normal member taking three arguments and returning a std::vector<double> value.
            /*! \return The one block of real read. */
        std::vector<double> read_one_block_real(std::ifstream&, std::string, bool);

            //! A normal member taking four arguments and returning a std::vector<double> value.
            /*! \return The one block of real read. */
        std::vector<double> read_one_block_real(std::ifstream&, std::string, bool, int);

            //! A normal member taking three arguments and returning a std::vector<std::string> value.
            /*! \return The one block of std::string read. */
        std::vector<std::string> read_one_block_string(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning an int value.
            /*! \return One int read. */
        int read_int(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning a double value.
            /*! \return One real read. */
        double read_real(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning a std::string value.
            /*! \return One std::string read. */
        std::string read_string(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning a std::vector<MOPC> value.
            /*! \return The MOPC block read. */
        std::vector<MOPC> read_MOPC_block(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning a std::vector<NCEG> value.
            /*! \return The NCEG block read. */
        std::vector<NCEG> read_NCEG_block(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning an AEDF value.
            /*! \return The AEDF block read. */
        AEDF read_AEDF_block(std::ifstream&, std::string, bool);

            //! A normal member taking three arguments and returning a void value.
            /*! Read the file and set parameters on it. */
        void read_file_wfx(std::ifstream&);

            //! A normal member taking five arguments and returning a void value.
            /*! Write one block of int */
        void write_one_block_int(std::ofstream&, std::vector<int>, std::string, bool, int);

            //! A normal member taking four arguments and returning a void value.
            /*! Write one block of real */
        void write_one_block_real(std::ofstream&, std::vector<double>, std::string, bool);

            //! A normal member taking five arguments and returning a void value.
            /*! Write one block of real */
        void write_one_block_real(std::ofstream&, std::vector<double>, std::string, bool, int);

            //! A normal member taking four arguments and returning a void value.
            /*! Write one block of std::string */
        void write_one_block_string(std::ofstream&, std::vector<std::string>, std::string, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write one matrix of std::string */
        void write_one_matrix_real(std::ofstream&, std::vector<std::vector<double>>, std::string, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write one int */
        void write_int(std::ofstream&, int, std::string, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write one real */
        void write_real(std::ofstream&, double, std::string, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write one std::string */
        void write_string(std::ofstream&, std::string ,std::string, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write MOPC block */
        void write_MOPC_block(std::ofstream&, std::vector<std::vector<MOPC>>, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write NCEG block */
        void write_NCEG_block(std::ofstream&, std::vector<NCEG>, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write AEDF block */
        void write_AEDF_block(std::ofstream&, AEDF, bool);

            //! A normal member taking four arguments and returning a void value.
            /*! Write all of the parameters in a wfx file */
        void write_file_wfx(std::ofstream&);
        void PrintData() { print_error("WFX::PrintData() is not implemented yet."); }
};

    //! A function taking three arguments and returning a long int value.
    /*! \return The position of a block and the number of elements in. */
long int LocaliseBlock(std::ifstream&, int&, std::string);

    //! A function taking four arguments and returning a long int value.
    /*! \return The position of a MO block, their number, and the number of elements in. */
long int LocaliseMO(std::ifstream&, int&, int&, std::string);

    //! A function taking one argument and returning a std::vector<int> value.
    /*! \return The Lx, Ly, and Lz values of one primitive type*/
std::vector<int> setLxyz(int);

#endif
