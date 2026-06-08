#ifndef CDFTT_DESCRIPTORS_H_INCLUDED
#define CDFTT_DESCRIPTORS_H_INCLUDED


/*
------------------------------------------------------------------------------------------------------------------------------
 Energies (hardness, mu, w, Xi, DEmin, wk-, wk+, hardnessk-, hardnessk+, hardnessk) are given in eV
 Softnesses (S, sk-, sk+) are given in eV^-1
------------------------------------------------------------------------------------------------------------------------------
    mu-         = eHOMO
    mu+         = eLUMO
    mu          = Chemical potential = (mu+ + mu-)/2
    hardness    = Chemical hardness = (mu+  -  mu-)
    Xi          = Electronegativity = -mu
    w           = Electrophilicity index = mu^2/(2 hardness) 
    w-          = propensity to donate electron = mu-^2/(2 hardness) 
    w+          = propensity to accept electron = mu+^2/(2 hardness) 
    S           = Global softness = 1/hardness
    Qmax        = Maximal electronic charge accepted by an electrophile = -mu/hardness
    DEmin       = Energy decrease if the electrophile take Qmax = -mu^2/(2 hardness) 
    fk-         = Local Fukui electrophilic attack
    fk+         = Local Fukui nucleophilic attack
    sk-         = Local softness electrophilic attack = S fk-
    sk+         = Local softness nucleophilic attack = S fk+
    wk-         = Local philicity index of electrophilic attack = w fk-
    wk+         = Local philicity index of nucleophilic attack = w fk+
    hardnessk-  = Local hardness = mu+ fk+ - mu- fk- - (mu+- mu-)*(fk+-fk-)
    hardnessk+  = Local hardness = mu+ fk+ - mu- fk- + (mu+- mu-)*(fk+-fk-)
    hardnessk   = Local hardness = mu+ fk+ - mu- fk-
    Deltafk     = Dual descripor = (fk+ - fk-) : 
                      Deltafk > 0 <=> site favored for a nucleophilic attack
                      Deltafk < 0 <=> site favored for an electrophilic attack
------------------------------------------------------------------------------------------------------------------------------
References:
- Revisiting the definition of local hardness and hardness kernel 
C. A. Polanco-Ramrez et al
Phys. Chem. Chem. Phys., 2017, 19, 12355-12364
DOI: 10.1039/c7cp00691h
- Applications of the Conceptual Density Functional Theory 
Indices to Organic Chemistry Reactivity
Luis R. Domingo, Mar Ríos-Gutiérrez and Patricia Pérez 
Molecules 2016, 21, 748; doi:10.3390/molecules21060748
- Electrodonating and Electroaccepting Powers
José L. Gazquez, André Cedillo, and Alberto Vela
J. Phys. Chem. A 2007, 111, 1966-1970, DOI: 10.1021/jp065459f
- Introducing “UCA-FUKUI” software: reactivity-index calculations
Jesús Sánchez-Márquez et al.
J Mol Model (2014) 20:2492, DOI 10.1007/s00894-014-2492-1
- Dual descriptor and molecular electrostatic potential: 
complementary tools for the study of the coordination 
chemistry of ambiphilic ligands
F.  Guégan et al.
Phys.Chem.Chem.Phys., 2014, 16 , 15558-15569, 
DOI: 10.1039/c4cp01613k
- New Dual Descriptor for Chemical Reactivity
Ch. Morell et al.
J. Phys. Chem. A 2005, 109, 205-212, DOI: 10.1021/jp046577a
------------------------------------------------------------------------------------------------------------------------------
*/

#include <fstream>
#include <vector>

#include <Common/PeriodicTable.h>
#include <Common/Structure.h>
#include <Cube/Grid.h>
#include <Cube/GridCP.h>
#include <Utils/Enums.hpp>


class Descriptors
{    
    private:
        Structure _str;
        bool _okCharge;
        
        //----------------------------------------------------------------------------------------------------//
        // GLOBAL DESCRIPTORS
        //----------------------------------------------------------------------------------------------------//

        double _mu;
        double _mup;
        double _mum;


        double _xi;
        double _hardness;
        double _w;
        double _wp;
        double _wm;
        double _S;
        double _Qmax;
        double _DEmin;

        //----------------------------------------------------------------------------------------------------//
        // LOCAL DESCRIPTORS
        //----------------------------------------------------------------------------------------------------//

        std::vector<double> _Q0;
        std::vector<double> _Qm;
        std::vector<double> _Qp;

        std::vector<double> _fk0;
        std::vector<double> _fkm;
        std::vector<double> _fkp;

        std::vector<double> _deltafk;
        std::vector<double> _wkm;
        std::vector<double> _wkp;
        std::vector<double> _Skm;
        std::vector<double> _Skp;
        std::vector<double> _Skfrac;
        std::vector<double> _hardnessk;
        std::vector<double> _hardnesskm;
        std::vector<double> _hardnesskp;


        //----------------------------------------------------------------------------------------------------//
        // PRIVATE METHODS
        //----------------------------------------------------------------------------------------------------//

        void sortCharges(const std::vector<double>& Q1, const std::vector<double>& Q2, const std::vector<double>& Q3, std::vector<double>& E, double I, double A);


    public:
        //----------------------------------------------------------------------------------------------------//
        // CONSTRUCTORS
        //----------------------------------------------------------------------------------------------------//

        //! Default Constructor
        /*! Sets all attributes to 0 and calls reset()*/
        Descriptors();

        //! Constructor
        /*! Build an object from charge std::vectors*/
        Descriptors(const Structure& S, std::vector<double> Q0, std::vector<double> Qm, std::vector<double> Qp, double I, double A);

        //! Constructor
        /*! Construct and computes descriptors calling compute_all_from_charge()*/
        Descriptors(const Structure& S, std::vector<double> Q0, std::vector<double> Qm, std::vector<double> Qp, std::vector<double> E);

        //! Constructor
        /*! Build an object from grids*/
        /*! Aimmethod specifies the grid based AIM computation used*/
        /*! Aimmethod = 0: On grid method*/
        /*! Aimmethod = 1: Near grid method, no refinement*/
        /*! Aimmethod = 2: Near grid method, with refinement*/
        Descriptors(const Grid& AIM0, const Grid& AIMM, const Grid& AIMP, double I, double A, PartitionMethod partitionMethod);

        //! Constructor
        /*! Build an object from cube files*/
        /*! Aimmethod specifies the grid based AIM computation used*/
        /*! Aimmethod = 0: On grid method*/
        /*! Aimmethod = 1: Near grid method, no refinement*/
        /*! Aimmethod = 2: Near grid method, with refinement*/
        Descriptors(std::ifstream& file0, std::ifstream& fileM, std::ifstream& fileP, double I, double A, PartitionMethod partitionMethod);

        //! Constructor
        /*! Construct and computes descriptors calling compute_all_from_grid()*/
        Descriptors(std::ifstream& file0, std::ifstream& fileM, std::ifstream& fileP, std::vector<double> E, PartitionMethod partitionMethod);
        //! Constructor
        /*! Build an object from FCHK file*/
        Descriptors(FCHK&, const PeriodicTable&);

        //! Constructor
        /*! Build an object from LOG file*/
        Descriptors(LOG&, const PeriodicTable&);

        //! Constructor
        /*! Build an object from MOLDEN or GAB files*/
        Descriptors(MOLDENGAB&, const PeriodicTable&);

        //! Constructor
        /*! Build an object from WFX file*/
        Descriptors(WFX&, const PeriodicTable&);



        //----------------------------------------------------------------------------------------------------//
        // GETTERS
        //----------------------------------------------------------------------------------------------------//

        const std::vector<double>& get_deltafk() const;

        //----------------------------------------------------------------------------------------------------//
        // OTHER PUBLIC METHODS
        //----------------------------------------------------------------------------------------------------//


        //! reset std::vectors
        /*! sets the std::vector attributes to size 0*/
        void reset();

        //! Compute chemical descriptors
        /*! Calculates all descriptors from partial charge std::vectors and stores the data in the class attributes.*/
        /*! Q0 no charge change*/
        /*! Qm electron added*/
        /*! Qp electron removed*/
        void compute_All_From_Charge(double I, double A);

        //! Compute chemical descriptors
        /*! Calculates all descriptors from Grids and stores the data in the class attributes.*/
        /*! AIM0 no charge change*/
        /*! AIMM electron added*/
        /*! AIMP electron removed*/
        /*! Aimmethod specifies the grid based AIM computation used*/
        /*! Aimmethod = 0: On grid method*/
        /*! Aimmethod = 1: Near grid method, no refinement*/
        /*! Aimmethod = 2: Near grid method, with refinement*/
        void compute_All_From_Grid(const Grid& AIM0, const Grid& AIMM, const Grid& AIMP, double I, double A, PartitionMethod partitionMethod);

        //! Compute chemical descriptors
        /*! Calculates all descriptors from cube files and stores the data in the class attributes.*/
        /*! file0 no charge change*/
        /*! fileM electron added*/
        /*! fileP electron removed*/
        /*! Aimmethod specifies the grid based AIM computation used*/
        /*! Aimmethod = 0: On grid method*/
        /*! Aimmethod = 1: Near grid method, no refinement*/
        /*! Aimmethod = 2: Near grid method, with refinement*/
        void compute_All_From_Cube(std::ifstream& file0, std::ifstream& fileM, std::ifstream& fileP, double I, double A, PartitionMethod partitionMethod);

        //! Compute chemical descriptors
        /*! Compute all descriptors given fukui and chemical potential*/
        void compute_all();

        //! Mu, Mum, Mup
        /*! Sets the values of mu given the ionization energy and the electron affinity*/
        void set_all_mu(const double& I, const double& A);

        //! fukui
        /*! Calculates and sets the values of the fukui functions*/
        void compute_fk_From_Charge();

        //! Compute charges
        /*! Compute and return partial charges given a grid and an integration method*/
        std::vector<double> compute_Charges_From_Grid(const Grid& AIM, PartitionMethod partitionMethod);

        //! Compute charges
        /*! Compute and return partial charges given a cube file and an integration method*/
        std::vector<double> compute_Charges_From_File(std::ifstream& file, PartitionMethod partitionMethod);

        //! Compute all descriptors
        /*! Calls compute charges to calculate all descriptors from cube files. Requires 3 cube files and a std::vector containing the energies of each cube files*/
        void compute_All_From_Cube(std::ifstream& file1, std::ifstream& file2, std::ifstream& file3, std::vector<double> E, PartitionMethod partitionMethod);

        //! Compute charges
        /*! Creates a Becke grid from a cube grid and computes Becke charges*/
        std::vector<double> compute_Charges_From_Becke(const Grid& grid);

        //! Compute all
        /*! compute all descriptors from 3 std::vectors of partial charges. Sorts the std::vector by total charge*/
        void compute_All_From_Charges(const Structure& S, std::vector<double> Q1, std::vector<double> Q2, std::vector<double> Q3, std::vector<double> E);



        /********************************************************************************************/

        //! print
        /*! overload of flux operator. Prints out the data in a table*/
        friend std::ostream& operator<<(std::ostream& flux, const Descriptors&);

        //! fukui
        /*! Calculates and sets the values of the fukui functions*/
        void compute_fk();

        /**
         * @brief Sets the values of mu and fukui functions.
         * 
         * @param[in] f Table of fukui function values (f-, f+) for each atom.
         * @param[in] homoEnergy HOMO energy.
         * @param[in] lumoEnergy LUMO energy.
         */
        void set_mu_fk_data(const std::vector<std::vector<double>>& f, double homoEnergy, double lumoEnergy);

        //! fukui
        /*! Sets the values of the descriptors functions*/
        void set_mu_fk_data(const std::vector<std::vector<double>>& data);
};

#endif
