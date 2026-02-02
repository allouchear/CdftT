#ifndef CDFTT_CONSTANTS_H
#define CDFTT_CONSTANTS_H

/*
static gdouble hbar = 6.62606891e-34/2/PI;// PRL 1998 NIST
static gdouble e = 1.602176462e-19;
static gdouble a0 = 0.5291772083e-10;
static gdouble me = 9.10938188e-31;
static gdouble c =299792458.0;
static gdouble OneUnderfourPIEpsilon0 = 1e-7*299792458.0*299792458.0;
static gdouble NA =6.0221415e23; //(6.0221415(10)x10^23; from NIST
static gdouble kb =1.3806505e-23; //1.380 6505(24) x 10^-23 J K-1  from NIST
*/

#define BSIZE 1024
#define RADTODEG   57.29577951308232090712
#define AUTODEB  2.54158059
#define PI   3.14159265358979323846
#define AMUTOAU 1822.88848121
#define AUTOCM1 219474.63633664 
#define AUTOKCAL 627.509544796
#define KCALTOKJ 4.18400000
/* double AUTOfs = a0^2 me hbar^-1*1e15 */
#define AUTOfs 0.0241888427177214
#define fsInAU (1/AUTOfs)
#define KbInAU 3.1668114e-6
/* double fsInAKMA = 1.0/sqrt(1e-10*1e-10*1.6605655e-27*6.022045e23/4184.0)/1e15;*/
#define fsInAKMA 0.020454828110640
#define kcalcmM1 349.75511054
/* Kb = 6.022045e23*1.380662e-23/4184.0 */
/* Kb in AKMA */
#define Kb  1.9871914e-3
/* hbar in AKMA */
/* 6.62606957e-34*1.43932636e+20*0.020454828110640*1e15/(2pi) */
#define hbar   (6.62606957e-34*1.43932636e20*0.020454828110640*1e15/2.0/PI)

#define DEGTORAD   0.017453293 


#define epsilon0 (8.854187817e-12)
#define PRECISION 1e-10



namespace Constants
{
    //----------------------------------------------------------------------------------------------------//
    // PHYSICS CONSTANTS
    //----------------------------------------------------------------------------------------------------//

    /** @brief Avogadro constant in reciprocal moles (mol⁻¹). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?na). */
    const double AVOGADRO_CONSTANT = 6.02214076e23;

    /** @brief Boltzmann constant in joules per kelvin (J/K). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?k). */
    const double BOLTZMANN_CONSTANT = 1.380649e-23;

    /** @brief Planck constant in joule seconds (J·s). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?h). */
    const double PLANCK_CONSTANT = 6.62607015e-34;

    /** @brief Speed of light in vacuum in meters per second (m/s). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?c). */
    const double SPEED_OF_LIGHT = 299792458.0;

    ///////////////
    // SHORTCUTS

    /** @brief Boltzmann constant in joules per kelvin (J/K). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?k). */
    const double KB = BOLTZMANN_CONSTANT;

    /** @brief Avogadro constant in reciprocal moles (mol⁻¹). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?na). */
    const double NA = AVOGADRO_CONSTANT;


    //----------------------------------------------------------------------------------------------------//
    // CONVERSION CONSTANTS
    //----------------------------------------------------------------------------------------------------//

    //////////////////
    // ENERGY UNITS

    /** @brief Conversion constant from Hartrees (atomic energy unit) to electron-volts (eV). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?eqRhrev). */
    const double HARTREE_TO_EV = 27.211386245981;

    /** @brief Conversion constant from electron-volts (eV) to Hartrees (atomic energy unit). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?evhr). */
    const double EV_TO_HARTREE = 3.6749322175665e-2;

    /** @brief Conversion constant from Hartrees (atomic energy unit) to Joules (J). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?hrj). */
    const double HARTREE_TO_JOULE = 4.3597447222060e-18;

    /** @brief Conversion constant from Joules (J) to Hartrees (atomic energy unit). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?jhr). */
    const double JOULE_TO_HARTREE = 2.2937122783969e17;

    //////////////////
    // LENGTH UNITS

    /** @brief Conversion constant from Bohr radius (atomic length unit) to Ångstrom (Å). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0). */
    const double BOHR_RADIUS_TO_ANGSTROM = 0.529177210544;

    /** @brief Conversion constant from Ångstrom (Å) to Bohr radius (atomic length unit). */
    const double ANGSTROM_TO_BOHR_RADIUS = 1.0 / BOHR_RADIUS_TO_ANGSTROM;

}





#endif // CDFTT_CONSTANTS_H
