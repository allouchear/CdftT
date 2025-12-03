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
#define BOHRTOANG   0.5291772083
#define ANGTOBOHR  (1.0/BOHRTOANG)
#define RADTODEG   57.29577951308232090712
#define AUTODEB  2.54158059
#define AUTOEV  27.21138469                         // TODO: remplacer par HARTREE_TO_EV
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
/* c en m/s */
#define slight (299792458.0) 
#define NAvogadro (6.02214129e23)
#define hPlank (6.6260693e-34)
#define epsilon0 (8.854187817e-12)
#define PRECISION 1e-10



namespace Constants
{
    /** @brief Conversion constant from Hartrees (atomic unit) to electron-volts (eV). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?eqRhrev). */
    const double HARTREE_TO_EV = 27.211386245981;

    /** @brief Conversion constant from electron-volts (eV) to Hartrees (atomic unit). Source: NIST (https://physics.nist.gov/cgi-bin/cuu/Value?eqRhrev). */
    const double EV_TO_HARTREE = 3.6749322175665e-2;
}





#endif // CDFTT_CONSTANTS_H
