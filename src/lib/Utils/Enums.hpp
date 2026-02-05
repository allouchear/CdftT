#ifndef CDFTT_ENUMS_HPP_INCLUDED
#define CDFTT_ENUMS_HPP_INCLUDED

#include <string>
#include <unordered_map>


//----------------------------------------------------------------------------------------------------//
// ENUM TYPES
//----------------------------------------------------------------------------------------------------//

enum class ELFMethod
{
    BECKE,
    SAVIN,
    UNKNOWN
};

enum class GridSize
{
    COARSE,
    CUSTOM,
    FINE,
    MEDIUM,
    UNKNOWN
};

enum class OrbitalType
{
    ALL,
    CUSTOM,
    HOMO,
    HOMO_LUMO,
    LUMO,
    OCCUPIED,
    VIRTUAL,
    UNKNOWN
};

enum class PartitionMethod
{
    AIM_NEAR_GRID,
    AIM_NEAR_GRID_REFINEMENT,
    AIM_ON_GRID,
    B2S,
    BBS,
    BECKE,
    FD,
    FMO,
    VDD,
    UNKNOWN
};

enum class RunType
{
    COMPUTE_DESCRIPTORS,
    COMPUTE_ENERGY_WITH_POINT_CHARGES,
    COMPUTE_GRID_DIFFERENCE,
    COMPUTE_INTEGRALS,
    COMPUTE_PARTIAL_CHARGES,
    CONVERT_ORBITALS,
    HELP,
    LAMBDA_DIAGNOSTIC,
    MAKE_DENSITY_CUBE,
    MAKE_ORBITALS_CUBE,
    MAKE_ELF_CUBE,
    UNKNOWN
};

enum class SpinType
{
    ALPHA = 0,
    BETA = 1,
    ALPHA_BETA,
    UNKNOWN
};

//----------------------------------------------------------------------------------------------------//
// ENUM <-> STRING/CHAR MAPS DECLARATION
//----------------------------------------------------------------------------------------------------//

namespace EnumConversionMaps
{
    extern std::unordered_map<ELFMethod, std::string> elfMethod_string;
    extern std::unordered_map<GridSize, std::string> gridSize_string;
    extern std::unordered_map<OrbitalType, std::string> orbitalType_string;
    extern std::unordered_map<PartitionMethod, std::string> partitionMethod_string;
    extern std::unordered_map<RunType, std::string> runType_string;
    extern std::unordered_map<SpinType, char> spinType_char;
    extern std::unordered_map<SpinType, std::string> spinType_string;
}


//----------------------------------------------------------------------------------------------------//
// ENUM <-> STRING/CHAR CONVERSION FUNCTIONS
//----------------------------------------------------------------------------------------------------//

template<typename T> char enum_to_char(const T& enumValue, const std::unordered_map<T, char>& enumToCharMap);
template<typename T> std::string enum_to_string(const T& enumValue, const std::unordered_map<T, std::string>& enumToStringMap);
template <typename T> T enum_from_char(const char charValue, const std::unordered_map<T, char>& enumToCharMap, const T& defaultValue);
template<typename T> T enum_from_string(const std::string& strValue, const std::unordered_map<T, std::string>& enumToStringMap, const T& defaultValue);

#include <Utils/Utils_enumConversion.tpp>


std::string to_string(ELFMethod method);
ELFMethod elfMethod_from_string(const std::string& strMethod);
std::string to_string(GridSize size);
GridSize gridSize_from_string(const std::string& strSize);
std::string to_string(OrbitalType type);
OrbitalType orbitalType_from_string(const std::string& strType);
std::string to_string(PartitionMethod method);
PartitionMethod partitionMethod_from_string(const std::string& strMethod);
std::string to_string(RunType runType);
RunType runType_from_string(const std::string& strRunType);
char to_char(SpinType spinType);
std::string to_string(SpinType spinType);
SpinType spinType_from_char(const char charSpinType);
SpinType spinType_from_string(const std::string& strSpinType);


#endif // CDFTT_ENUMS_HPP_INCLUDED