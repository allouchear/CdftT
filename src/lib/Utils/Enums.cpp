#include <string>
#include <unordered_map>

#include <Utils/Enums.hpp>
#include <Utils/Utils.h>


//----------------------------------------------------------------------------------------------------//
// ENUM <-> STRING/CHAR MAPS DEFINITIONS
//----------------------------------------------------------------------------------------------------//

namespace EnumConversionMaps
{
    std::unordered_map<ELFMethod, std::string> elfMethod_string =
    {
        { ELFMethod::BECKE, "Becke" },
        { ELFMethod::SAVIN, "Savin" },
        { ELFMethod::UNKNOWN, "UNKNOWN" }
    };

    std::unordered_map<GridSize, std::string> gridSize_string =
    {
        { GridSize::COARSE, "Coarse" },
        { GridSize::CUSTOM, "Custom" },
        { GridSize::FINE, "Fine" },
        { GridSize::MEDIUM, "Medium" },
        { GridSize::UNKNOWN, "UNKNOWN" }
    };

    std::unordered_map<HFType, std::string> hfType_string =
    {
        { HFType::RHF, "RHF" },
        { HFType::UHF, "UHF" },
        { HFType::UNKNOWN, "UNKNOWN" }
    };

    std::unordered_map<OrbitalType, std::string> orbitalType_string =
    {
        { OrbitalType::ALL, "All" },
        { OrbitalType::CUSTOM, "Custom" },
        { OrbitalType::HOMO, "HOMO" },
        { OrbitalType::HOMO_LUMO, "HOMO-LUMO" },
        { OrbitalType::LUMO, "LUMO" },
        { OrbitalType::OCCUPIED, "Occupied" },
        { OrbitalType::VIRTUAL, "Virtual" },
        { OrbitalType::UNKNOWN, "UNKNOWN" }
    };

    std::unordered_map<PartitionMethod, std::string> partitionMethod_string =
    {
        { PartitionMethod::AIM_ON_GRID, "On-Grid" },
        { PartitionMethod::AIM_NEAR_GRID, "Near-Grid" },
        { PartitionMethod::AIM_NEAR_GRID_REFINEMENT, "Near-Grid-Refinement" },
        { PartitionMethod::B2S, "B2S" },
        { PartitionMethod::BBS, "BBS" },
        { PartitionMethod::BECKE, "Becke" },
        { PartitionMethod::FD, "FD" },
        { PartitionMethod::FMO, "FMO" },
        { PartitionMethod::VDD, "VDD" },
        { PartitionMethod::UNKNOWN, "UNKNOWN" }
    };

    std::unordered_map<RunType, std::string> runType_string =
    {
        { RunType::COMPUTE_DESCRIPTORS, "ComputeDescriptors" },
        { RunType::COMPUTE_ENERGY_WITH_POINT_CHARGES, "ComputeEnergyWithPointCharges" },
        { RunType::COMPUTE_GRID_DIFFERENCE, "ComputeGridDifference" },
        { RunType::COMPUTE_INTEGRALS, "ComputeIntegrals" },
        { RunType::COMPUTE_PARTIAL_CHARGES, "ComputePartialCharges" },
        { RunType::CONVERT_ORBITALS, "ConvertOrbitals" },
        { RunType::HELP, "Help" },
        { RunType::LAMBDA_DIAGNOSTIC, "LambdaDiagnostic" },
        { RunType::MAKE_DENSITY_CUBE, "MakeDensityCube" },
        { RunType::MAKE_ORBITALS_CUBE, "MakeOrbitalsCube" },
        { RunType::MAKE_ELF_CUBE, "MakeELFCube" },
        { RunType::UNKNOWN, "UNKNOWN" }
    };

    std::unordered_map<SpinType, char> spinType_char =
    {
        { SpinType::ALPHA, 'A' },
        { SpinType::BETA, 'B' },
        { SpinType::UNKNOWN, 'U' }
    };

    std::unordered_map<SpinType, std::string> spinType_string =
    {
        { SpinType::ALPHA, "Alpha" },
        { SpinType::BETA, "Beta" },
        { SpinType::ALPHA_BETA, "Alpha-Beta" },
        { SpinType::UNKNOWN, "UNKNOWN" }
    };
}


//----------------------------------------------------------------------------------------------------//
// ENUM <-> STRING/CHAR CONVERSION FUNCTIONS
//----------------------------------------------------------------------------------------------------//

std::string to_string(ELFMethod method)
{
    return enum_to_string(method, EnumConversionMaps::elfMethod_string);
}

ELFMethod elfMethod_from_string(const std::string& strMethod)
{
    return enum_from_string(strMethod, EnumConversionMaps::elfMethod_string, ELFMethod::UNKNOWN);
}


std::string to_string(GridSize size)
{
    return enum_to_string(size, EnumConversionMaps::gridSize_string);
}

GridSize gridSize_from_string(const std::string& strSize)
{
    return enum_from_string(strSize, EnumConversionMaps::gridSize_string, GridSize::UNKNOWN);
}


std::string to_string(HFType hfType)
{
    return enum_to_string(hfType, EnumConversionMaps::hfType_string);
}

HFType hfType_from_string(const std::string& strHfType)
{
    return enum_from_string(strHfType, EnumConversionMaps::hfType_string, HFType::UNKNOWN);
}


std::string to_string(OrbitalType orbitalType)
{
    return enum_to_string(orbitalType, EnumConversionMaps::orbitalType_string);
}

OrbitalType orbitalType_from_string(const std::string& strOrbitalType)
{
    return enum_from_string(strOrbitalType, EnumConversionMaps::orbitalType_string, OrbitalType::UNKNOWN);
}


std::string to_string(PartitionMethod method)
{
    return enum_to_string(method, EnumConversionMaps::partitionMethod_string);
}

PartitionMethod partitionMethod_from_string(const std::string& strMethod)
{
    return enum_from_string(strMethod, EnumConversionMaps::partitionMethod_string, PartitionMethod::UNKNOWN);
}


std::string to_string(RunType runType)
{
    return enum_to_string(runType, EnumConversionMaps::runType_string);
}

RunType runType_from_string(const std::string& strRunType)
{
    return enum_from_string(strRunType, EnumConversionMaps::runType_string, RunType::UNKNOWN);
}


char to_char(SpinType spinType)
{
    return enum_to_char(spinType, EnumConversionMaps::spinType_char);
}

std::string to_string(SpinType spinType)
{
    return enum_to_string(spinType, EnumConversionMaps::spinType_string);
}

SpinType spinType_from_char(const char charSpinType)
{
    return enum_from_char(charSpinType, EnumConversionMaps::spinType_char, SpinType::UNKNOWN);
}

SpinType spinType_from_string(const std::string& strSpinType)
{
    return enum_from_string(strSpinType, EnumConversionMaps::spinType_string, SpinType::UNKNOWN);
}
