#include <cctype>
#include <string>
#include <unordered_map>

#include "../Utils/Utils.h"


template<typename T>
char enum_to_char(const T& enumValue, const std::unordered_map<T, char>& enumToCharMap)
{
    char charValue = 'U';

    if (const auto& it = enumToCharMap.find(enumValue); it != enumToCharMap.end())
    {
        charValue = it->second;
    }

    return charValue;
}

template<typename T>
std::string enum_to_string(const T& enumValue, const std::unordered_map<T, std::string>& enumToStringMap)
{
    std::string strMethod = "UNKNOWN";

    if (const auto& it = enumToStringMap.find(enumValue); it != enumToStringMap.end())
    {
        strMethod = it->second;
    }

    return strMethod;
}


template <typename T>
T enum_from_char(const char charValue, const std::unordered_map<T, char>& enumToCharMap, const T& defaultValue)
{
    T enumValue = defaultValue;

    char upperCharValue = std::toupper(charValue);
    for (const auto& pair : enumToCharMap)
    {
        if (std::toupper(pair.second) == upperCharValue)
        {
            enumValue = pair.first;
            break;
        }
    }

    return enumValue;
}

template<typename T>
T enum_from_string(const std::string& strValue, const std::unordered_map<T, std::string>& enumToStringMap, const T& defaultValue)
{
    T enumValue = defaultValue;

    std::string upperStrValue = to_upper(strValue);
    for (const auto& pair : enumToStringMap)
    {
        if (to_upper(pair.second) == upperStrValue)
        {
            enumValue = pair.first;
            break;
        }
    }

    return enumValue;
}