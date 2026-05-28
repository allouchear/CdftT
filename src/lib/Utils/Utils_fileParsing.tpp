#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

/*
template <typename T>
bool readOneType(std::ifstream &inputFile, const std::string &tag, T &x)
{
    std::string value;
    if (readOneString(inputFile, tag, value))
    {
        size_t i = 0;
        while (i < value.size())
        {
            if (value[i] == ',' or value[i] == ';')
            {
                value[i] = ' ';
            }
            i++;
        }
        std::stringstream ss(value);
        ss >> x;
        if (ss.fail())
        {
            return false;
        }
        return true;
    }
    return false;
}

template <typename T>
bool readListType(std::ifstream &inputFile, const std::string &tag, std::vector<T> &x)
{
    std::string value;
    if (readOneString(inputFile, tag, value))
    {
        for (char &ch : value)
        {
            if (ch == ',' or ch == ';')
            {
                ch = ' ';
            }
        }
        T a;
        std::stringstream ss(value);
        while (!ss.eof())
        {
            ss >> a;
            if (ss.fail() and !ss.eof())
            {
                return false;
            }
            x.push_back(a);
            if (ss.peek() == EOF)
            {
                break;
            }
        }
        return true;
    }
    return false;
}
*/

template <typename T>
bool readOneType(std::ifstream &inputFile, const std::string &tag, T &x)
{
    std::vector<T> values;
    bool ok = readListType(inputFile, tag, values);

    if (ok)
    {
        // Here we allow to read a list of values but only keep the first one, in case the user specified several values by mistake
        // This behaviour can be changed if so desired, just check that values contains only one element and set ok to false otherwise.

        if (!values.empty()) // values should never be empty since readListType() should return false if no value is read, but just in case
        {
            x = values[0];
        }
        else 
        {
            ok = false;
        }
    }

    return ok;
}

template <typename T, const int N>
bool readOneTypeArray(std::ifstream &inputFile, const std::string &tag, std::array<T, N> &x)
{
    std::vector<std::array<T, N>> values;
    bool ok = readListTypeArray(inputFile, tag, values);

    if (ok)
    {
        // Here we allow to read a list of values but only keep the first one, in case the user specified several values by mistake
        // This behaviour can be changed if so desired, just check that values contains only one element and set ok to false otherwise.

        if (!values.empty()) // values should never be empty since readListType() should return false if no value is read, but just in case
        {
            x = values[0];
        }
        else 
        {
            ok = false;
        }
    }

    return ok;
}

template <typename T>
bool readListType(std::ifstream &inputFile, const std::string &tag, std::vector<T> &x)
{
    bool ok = true;

    std::string value;
    if (readOneString(inputFile, tag, value))
    {
        std::regex separatorRegex("\\s*([^,;]+?)\\s*(?=[,;]|$)");

        std::sregex_iterator separatorRegexBegin(value.begin(), value.end(), separatorRegex);
        std::sregex_iterator separatorRegexEnd = std::sregex_iterator();

        if (separatorRegexBegin == separatorRegexEnd)
        {
            ok = false;
        }

        for (std::sregex_iterator it = separatorRegexBegin; ok && it != separatorRegexEnd; ++it)
        {
            std::smatch separatorRegexMatch = *it;

            T match;
            std::stringstream buffer(separatorRegexMatch[1].str());

            buffer >> match;
            if (buffer.fail() and !buffer.eof())
            {
                ok = false;
            }
            else
            {
                x.push_back(match);
            }
        }

        if (!ok)
        {
            std::stringstream errorMessage;
            errorMessage << "Error in readListType(): incorrect value or format for the \"" << tag << "\" parameter." << std::endl;
            errorMessage << "Please check documentation and the \"" << tag << "\" parameter value in the input file.";

            print_error(errorMessage.str());

            std::exit(1);
        }
    }
    else
    {
        ok = false;
    }

    return ok;
}

template<typename T, const int N>
bool readListTypeArray(std::ifstream& inputFile, const std::string& tag, std::vector<std::array<T, N>>& x)
{
    std::string value;
    bool ok = readOneString(inputFile, tag, value);

    std::cout << "Value read for tag \"" << tag << "\": " << value << std::endl;

    if (ok)
    {
        // Get the regex pattern to read read an array of N elements
        std::string tupleRegexStr = makeTupleRegex(N);

        // Check that the value matches the expected format for a list of arrays, and if so we read the arrays
        std::regex validationRegex("^\\s*" + tupleRegexStr + "(?:\\s*[,;]\\s*" + tupleRegexStr + ")*\\s*$");
        std::smatch validationRegexMatch;
        if (std::regex_search(value, validationRegexMatch, validationRegex))
        {
            std::regex tupleRegex(tupleRegexStr);

            std::sregex_iterator tupleRegexBegin(value.begin(), value.end(), tupleRegex);
            std::sregex_iterator tupleRegexEnd = std::sregex_iterator();

            if (tupleRegexBegin == tupleRegexEnd)
            {
                ok = false;
            }

            for (std::sregex_iterator it = tupleRegexBegin; ok && it != tupleRegexEnd; ++it)
            {
                std::smatch tupleRegexMatch = *it;

                std::array<T, N> array;

                for (int i = 0; ok && i < N; ++i)
                {
                    T match;
                    std::stringstream buffer(tupleRegexMatch[i + 1].str()); // i+1 because the first match (index 0) is the whole array, while the N next matches (index 1 to N included) are the elements of the array

                    buffer >> match;
                    if (buffer.fail() and !buffer.eof())
                    {
                        ok = false;
                    }
                    else
                    {
                        array[i] = match;
                    }
                }

                if (ok)
                {
                    x.push_back(array);
                }
            }
        }
        else
        {
            ok = false;
        }

        if (!ok)
        {
            std::stringstream errorMessage;
            errorMessage << "Error in readListType(): incorrect value or format for the \"" << tag << "\" parameter." << std::endl;
            errorMessage << "Please check documentation and the \"" << tag << "\" parameter value in the input file.";

            print_error(errorMessage.str());

            std::exit(1);
        }
    }

    return ok;
}
