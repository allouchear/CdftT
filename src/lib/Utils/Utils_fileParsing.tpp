#include <fstream>
#include <sstream>
#include <string>
#include <vector>


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
