#ifndef neuron_types_h
#define neuron_types_h

#include <map>
#include <utility>
#include <string>
#include <ostream>

namespace neuron
{

// map for translating e type string into an integer value
class e_type_map
{
public:
    e_type_map();

    // returns the code for the given string
    // -1 if the string is not in the map
    long operator[](const char *str);

    static void print(std::ostream &os);

private:
    std::map<std::string, long> tmap;
};

// map for translating m type string into an integer value
class m_type_map
{
public:
    m_type_map();

    // returns the code for the given string
    // -1 if the string is not in the map
    long operator[](const char *str);

    static void print(std::ostream &os);

private:
    std::map<std::string, long> tmap;
};

};

#endif
