#include "neuron_types.h"

namespace neuron
{
namespace {
const char *m_type[] = {
    "L1_DAC", "L1_DLAC", "L1_HAC", "L1_NGC-DA",
    "L1_NGC-SA", "L1_SLAC", "L23_BP", "L23_BTC",
    "L23_ChC", "L23_DBC", "L23_LBC", "L23_MC",
    "L23_NBC", "L23_NGC", "L23_PC", "L23_SBC",
    "L4_BP", "L4_BTC", "L4_ChC", "L4_DBC",
    "L4_LBC", "L4_MC", "L4_NBC", "L4_NGC",
    "L4_PC", "L4_SBC", "L4_SP", "L4_SS",
    "L5_BP", "L5_BTC", "L5_ChC", "L5_DBC",
    "L5_LBC", "L5_MC", "L5_NBC", "L5_NGC",
    "L5_SBC", "L5_STPC", "L5_TTPC1", "L5_TTPC2",
    "L5_UTPC", "L6_BP", "L6_BPC", "L6_BTC",
    "L6_ChC", "L6_DBC", "L6_IPC", "L6_LBC",
    "L6_MC", "L6_NBC", "L6_NGC", "L6_SBC",
    "L6_TPC_L1", "L6_TPC_L4", "L6_UTPC"
    };

const char *e_type[] = {
    "bAC", "bIR", "bNAC", "bSTUT",
    "cAC", "cADpyr", "cIR", "cNAC",
    "cSTUT", "dNAC", "dSTUT"
    };
}

// --------------------------------------------------------------------------
e_type_map::e_type_map()
{
    // construct the map
    size_t n = sizeof(e_type)/sizeof(char*);
    for (size_t i = 0; i < n; ++i)
       tmap[e_type[i]] = i;
}

// --------------------------------------------------------------------------
long e_type_map::operator[](const char *str)
{
    std::map<std::string, long>::iterator it = tmap.find(str);
    if (it == tmap.end())
        return -1;
    return it->second;
}

// --------------------------------------------------------------------------
void e_type_map::print(std::ostream &os)
{
    os << "e types:" << std::endl
        << "=========================" << std::endl;
    size_t n = sizeof(e_type)/sizeof(char*);
    for (size_t i = 0; i < n; ++i)
        os << e_type[i] << " = " << i << std::endl;
    os << "=========================" << std::endl;
}

// --------------------------------------------------------------------------
m_type_map::m_type_map()
{
    // construct the map
    size_t n = sizeof(m_type)/sizeof(char*);
    for (size_t i = 0; i < n; ++i)
       tmap[m_type[i]] = i;
}

// --------------------------------------------------------------------------
long m_type_map::operator[](const char *str)
{
    std::map<std::string, long>::iterator it = tmap.find(str);
    if (it == tmap.end())
        return -1;
    return it->second;
}

// --------------------------------------------------------------------------
void m_type_map::print(std::ostream &os)
{
    os << "e types:" << std::endl
        << "=========================" << std::endl;
    size_t n = sizeof(m_type)/sizeof(char*);
    for (size_t i = 0; i < n; ++i)
        os << m_type[i] << " = " << i << std::endl;
    os << "=========================" << std::endl;
}

}
