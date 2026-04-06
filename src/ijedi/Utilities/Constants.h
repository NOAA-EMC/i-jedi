#pragma once

// -------------------------------------------------------------------------------------------------

#include <string>
#include <vector>

// -------------------------------------------------------------------------------------------------

namespace ijedi
{

    // Interfaces for getting constants
    double getConstant(const std::string constName);
    std::vector<std::string> getAllConstantsNames();

    // Function for accessing the constants from Fortran
    extern "C"
    {
        void getConstantF(const char constNameC[], double &constValueC);
    }
}  // namespace ijedi

// -------------------------------------------------------------------------------------------------
