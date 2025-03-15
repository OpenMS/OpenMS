// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DOXYGEN_WARNING_CHECKER_H
#define DOXYGEN_WARNING_CHECKER_H

#include <iostream>
#include <string>
#include <fstream>

class DoxygenWarningChecker
{
public:
    static void printUsage(const std::string& program_name)
    {
        std::cerr << "Usage:\n   " << program_name << " <path to doxygen-error.log> <doxygen version to print>\n";
    }

    static void printVersion(const std::string& version)
    {
        std::cout << "Doxygen version: " << version << std::endl;
    }

    static bool checkErrors(const std::string& file_path)
    {
        std::ifstream is(file_path);
        if (!is)
        {
            std::cerr << "Error: File '" << file_path << "' cannot be opened.\n";
            return false;
        }

        std::cout << "Opening '" << file_path << "' to check for doxygen errors...\n"
                  << "----------- ERRORS/WARNINGS -----------------" << std::endl;

        int line_count = 0, error_count = 0;
        for (std::string line; std::getline(is, line);)
        {
            if (line.empty()) continue;
            ++line_count;

            // Ignore certain warnings
            if (line.find("Consider increasing DOT_GRAPH_MAX_NODES") != std::string::npos) continue;

            std::cerr << line << '\n';
            ++error_count;
        }

        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "Skipped over " << line_count - error_count << " lines with unavoidable warnings";
        if (error_count)
        {
            std::cerr << "\n\nFound " << error_count << " Doxygen warnings. See above. Please fix them.\n";
            return false;
        }

        return true;
    }
};

#endif // DOXYGEN_WARNING_CHECKER_H
