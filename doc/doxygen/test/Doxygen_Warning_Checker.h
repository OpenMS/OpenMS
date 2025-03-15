// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#ifndef DOXYGEN_WARNING_CHECKER_H
#define DOXYGEN_WARNING_CHECKER_H

#include <iostream>
#include <string>
#include <fstream>

/**
 * @class DoxygenWarningChecker
 * @brief A utility class to analyze and report Doxygen warnings from log files.
 */
class DoxygenWarningChecker
{
public:
    /**
     * @brief Prints usage instructions for the program.
     * @param program_name Name of the executable.
     */
    static void printUsage(const std::string& program_name)
    {
        std::cerr << "Usage:\n   " << program_name << " <path to doxygen-error.log> <doxygen version to print>\n";
    }

    /**
     * @brief Prints the specified Doxygen version.
     * @param version The Doxygen version provided as an argument.
     */
    static void printVersion(const std::string& version)
    {
        std::cout << "Doxygen version: " << version << std::endl;
    }

    /**
     * @brief Reads and analyzes the specified Doxygen log file for errors/warnings.
     * @param file_path Path to the Doxygen log file.
     * @return Returns true if no critical warnings are found; otherwise, false.
     */
    static bool checkErrors(const std::string& file_path)
    {
        std::ifstream is(file_path);
        if (!is)
        {
            std::cerr << "Error: File '" << file_path << "' cannot be opened.\n";
            return false;
        }

        std::cout << "Opening '" << file_path << "' to check for Doxygen errors...\n"
                  << "----------- ERRORS/WARNINGS -----------------" << std::endl;

        int line_count = 0, error_count = 0;

        // Process each line of the log file
        for (std::string line; std::getline(is, line);)
        {
            if (line.empty()) continue; // Ignore empty lines
            ++line_count;

            // Ignore specific known non-critical warnings
            if (line.find("Consider increasing DOT_GRAPH_MAX_NODES") != std::string::npos) continue;

            // Print and count the warnings that need attention
            std::cerr << line << '\n';
            ++error_count;
        }

        // Summary of warnings found
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "Skipped over " << line_count - error_count << " lines with unavoidable warnings.\n";

        if (error_count)
        {
            std::cerr << "\n\nFound " << error_count << " Doxygen warnings. See above. Please fix them.\n";
            return false;
        }

        return true;
    }
};

#endif // DOXYGEN_WARNING_CHECKER_H
