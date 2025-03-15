// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include "Doxygen_Warning_Checker.h"

/**
 * @brief Entry point for the Doxygen warning checker.
 * @param argc Argument count (expected: 3)
 * @param argv Argument values: 
 *        - argv[1]: Path to the Doxygen error log file.
 *        - argv[2]: Doxygen version to be printed.
 * @return 0 if no critical warnings are found; 1 otherwise.
 */
int main(int argc, char** argv)
{
    // Validate command-line arguments
    if (argc != 3)
    {
        DoxygenWarningChecker::printUsage(argv[0]);
        return 1;
    }

    std::cout << "Note: Please make sure to run the 'doc' target before running this test, "
                 "so the 'doxygen-error.log' is up to date.\n";

    // Print the provided Doxygen version
    DoxygenWarningChecker::printVersion(argv[2]);

    // Analyze the Doxygen error log file
    if (!DoxygenWarningChecker::checkErrors(argv[1]))
    {
        return 1; // Return error status if warnings are found
    }

    return 0; // Return success status if no critical warnings
}
