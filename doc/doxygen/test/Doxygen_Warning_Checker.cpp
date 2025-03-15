// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include "Doxygen_Warning_Checker.h"

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        DoxygenWarningChecker::printUsage(argv[0]);
        return 1;
    }

    std::cout << "Note: Please make sure to run the 'doc' target before running this test, so the 'doxygen-error.log' is up to date.\n";
    
    DoxygenWarningChecker::printVersion(argv[2]);

    if (!DoxygenWarningChecker::checkErrors(argv[1]))
    {
        return 1;
    }

    return 0;
}
