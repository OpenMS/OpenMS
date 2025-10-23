// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Kohlbacher $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Types.h>
#include <clocale>
#include <string>

namespace OpenMS::Internal
{
  // Helper function to initialize OpenMS locale while preserving the environment locale
  // This is important for Python bindings to not affect Python's locale settings
  const char* initializeOpenMSLocale()
  {
    // Save the current locale
    const char* current_locale = setlocale(LC_ALL, nullptr);
    std::string saved_locale;
    if (current_locale != nullptr)
    {
      // Make a copy of the current locale string since setlocale may overwrite it
      saved_locale = current_locale;
    }
    
    // Set locale to "C" for OpenMS operations and store the result
    const char* c_locale = setlocale(LC_ALL, "C");
    
    // Restore the original locale to avoid affecting the calling environment
    if (!saved_locale.empty())
    {
      setlocale(LC_ALL, saved_locale.c_str());
    }
    
    // Return "C" which OpenMS will use internally
    return c_locale;
  }
  
  const char * OpenMS_locale = initializeOpenMSLocale();  
} // OpenMS // Internal
