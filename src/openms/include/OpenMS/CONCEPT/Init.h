// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Exception.h>

/**
    @brief Initialization procedures

    Part of the library needs to be initialized to be used properly, these
    procedures need to be run at startup (and preferably only once).

*/

namespace OpenMS::Internal
{
  /**
    @brief Calls xercesc::XMLPlatformUtils::Initialize() under a process-wide lock.

    Xerces counts its initialisations in a plain integer and requires callers to serialise
    Initialize() and Terminate() (its threading FAQ). OpenMS initialises Xerces once when the
    library is loaded; every parser calls this again before parsing, also from several threads
    at once (e.g. the chunks of a parallel mzML read), so the calls are serialised here.
    Throws xercesc::XMLException as Initialize() does.
  */
  OPENMS_DLLAPI void xercesInitialize();

  /// Calls xercesc::XMLPlatformUtils::Terminate() under the lock of xercesInitialize()
  OPENMS_DLLAPI void xercesTerminate();
}
