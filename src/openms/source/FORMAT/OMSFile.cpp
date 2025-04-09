// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/FORMAT/OMSFileLoad.h>
#include <OpenMS/FORMAT/OMSFileStore.h>

#include <fstream>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS
{
  void OMSFile::store(const String& filename, const IdentificationData& id_data)
  {
    OpenMS::Internal::OMSFileStore helper(filename, log_type_);
    helper.store(id_data);
  }

  void OMSFile::store(const String& filename, const FeatureMap& features)
  {
    OpenMS::Internal::OMSFileStore helper(filename, log_type_);
    helper.store(features);
  }

  void OMSFile::store(const String& filename, const ConsensusMap& consensus)
  {
    OpenMS::Internal::OMSFileStore helper(filename, log_type_);
    helper.store(consensus);
  }

  void OMSFile::load(const String& filename, IdentificationData& id_data)
  {
    OpenMS::Internal::OMSFileLoad helper(filename, log_type_);
    helper.load(id_data);
  }

  void OMSFile::load(const String& filename, FeatureMap& features)
  {
    OpenMS::Internal::OMSFileLoad helper(filename, log_type_);
    helper.load(features);
  }

  void OMSFile::load(const String& filename, ConsensusMap& consensus)
  {
    OpenMS::Internal::OMSFileLoad helper(filename, log_type_);
    helper.load(consensus);
  }

  void OMSFile::exportToJSON(const String& filename_in, const String& filename_out)
  {
    OpenMS::Internal::OMSFileLoad helper(filename_in, log_type_);
    ofstream output(filename_out.c_str());
    if (output.is_open())
    {
      helper.exportToJSON(output);
    }
    else
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename_out);
    }
  }

}
