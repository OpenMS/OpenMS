// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;

namespace OpenMS
{
  SearchEngineBase::SearchEngineBase(const String& tool_name, const String& tool_description, bool official, const std::vector<Citation>& citations, bool toolhandler_test) :
    TOPPBase(tool_name, tool_description, official, citations, toolhandler_test)
  {
  }

  SearchEngineBase::~SearchEngineBase() {}

  
  String SearchEngineBase::getRawfileName(int ms_level) const
  {
    String inputfile_name = getStringOption_("in");
    FileHandler fh;
    auto type = fh.getType(inputfile_name);
    switch (type)
    {
      case FileTypes::MZML:
      {
        MzMLFile mzml;
        mzml.getOptions().setMSLevels({ ms_level }); // only query MS2 (or whatever ms_level is)
        const auto& centroid_info = mzml.getCentroidInfo(inputfile_name);
        const auto& lvl_info = centroid_info.find(ms_level);
        if (lvl_info == centroid_info.end())
        {
          throw Exception::FileEmpty(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: No MS" + String(ms_level) + " spectra in input file.");
        }

        if (lvl_info->second.count_profile > 0)
        {
          if (getFlag_("force"))
          {
            OPENMS_LOG_WARN << "Warning: Profile data found, but centroid MS spectra required. "
                               "Since '-force' flag is in effect, we will continue, but results are likely bogus." << std::endl;
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error: Profile data provided but centroided MS" + String(ms_level) + " spectra required. To enforce processing (unwise!) of the data enable the -force flag (results will be bogus!).");
          }
        }
        if (lvl_info->second.count_centroided == 0)
        {
          if (getFlag_("force"))
          {
            OPENMS_LOG_WARN << "Warning: No centroided MS" + String(ms_level) + " were found, but are required. "
                               "Since '-force' flag is in effect, we will continue, but results might be bogus." << std::endl;
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Error: No centroided MS" + String(ms_level) + " spectra were found, but are required. To enforce processing of the data enable the -force flag (results will likely be bogus!).");
          }
        }
        // do no check for UNKNOWN, since it does not really tell much (UNKNOWN can only occur if meta data is missing and our peak type estimation fails (which only happens for (almost) empty spectra))
      }
      case FileTypes::MGF:
        // no warning required. MGF files should be centroided by definition
        break;
      default:
        OPENMS_LOG_WARN << "Warning: make sure that MS" << ms_level << " spectra in '" << inputfile_name << "' are centroided. Otherwise the results may be undefined!";
    }

      
    return inputfile_name;
  }

  String SearchEngineBase::getDBFilename(String db) const
  {
    String db_name(db.empty() ? getStringOption_("database") : db);
    if (!File::readable(db_name))
    {
      db_name = File::findDatabase(db_name);
    }
    return db_name;
  }

}
