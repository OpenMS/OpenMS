// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/QC/QCBase.h>
#include <cstdio>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
// We do not want this class to show up in the docu:


class TOPPQualityControl : public TOPPBase
{
public:
  TOPPQualityControl() : TOPPBase("QualityControl", "Does quality control for various input file types.", false)
  {
  }
protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in_raw","<file>",{},"MzML input", false);
    setValidFormats_("in_raw", {"mzML"});
    registerInputFileList_("in_postFDR","<file>",{},"featureXML input", false);
    setValidFormats_("in_postFDR", {"featureXML"});
    registerInputFile_("in_con","<file>","","Contaminant database input", false);
    setValidFormats_("in_con", {"fasta"});
    //possible additions:
    //"mzData,mzXML,dta,dta2d,mgf,featureXML,consensusXML,idXML,pepXML,fid,mzid,trafoXML,fasta"

  }
  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    //
    // Read input, check for same length and get that length
    QCBase::Status status;
    UInt64 number_exps(0);
    StringList in_raw = updateFileStatus_(status, number_exps, "in_raw", QCBase::Requires::RAWMZML);
    StringList in_postFDR = updateFileStatus_(status, number_exps, "in_postFDR", QCBase::Requires::POSTFDRFEAT);

    // load databases and other single file inputs
    String in_con = getStringOption_("in_con");
    FASTAFile fasta_file;
    vector<FASTAFile::FASTAEntry> contaminants;
    if (!in_con.empty())
    {
      fasta_file.load(in_con,contaminants);
      status |= QCBase::Requires::CONTAMINANTS;
    }

    // Loop through file lists
    for (Size i = 0; i < number_exps; ++i)
    {
      //-------------------------------------------------------------
      // reading input
      //-------------------------------------------------------------
      MzMLFile mzml_file;
      PeakMap exp;
      if (!in_raw.empty())
      {
        mzml_file.load(in_raw[i],exp);
      }

      FeatureXMLFile fxml_file;
      FeatureMap fmap;
      if (!in_postFDR.empty())
      {
        fxml_file.load(in_postFDR[i], fmap);
      }
      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
    }
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    return EXECUTION_OK;
  }

private:
  StringList updateFileStatus_(QCBase::Status& status, UInt64& number_exps, const String& port, const QCBase::Requires& req)
  {
    // since files are optional, leave function if non are provided by the user
    StringList files = getStringList_(port);
    if (!files.empty())
    {
      if (number_exps == 0) number_exps = files.size(); // Number of experiments is determined from first non empty file list.
      if (number_exps != files.size()) // exit if any file list has different length
      {
        cerr << port + ": invalid number of files. Expected were " << number_exps << ".\n";
        exit(ILLEGAL_PARAMETERS);
      }
      status |= req;
    }
    return files;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPQualityControl tool;
  return tool.main(argc, argv);
}