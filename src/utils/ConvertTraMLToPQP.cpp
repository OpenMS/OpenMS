// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPReader.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_ConvertTraMLToPQP ConvertTraMLToPQP

  @brief Converts TraML files to OpenSWATH transition PQP files

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPConvertTraMLToPQP : public TOPPBase
{
public:

  TOPPConvertTraMLToPQP() :
  TOPPBase("ConvertTraMLToPQP", "Converts a TraML file to an OpenSWATH transition PQP file")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input TraML file");
    setValidFormats_("in", ListUtils::create<String>("TraML"));

    registerOutputFile_("out", "<file>", "", "Output OpenSWATH transition PQP file");
    setValidFormats_("out", ListUtils::create<String>("pqp"));
  }

  ExitCodes main_(int, const char **)
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    const char * tr_file = out.c_str();

    TraMLFile traml;
    TargetedExperiment targeted_exp;

    std::cout << "Reading " << in << std::endl;
    traml.load(in, targeted_exp);
    TransitionPQPReader pqp_reader = TransitionPQPReader();
    pqp_reader.setLogType(log_type_);
    pqp_reader.convertTargetedExperimentToPQP(tr_file, targeted_exp);
    std::cout << "Writing " << out << std::endl;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPConvertTraMLToPQP tool;
  return tool.main(argc, argv);
}

/// @endcond
