// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace std;
using namespace OpenMS;

/**
    @page UTILS_NucleotideID NucleotideID

    @brief TODO

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NucleotideID.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NucleotideID.html
*/

class TOPPNucleotideID :
  public TOPPBase
{

public:
  TOPPNucleotideID() :
    TOPPBase("NucleotideID", "Tool for nucleotide chain identification.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in_fasta", "<file>", "", "Fasta file containing the nucleotide library\n");
    setValidFormats_("in_fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_db_mapping", "<file>", "", "mapping file used in AccurateMassSearch\n");
    setValidFormats_("out_db_mapping", ListUtils::create<String>("csv"));

// TODO: other AccurateMassSearch files
  }


  ExitCodes main_(int, const char**)
  {
    // TODO: CODE
    //
    //
    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPNucleotideID tool;
  return tool.main(argc, argv);
}
