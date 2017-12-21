// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDRipper.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QDir>

using std::vector;
using std::pair;
using std::map;

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_IDRipper IDRipper

  @brief IDRipper splits the protein/peptide identifications of an idXML file into several idXML files according their annotated file origin.

  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDRipper\f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDFilter</td>
  <td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDMapper</td>
  </tr>
  </table>
  </center>

  <B>Example</B>

  <p>Assuming each peptide identification in a given idXML file is annotated with its file origin (e.g. IDRipper_test.idXML) :</p>

  @p <tt><userParam type="string" name="file_origin" value="IDMerger1_test.idXML"/></tt> or <br />
  @p <tt><userParam type="string" name="file_origin" value="IDMerger2_test.idXML"/></tt>

  <p>Obviously the file contains protein/peptide identifications of IDMerger1_test.idXML and IDMerger2_test.idXML.</p>

  <p>Calling IDRipper with an input file (here: @p -in IDRipper_test.idXML) and an output directory (via @p out or @p out_path) will
  result in two idXML files stored in the specified directory and named according to their file origin.</p>

  <p>In theory, merging files with @p IDMerger and splitting the resulting file with @p IDRipper will result in the original input files.

  <B>NOTE: The meta value "file_origin" is removed by the <tt>IDSplitter</tt>!</B>

  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_IDRipper.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_IDRipper.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDRipper :
  public TOPPBase
{
public:
  TOPPIDRipper() :
    TOPPBase("IDRipper", "Split protein/peptide identification file into several files according to annotated file origin.")
  {

  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file, in which the protein/peptide identifications must be tagged with 'file_origin'");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "The path to this file is used as the output directory.", false, false);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerStringOption_("out_path", "<file>", "", "Directory for the output files after ripping according to 'file_origin'. If 'out_path' is set, 'out' is ignored.", false, false);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String file_name = getStringOption_("in");
    String out_dir = getStringOption_("out");
    String out_dir_ = getStringOption_("out_path");
    String output_directory;

    //if neither 'out' nor 'out_dir' is set throw an exception
    if (out_dir.empty() && out_dir_.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Please specify an output directory! There are two options to do so. Use 'out' to specify the directory and basename of the resulting files, or use 'out_path' to specify a path");
    }

    QString dir_path = (!out_dir.empty() ?
                        QFileInfo(out_dir.toQString()).absolutePath() :
                        QFileInfo(out_dir_.toQString()).absolutePath());

    if (!QDir(dir_path).exists())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Specified path does not exist");
    }
    output_directory = dir_path.toStdString();


    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    IdXMLFile().load(file_name, proteins, peptides);

    //ensure protein and peptide identifications are presented, otherwise we don't have to rip anything anyhow
    if (proteins.empty() || peptides.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "idXML file has to store protein and peptide identifications!");
    }

    map<String, pair<vector<ProteinIdentification>, vector<PeptideIdentification> > > ripped;

    // rip the idXML-file into several idXML according to the annotated file origin
    IDRipper ripper;
    ripper.rip(ripped, proteins, peptides);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    map<String, pair<vector<ProteinIdentification>, vector<PeptideIdentification> > >::iterator it;
    for (it = ripped.begin(); it != ripped.end(); ++it)
    {
      QString output = output_directory.toQString();
      // create full absolute path with filename
      String out = QDir::toNativeSeparators(output.append(QString("/")).append(it->first.toQString())).toStdString();
      LOG_INFO << "Storing file: '" << out << "'." << std::endl;

      QDir dir(output_directory.toQString());
      if (!dir.exists())
      {
        if (!File::writable(output_directory))
        {
          LOG_WARN << "Warning: Cannot create folder: '" << output_directory << "'." << std::endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }
        dir.mkpath(".");
      }
      IdXMLFile().store(out, it->second.first, it->second.second);
    }
    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIDRipper tool;
  return tool.main(argc, argv);
}

/// @endcond
