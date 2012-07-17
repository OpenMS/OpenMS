// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Immanuel Luhn$
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

    @brief IDRipper splits the protein/peptide identification of an idXML file into several idXML files according their annotated file origin.

    <CENTER>
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
    </CENTER>
  <CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDMerger \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
    </tr>
  </table>
  </CENTER>

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_IDRipper.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_IDRipper.html

   <B>Example</B>

   <p>Assuming each peptide identification in a given idXML-file is annotated with its file origin (e.g. IDRipper_test.idXML) :</p>

   @p <UserParam type="string" name="file_origin" value="IDMerger1_test.idXML"/> or <br />
   @p <UserParam type="string" name="file_origin" value="IDMerger2_test.idXML"/>

   <p>Obviously the file contains protein/peptide identifications of IDMerger1_test.idXML and IDMerger2_test.idXML.</p>

   <p>Calling IDRipper with an input file (here: @p -in IDRipper_test.idXML) and an output directory (via @p out or @p out_path) will
   result in two idXML-files stored in the specified directory and named according their file origin.</p>

  <p>In theory, merging files with @p IDMerger and rip the resulting file with @p IDSplitter will result in the original input files.

  <B>NOTE: The meta value file origin is removed by the @p IDSplitter!!</B>

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDRipper
    : public TOPPBase
{
 public:
    TOPPIDRipper()
        : TOPPBase("IDRipper","Split protein/peptide identification file into several files according annotated file origin.")
    {

    }

 protected:

    void registerOptionsAndFlags_()
    {
        registerInputFile_("in","<file>","","IdXML-file, whereas the protein/peptide identifications must be tagged with file_origin");
        setValidFormats_("in",StringList::create("idXML"));
        registerOutputFile_("out","<file>","","The path to the file is used as the output directory.",false,false);
        setValidFormats_("out",StringList::create("idXML"));
        registerStringOption_("out_path","<file>","","Directory for the IdXML-files after ripping according file_origin tag. If out_path is set, out is ignored.",false,false);
    }

    ExitCodes main_(int , const char**)
    {
        //-------------------------------------------------------------
        // parameter handling
        //-------------------------------------------------------------

        String file_name = getStringOption_("in");
        String out_dir = getStringOption_("out");
        String out_dir_ = getStringOption_("out_path");
        String output_directory;

        //if neither 'out' nor 'out_dir' is set throw an exception
        if ( out_dir.empty() && out_dir_.empty())
        {
            throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__, "Please specify an output directory! There are two options to do so. Use 'out' to specify the directory and basename of the resulting files, or use 'out_path' to specify a path");
        }

        QString dir = (!out_dir.empty()) ?
                QFileInfo(out_dir.toQString()).absolutePath()
              : QFileInfo(out_dir_.toQString()).absolutePath();

        if ( ! QDir(dir).exists() )
          throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__, "Specified path does not exist");
        output_directory = dir.toStdString();


        //-------------------------------------------------------------
        // calculations
        //-------------------------------------------------------------

        vector<ProteinIdentification> proteins;
        vector<PeptideIdentification> peptides;
        IdXMLFile().load(file_name, proteins, peptides);

        //ensure protein and peptide identifications are presented, otherwise we don't have to rip anything anyhow
        if ( proteins.empty() || peptides.empty())
        {
          throw Exception::Precondition(__FILE__,__LINE__, __PRETTY_FUNCTION__, "idXML file has to store protein and peptide identifications!");
        }

        map<String, pair< vector<ProteinIdentification>, vector<PeptideIdentification> > > ripped;

        // rip the idXML-file into several idXML according to the annotated file origin
        IDRipper ripper;
        ripper.rip(ripped, proteins, peptides);

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------

        map<String,pair< vector<ProteinIdentification>,vector<PeptideIdentification> > >::iterator it;
        for (it = ripped.begin(); it != ripped.end(); ++it )
        {
          QString output = output_directory.toQString();
          // create full absolute path with filename
          String out = QDir::toNativeSeparators(output.append(QString("/")).append(it->first.toQString())).toStdString();
          IdXMLFile().store(out, it->second.first, it->second.second);
        }
        return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
    TOPPIDRipper tool;
    return tool.main(argc,argv);
}

/// @endcond
