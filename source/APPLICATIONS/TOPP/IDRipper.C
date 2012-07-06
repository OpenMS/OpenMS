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
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/ID/IDRipper.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QDir>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IDRipper IDRipper

    @brief Splits one idXML file into several idXML files according to file origin.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDRipper\f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMerger</td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IsoformResolver</td>
        </tr>
    </table>
    </CENTER>

    In general, the number of idXML files that can result from splitting is not limited.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_IDRipper.cli

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDRipper
    : public TOPPBase
{
 public:
    TOPPIDRipper()
        : TOPPBase("IDRipper","Split one protein/peptide identification file into several files.")
    {

    }

 protected:

    void registerOptionsAndFlags_()
    {
        registerInputFile_("in","<file>","","IdXML-file, whereas the protein/peptide identifications must be tagged with file_origin");
        setValidFormats_("in",StringList::create("idXML"));
        registerOutputFile_("out","<file>","","output directory for files",false,false);
        setValidFormats_("out",StringList::create("idXML"));
        registerStringOption_("out_path","<file>","","Directory for the IdXML-files after ripping according file_origin tag. If out_path is set, out is ignored.",false,true);
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

        // if the output-file is set use the path as directory
        if ( !out_dir.empty())
        {
          //if ( ! QDir(out_dir.toQString()).exists() )
           // throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__, "Specified path does not exist");
          QDir dir = QFileInfo(out_dir.toQString()).dir();
          output_directory = dir.dirName().toStdString();
        }
        // otherwise use specified output directory
        else
        {
          //if ( ! QDir(out_dir_.toQString()).exists() )
           // throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__, "Specified path does not exist");
          QDir dir = QFileInfo(out_dir_.toQString()).dir();
          output_directory = dir.dirName().toStdString();
        }

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
          //@TODO absolute relative path
          QString output = output_directory.toQString();
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
