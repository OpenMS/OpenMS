// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>
#include <OpenMS/FORMAT/ToolDescriptionFile.h>
#include <QFileInfo>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <QFileInfo>
#include <QDir>
#include <QtCore/QProcess>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_GenericWrapper GenericWrapper

	@brief Allows generically the wrapping of external tools.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ GenericWrapper \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any file the external tool can read </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool reading the output format </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter (to produce pepXML) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> \f$ \longrightarrow \f$ GenericWrapper (type 'ProteinProphet') \f$ \longrightarrow \f$</td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter (protXML to idXML) </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> RAW file </td>
			<td VALIGN="middle" ROWSPAN=1> \f$ \longrightarrow \f$ GenericWrapper (type 'RAWFileConvert') \f$ \longrightarrow \f$</td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool accepting mzML </td>
		</tr>

	</table>
</CENTER>

  This tool is a wrapper to call external (non-OpenMS) executables/scripts.
  Each supported tool is represented by a certain <tt>type</tt>.
  Each type exposes certain parameters which you can set (usually at least a <tt>in</tt> and <tt>out</tt>).
  
  To obtain support for more external programs, visit the OpenMS website or (if you cannot find your tool there) ask on the OpenMS mailing list.
  
  <b>The following section is for experts only, who want to add their own external tool:</b>

  Each external tool is configured via a wrapper XML file in 'OpenMS/share/OpenMS/TOOLS/EXTERNAL'. All files have the ending .ttd (TOPP tool description).
  You can add one or more wrappers (i.e. types) per file, but we recommend one. The filename does not really matter, but it should be descriptive. 
  
  The ttd file has the following structure:
  <table>

  <tr><th>type</th><td>
  The name of the type which is added to list of valid GenericWrapper types. It should be unique, otherwise you get a fatal error. 
  </td></tr>

  <tr><th>category</th><td>
  Category for TOPPAS. 
  </td></tr>

  <tr><th>cloptions</th><td>
  Command line options (arguments) appended to the executable.
  This string might contain placeholders of the form "%&lt;i&gt;"
  where each placeholder will be substituted with a value that is determined in the 
  mappings section (see below).

  Example:
@code
<cloptions>-o "%1" --mzML "%2"</cloptions>
@endcode
  </td></tr>

  <tr><th>path</th><td>
  Path (can be relative) to the executable that is executed.
  </td></tr>

  <tr><th>mappings</th><td>
  Used to replace placeholders with input parameters.
  The mapping id corresponds to the placeholder in <tt>cloptions</tt>.
  The template used as starting string is given in <tt>cl</tt>.
  All tokens therein will be replaced and the result will be patched into the <tt>cloptions</tt> string.
  Allowed tokens are:
  <ul>
  <li>%TMP  --> The current temp directory, fetched using File::getTempDirectory()
  <li>%BASENAME[file] --> the basename of a file, e.g. c:\tmp\myfile.mzML gives 'myfile'
  <li>%%<param> --> any param registered in the ini_param section, e.g. '%%in'
  </ul>

  Example:
@code
  <mapping id="2" cl="-output_file %BASENAME[%%in].mgf -temp_dir %TMP -depth 3" />
@endcode
  </td></tr>

  <tr><th>ini_param</th><td>
  Contains part of a normal INI file with describes the parameters. Valid tags are those that are in the ParamXML scheme below 'NODE', e.g. 'ITEM'.
  Example:
@code
  <ITEM name="out" value="" type="string" description="output XML file containg regression line and confidence interval" tags="output file" />
  <ITEM name="mz_tolerance" value="1" type="float" description="Tolerance in m/z dimension" />
@endcode
  </td></tr>

  </table>



	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_GenericWrapper.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

    const String getExternalToolsPath()
    {
      return File::getOpenMSDataPath() + "/TOOLS/EXTERNAL";
    }

    std::vector<Internal::ToolDescription> tools_external;

    const QStringList getExternalToolConfigFiles()
    {
      // look at environment in addition?!
      QStringList paths;
      paths << getExternalToolsPath().toQString(); // hardcoded OpenMS path
      
      QStringList all_files;
      for (int p=0;p<paths.size();++p)
      {
        QDir dir(paths[p], "*.ttd");
        QStringList files = dir.entryList();
        for (int i=0;i<files.size();++i)
        {
          files[i] = dir.absolutePath()+QDir::separator()+files[i];
        }
        all_files << files;
      }
      //StringList list = StringList::create(getExternalToolsPath() + "/" + "msconvert.ttd");
      return all_files;
    }

    const StringList getExternalToolNames()
    {
      StringList names;
      std::set<String> unique_check;
      for (Size i=0;i<tools_external.size();++i)
      {
        if (unique_check.count(tools_external[i].type[0])>0)
        {
          LOG_ERROR << "Type '" << tools_external[i].type[0] << "' exists at least twice and is ambiguous. Please fix!\n";
          exit(TOPPBase::INCOMPATIBLE_INPUT_DATA);
        }
        unique_check.insert(tools_external[i].type[0]);
        names.push_back(tools_external[i].type[0]);
      }
      return names;
    }

    static void loadToolConfig()
    {
      QStringList files = getExternalToolConfigFiles();
      for (int i=0;i<files.size();++i)
      {
        ToolDescriptionFile tdf;
        std::vector<Internal::ToolDescription> tools;
        tdf.load(String(files[i]), tools);
        tools_external.insert(tools_external.end(),tools.begin(),tools.end()); // append
      }
    }


    void createFragment(String& fragment, const Param& param)
    {

      //std::cout << " fragment '" << fragment << "' --> '";
      // e.g.:  -input %BASENAME[%%in].mzML
      // iterate through all input params and replace with values:
      for (Param::ParamIterator it=param.begin(); it!=param.end(); ++it)
      {
        fragment.substitute("%%" + it->name, it->value);
      }
      if (fragment.hasSubstring("%%")) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid '%%' found in '" + fragment + "' after replaceing all parameters!", fragment);
      
      // %TMP% replace:
      fragment.substitute("%TMP", File::getTempDirectory());

      // %BASENAME% replace
      QRegExp rx("%BASENAME\\[(.*)\\]");
      rx.setMinimal(true);
      int pos = 0;
      QString t_tmp = fragment.toQString();
      while ((pos = rx.indexIn(t_tmp, pos+1)) != -1)
      {
        //std::cout << "match @ " << pos << "\n";
        String value = rx.cap(1); // paramname (hopefully)
        // replace in fragment:
        QFileInfo qfi(value.toQString());
        fragment.substitute("%BASENAME["+value+"]", qfi.baseName());
      }
      
      if (fragment.has('%')) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Mapping still contains a '%' after substitution! Did you use % instead of %%?", fragment);

      //std::cout << fragment << "'\n";
    }



class TOPPGenericWrapper
	: public TOPPBase
{
	public:
		TOPPGenericWrapper()
			: TOPPBase("GenericWrapper", "Allows the generic wrapping of external tools.")
		{
		}

	protected:
  
 
		void registerOptionsAndFlags_()
		{
      registerSubsection_("tool","tool specific parameters");
      registerStringOption_("type", "", "", "Which external tool configuration to load?! See '" + getExternalToolsPath() + "'.", true, false);
      setValidStrings_("type", getExternalToolNames());
		}

	  Param getSubsectionDefaults_(const String& section) const
	  {
      String type = getStringOption_("type");
      // find params for 'type'
      for (Size i=0;i<tools_external.size();++i)
      {
        //std::cout << "type: " << type << "  --- toolname: " << tools_external[i].type << "\n";
        if (type == tools_external[i].type[0]) 
        {
          //std::cout << "found: " << tools_external[i].param << "\n";
          return tools_external[i].param;
        }
      }
		  return Param();
	  }

		ExitCodes main_(int , const char**)
		{
      // find the config for the tool:
      String type = getStringOption_("type");
      Internal::ToolDescription td;
      for (Size i=0;i<tools_external.size();++i)
      {
        if (type == tools_external[i].type[0]) 
        {
          td = tools_external[i];
        }
      }

      String command_args = td.commandline;
      // check for double spaces and warn
      if (command_args.hasSubstring("  "))
      {
        LOG_WARN << "Commandline contains double spaces, which is not allowed. Condensing...\n";
        while (command_args.hasSubstring("  "))
        {
          command_args.substitute("  "," ");
        }
        LOG_WARN << "result: " << command_args << "\n";
      }

      ///// construct the command line:
      // go through mappings:
      for (std::map<Int, String>::const_iterator it = td.tr_table.mapping.begin(); it != td.tr_table.mapping.end(); ++it)
      {
        //std::cout << "mapping #" << it->first << "\n";
        String fragment = it->second;
        // fragment's placeholder evaluation:
        createFragment(fragment, this->getParam_().copy("tool:",true));

        // replace fragment in cl
        //std::cout << "replace : " << "%"+String(it->first) << " with '" << fragment << "\n";
        command_args.substitute("%"+String(it->first), fragment);
      }

      QProcess builder;
      builder.setProcessChannelMode(QProcess::MergedChannels);
      String call = td.path + " " + command_args;

      writeDebug_("call command: " + call, 1);

      builder.start(call.toQString());

      if (!builder.waitForFinished())
      {
        LOG_ERROR << ("External tool returned with non-zero exit code ("+String(builder.exitCode())+"). Aborting ...");
        LOG_ERROR << ("External tool output:\n"+ String(QString(builder.readAll())));
        return EXTERNAL_PROGRAM_ERROR;
      }
      else
      {
        LOG_ERROR << ("External tool output:\n"+ String(QString(builder.readAll())));
      }

      
      // post processing (file moving via 'file' command)
      if (td.tr_table.post_move.location != "")
      {
        // find target param:
        Param p = this->getParam_().copy("tool:",true);
        String target = td.tr_table.post_move.target;
        if (!p.exists(target)) throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Cannot find target parameter '" + target + "' being mapped from external tools output!", target);
        String source = td.tr_table.post_move.location;
        // fragment's placeholder evaluation:
        createFragment(source, p);
        // check if target already exists:
        String target_file = (String)p.getValue(target);
        if (File::exists(target_file))
        {
          if (!File::remove(target_file))
          {
            LOG_ERROR << "Cannot remove conflicting file '"+target_file+"'. Check permissions! Aborting ...";
            return CANNOT_WRITE_OUTPUT_FILE;
          }
        }
        // move to target
        writeDebug_(String("moving '") + source + "' to '" + target_file + "'", 1);
        bool move_ok = QFile::rename(source.toQString(),target_file.toQString());
        if (!move_ok)
        {
          LOG_ERROR << "Moving the target file '"+target_file+"' from '" + source + "' failed! Aborting ...";
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }
      return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
  loadToolConfig();

	TOPPGenericWrapper tool;
	return tool.main(argc,argv);
}

