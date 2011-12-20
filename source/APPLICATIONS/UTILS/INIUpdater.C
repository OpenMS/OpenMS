// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/APPLICATIONS/INIUpdater.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/VISUAL/TOPPASScene.h>


#include <QtGui/QApplication>
#include <QFileInfo>
#include <QFile>
#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_INIUpdater INIUpdater
	@brief Update INI and TOPPAS files from previous versions of OpenMS/TOPP

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ INIUpdater \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
		</tr>
	</table>
</CENTER>

	This tool can update old INI files and make them
   - compatible to new versions of %OpenMS
   - show new parameters introduced with a new %OpenMS version
   - delete old parameters which no longer have any effect

  The new INI files can be created in-place (with -i option), which will overwrite the
  existing file, but create a backup copy with [filename]_[version].ini,
  e.g.
  @code
  INIUpdater -in FileFilter.ini -i
  @endcode
  will create a file <tt>FileFilter_1.8.ini</tt> if the old ini version was 1.8.

  No backup will be created if -out is used, as the original files are not touched (unless you name them
  the same).

  

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_INIUpdater.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

using namespace OpenMS;


class TOPPINIUpdater
	: public TOPPBase
{
	public:
	TOPPINIUpdater()
		: TOPPBase("INIUpdater", "Update INI and TOPPAS files to new OpenMS version.", false),
      failed_(),
      tmp_files_()
	{
	}
  	
	protected:

	virtual void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>", StringList(), "INI/TOPPAS files that need updating.");
		setValidFormats_("in", StringList::create("ini,toppas")); 

		registerFlag_("i", "in-place: Override given INI/TOPPAS files with new content (not compatible with -out)");
		
		registerOutputFileList_("out","<files>", StringList(), "Optional list of output files (not compatible with -i).", false, false);
    setValidFormats_("out", StringList::create("ini,toppas"));
	}

  void updateTOPPAS(const String& infile,const String& outfile)
  {
    Int this_instance = getIntOption_("instance");
    INIUpdater updater;
    String tmp_ini_file = File::getTempDirectory() + "/" + File::getUniqueName() + "_INIUpdater.ini";
    tmp_files_.push_back(tmp_ini_file);

    String path = File::getExecutablePath();

    Param p;
    p.load(infile);

    // get version of TOPPAS file
    String version = "Unknown";
    if (!p.exists("info:version"))
    {
      writeLog_("No OpenMS version information found in file " + infile + "! Assuming OpenMS 1.8 and below.");
      version = "1.8.0";
    }
    else
    {
      version = p.getValue("info:version");
      // TODO: return on newer version?!
    }

    Int vertices = p.getValue("info:num_vertices");

    // update sections
    writeDebug_("#Vertices: " + vertices, 1);
    bool update_success = true;
    for (Int v=0; v<vertices; ++v)
    {
      String sec_inst = "vertices:" + String(v) + ":";
      // check for default instance
      if (!p.exists(sec_inst + "toppas_type"))
      {
        writeLog_("Update for file " + infile + " failed because the vertex #" + String(v) + " does not have a 'toppas_type' node. Check INI file for corruption!");
        update_success = false;
        break;
      }

      if (p.getValue(sec_inst + "toppas_type") != "tool")
      { // not a tool (but input/output/merge node)
        continue;
      }

      if (!p.exists(sec_inst + "tool_name"))
      {
        writeLog_("Update for file " + infile + " failed because the vertex #" + String(v) + " does not have a 'tool_name' node. Check INI file for corruption!");
        update_success = false;
        break;
      }
      
      String old_name = p.getValue(sec_inst + "tool_name");
      String new_tool;
      String ttype;
      // find mapping to new tool (might be the same name)
      if (p.exists(sec_inst + "tool_type")) ttype = p.getValue(sec_inst + "tool_type");
      if (!updater.getNewToolName(old_name, ttype, new_tool))
      {
        String type_text = ((ttype == "") ? "" : " with type '" + ttype + "' ");
        writeLog_("Update for file " + infile + " failed because the tool '" + old_name + "'" + type_text + "is unknown. TOPPAS file seems to be corrupted!");
        update_success = false;
        break;
      }

      // set new tool name
      p.setValue(sec_inst + "tool_name", new_tool);
      // delete TOPPAS type
      if (new_tool != "GenericWrapper")
      {
        p.setValue(sec_inst + "tool_type", "");
      }

      // get defaults of new tool by calling it
      int call = system(String("\"" + path + "/" + new_tool + "\"" + " -write_ini " + tmp_ini_file + " -instance " + String(this_instance)).c_str());
      if (call)
      {
        writeLog_("Update for file " + infile + " failed because the tool '" + new_tool + "' returned with an error! Check if the tool works properly.");
        update_success = false;
        break;
      }

      // update defaults with old values
      Param new_param;
      new_param.load(tmp_ini_file);
      new_param = new_param.copy(new_tool + ":1" , true);
      Param old_param = p.copy(sec_inst + "parameters", true);
      new_param.update(old_param, true, false);
      // push back changes
      p.remove(sec_inst + "parameters:");
      p.insert(sec_inst + "parameters", new_param);
    }

    if (!update_success)
    {
      failed_.push_back(infile);
      return;
    }

    p.store(tmp_ini_file);

    // update internal structure (e.g. edges format changed from 1.8 to 1.9)
    int argc = 1;
    const char* c = "IniUpdater";
    const char** argv = &c;
    QApplication app(argc, const_cast<char**>(argv), false);
    String tmp_dir = File::getTempDirectory() + "/" + File::getUniqueName();
    QDir d;
    d.mkpath(tmp_dir.toQString());
    TOPPASScene ts(0, tmp_dir.toQString(), false);
    p.store(tmp_ini_file);
    ts.load(tmp_ini_file);
    ts.store(tmp_ini_file);
    p.load(tmp_ini_file);

    // STORE
    if (outfile.empty())
    { // create a backup
      QFileInfo fi(infile.toQString());
      String new_name = String(fi.path()) + "/" + fi.completeBaseName() + "_v" + version + ".toppas";
      QFile::rename(infile.toQString(), new_name.toQString());
      // write new file
      p.store(infile);
    }
    else
    {
      p.store(outfile);
    }
  }

  void updateINI(const String& infile,const String& outfile)
  {
    Int this_instance = getIntOption_("instance");
    INIUpdater updater;
    String tmp_ini_file = File::getTempDirectory() + "/" + File::getUniqueName() + "_INIUpdater.ini";
    tmp_files_.push_back(tmp_ini_file);

    String path = File::getExecutablePath();

    Param p;
    p.load(infile);
    // get sections (usually there is only one - or the user has merged INI files manually)
    StringList sections = updater.getToolNamesFromINI(p);

    if (sections.empty())
    {
      writeLog_("Update for file " + infile + " failed because tool section does not exist. Check INI file for corruption!");
      failed_.push_back(infile);
      return;
    }

    // get version of first section
    String version = "Unknown";
    if (!p.exists(sections[0] + ":version"))
    {
      writeLog_("No OpenMS version information found in file " + infile + "! Cannot update!");
      failed_.push_back(infile);
      return;
    }
    else
    {
      version = p.getValue(sections[0] + ":version");
      // TODO: return on newer version?!
    }
    

    // update sections
    writeDebug_("Section names: " + sections.concatenate(", "), 1);
    bool update_success = true;
    for (Size s=0; s<sections.size(); ++s)
    {
      String sec_inst = sections[s] + ":" + String(this_instance) + ":";
      // check for default instance
      if (!p.exists(sec_inst + "debug"))
      {
        writeLog_("Update for file " + infile + " failed because the instance section '" + sec_inst + "' does not exist. Use -instance or check INI file for corruption!");
        update_success = false;
        break;
      }
      String new_tool;
      String ttype;
      // find mapping to new tool (might be the same name)
      if (p.exists(sec_inst + "type")) ttype = p.getValue(sec_inst + "type");
      if (!updater.getNewToolName(sections[s], ttype, new_tool))
      {
        String type_text = ((ttype == "") ? "" : " with type '" + ttype + "' ");
        writeLog_("Update for file " + infile + " failed because the tool '" + sections[s] + "'" + type_text + "is unknown. TOPPAS file seems to be corrupted!");
        update_success = false;
        break;
      }
      // get defaults of new tool by calling it
      int call = system(String("\"" + path + "/" + new_tool + "\"" + " -write_ini " + tmp_ini_file + " -instance " + String(this_instance)).c_str());
      if (call)
      {
        writeLog_("Update for file " + infile + " failed because the tool '" + new_tool + "' returned with an error! Check if the tool works properly.");
        update_success = false;
        break;
      }

      // update defaults with old values
      Param new_param;
      new_param.load(tmp_ini_file);
      new_param = new_param.copy(new_tool , true);
      Param old_param = p.copy(sections[s], true);
      new_param.update(old_param, true, false);
      // push back changes
      p.remove(sections[s] + ":");
      p.insert(new_tool   , new_param);
    }

    if (!update_success)
    {
      failed_.push_back(infile);
      return;
    }

    // STORE
    if (outfile.empty())
    { // create a backup
      QFileInfo fi(infile.toQString());
      String new_name = String(fi.path()) + "/" + fi.completeBaseName() + "_v" + version + ".ini";
      QFile::rename(infile.toQString(), new_name.toQString());
      std::cerr << "new name: " << new_name << "\n";
      // write new file
      p.store(infile);
    }
    else
    {
      p.store(outfile);
    }
  }

	ExitCodes main_(int, const char**)
	{
    StringList in  = getStringList_("in");
		StringList out = getStringList_("out");
    bool inplace = getFlag_("i");

    // consistency checks
    if (out.empty() && !inplace)
    {
      writeLog_("Cannot write output files, as neither -out nor -i are given. Use either of them, but not both!");
      return ILLEGAL_PARAMETERS;
    }
    if (out.size()>0 && inplace)
    {
      writeLog_("Two incompatible arguments given (-out and -i). Use either of them, but not both!");
      return ILLEGAL_PARAMETERS;
    }

    if (!inplace && out.size() != in.size())
    {
      writeLog_("Output and input file list length must be equal!");
      return ILLEGAL_PARAMETERS;
    }
    
    // do the conversion!
    FileHandler fh;
    for (Size i=0; i<in.size(); ++i)
    {
      FileTypes::Type f_type = fh.getType(in[i]);
      if (f_type == FileTypes::INI) updateINI(in[i], inplace ? "" : out[i]);
      else if (f_type == FileTypes::TOPPAS) updateTOPPAS(in[i], inplace ? "" : out[i]);
    }
    
    for (Size i=0; i<tmp_files_.size(); ++i)
    {
      // clean up
      File::remove(tmp_files_[i]);
    }


    if (failed_.size()>0)
    {
      writeLog_("The following INI/TOPPAS files could not be updated:\n  " + failed_.concatenate("\n  "));
      return INPUT_FILE_CORRUPT;
    }

    return EXECUTION_OK;

  }

  StringList failed_; // list of failed INI/TOPPAS files

  StringList tmp_files_;

};


/// @endcond

int main(int argc, const char** argv)
{
  TOPPINIUpdater tool;
  return tool.main(argc, argv);
}

