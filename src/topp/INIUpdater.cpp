// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/INIUpdater.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/VISUAL/TOPPASScene.h>


#include <QApplication>
#include <QFileInfo>
#include <QFile>
#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_INIUpdater INIUpdater
@brief Update INI and TOPPAS files from previous versions of OpenMS/TOPP

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
@verbinclude TOPP_INIUpdater.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_INIUpdater.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

using namespace OpenMS;


class TOPPINIUpdater :
  public TOPPBase
{
public:
  TOPPINIUpdater() :
    TOPPBase("INIUpdater", "Update INI and TOPPAS files to new OpenMS version."),
    failed_(),
    tmp_files_()
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "INI/TOPPAS files that need updating.");
    setValidFormats_("in", ListUtils::create<String>("ini,toppas"));

    registerFlag_("i", "in-place: Override given INI/TOPPAS files with new content (not compatible with -out)");

    registerOutputFileList_("out", "<files>", StringList(), "Optional list of output files (not compatible with -i).", false, false);
    setValidFormats_("out", ListUtils::create<String>("ini,toppas"));
  }

  void updateTOPPAS(const String& infile, const String& outfile)
  {
    Int this_instance = getIntOption_("instance");
    INIUpdater updater;
    String tmp_ini_file = File::getTempDirectory() + "/" + File::getUniqueName() + "_INIUpdater.ini";
    tmp_files_.push_back(tmp_ini_file);

    String path = File::getExecutablePath();

    ParamXMLFile paramFile;
    Param p;
    paramFile.load(infile, p);

    // get version of TOPPAS file
    String version = "Unknown";
    if (!p.exists("info:version"))
    {
      writeLogWarn_("No OpenMS version information found in file " + infile + "! Assuming OpenMS 1.8 and below.");
      version = "1.8.0";
    }
    else
    {
      version = p.getValue("info:version").toString();
      // TODO: return on newer version?!
    }

    Int vertices = p.getValue("info:num_vertices");

    // update sections
    writeDebug_("#Vertices: " + String(vertices), 1);
    bool update_success = true;
    for (Int v = 0; v < vertices; ++v)
    {
      String sec_inst = "vertices:" + String(v) + ":";
      // check for default instance
      if (!p.exists(sec_inst + "toppas_type"))
      {
        writeLogWarn_("Update for file " + infile + " failed because the vertex #" + String(v) + " does not have a 'toppas_type' node. Check INI file for corruption!");
        update_success = false;
        break;
      }

      if (p.getValue(sec_inst + "toppas_type") != "tool") // not a tool (but input/output/merge node)
      {
        continue;
      }

      if (!p.exists(sec_inst + "tool_name"))
      {
        writeLogWarn_("Update for file " + infile + " failed because the vertex #" + String(v) + " does not have a 'tool_name' node. Check INI file for corruption!");
        update_success = false;
        break;
      }

      String old_name = p.getValue(sec_inst + "tool_name").toString();
      String new_tool;
      String ttype;
      // find mapping to new tool (might be the same name)
      if (p.exists(sec_inst + "tool_type")) ttype = p.getValue(sec_inst + "tool_type").toString();
      if (!updater.getNewToolName(old_name, ttype, new_tool))
      {
        String type_text = ((ttype.empty()) ? "" : " with type '" + ttype + "' ");
        writeLogWarn_("Update for file " + infile + " failed because the tool '" + old_name + "'" + type_text + "is unknown. TOPPAS file seems to be corrupted!");
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
      QProcess pr;
      QStringList arguments;
      arguments << "-write_ini";
      arguments << tmp_ini_file.toQString();
      arguments << "-instance";
      arguments << String(this_instance).toQString();
      pr.start((path + "/" + new_tool).toQString(), arguments);
      if (!pr.waitForFinished(-1))
      {
        writeLogWarn_("Update for file " + infile + " failed because the tool '" + new_tool + "' returned with an error! Check if the tool works properly.");
        update_success = false;
        break;
      }

      // update defaults with old values
      Param new_param;
      paramFile.load(tmp_ini_file, new_param);
      new_param = new_param.copy(new_tool + ":1", true);
      Param old_param = p.copy(sec_inst + "parameters", true);
      new_param.update(old_param);
      // push back changes
      p.remove(sec_inst + "parameters:");
      p.insert(sec_inst + "parameters", new_param);
    }

    if (!update_success)
    {
      failed_.push_back(infile);
      return;
    }

    paramFile.store(tmp_ini_file, p);

    // update internal structure (e.g. edges format changed from 1.8 to 1.9)
    int argc = 1;
    const char* c = "IniUpdater";
    const char** argv = &c;

    QApplication app(argc, const_cast<char**>(argv), false);
    String tmp_dir = File::getTempDirectory() + "/" + File::getUniqueName();
    QDir d;
    d.mkpath(tmp_dir.toQString());
    TOPPASScene ts(nullptr, tmp_dir.toQString(), false);
    paramFile.store(tmp_ini_file, p);
    ts.load(tmp_ini_file);
    ts.store(tmp_ini_file);
    paramFile.load(tmp_ini_file, p);

    // STORE
    if (outfile.empty()) // create a backup
    {
      QFileInfo fi(infile.toQString());
      String new_name = String(fi.path()) + "/" + fi.completeBaseName() + "_v" + version + ".toppas";
      QFile::rename(infile.toQString(), new_name.toQString());
      // write new file
      paramFile.store(infile, p);
    }
    else
    {
      paramFile.store(outfile, p);
    }
  }

  void updateINI(const String& infile, const String& outfile)
  {
    Int this_instance = getIntOption_("instance");
    INIUpdater updater;
    String tmp_ini_file = File::getTempDirectory() + "/" + File::getUniqueName() + "_INIUpdater.ini";
    tmp_files_.push_back(tmp_ini_file);

    String path = File::getExecutablePath();

    Param p;
    ParamXMLFile paramFile;
    paramFile.load(infile, p);
    // get sections (usually there is only one - or the user has merged INI files manually)
    StringList sections = updater.getToolNamesFromINI(p);

    if (sections.empty())
    {
      writeLogWarn_("Update for file " + infile + " failed because tool section does not exist. Check INI file for corruption!");
      failed_.push_back(infile);
      return;
    }

    // get version of first section
    String version_old = "Unknown";
    if (!p.exists(sections[0] + ":version"))
    {
      writeLogWarn_("No OpenMS version information found in file " + infile + "! Cannot update!");
      failed_.push_back(infile);
      return;
    }
    else
    {
      version_old = p.getValue(sections[0] + ":version").toString();
      // TODO: return on newer version?!
    }


    // update sections
    writeDebug_("Section names: " + ListUtils::concatenate(sections, ", "), 1);
    bool update_success = true;
    for (Size s = 0; s < sections.size(); ++s)
    {
      String sec_inst = sections[s] + ":" + String(this_instance) + ":";
      // check for default instance
      if (!p.exists(sec_inst + "debug"))
      {
        writeLogWarn_("Update for file '" + infile + "' failed because the instance section '" + sec_inst + "' does not exist. Use -instance or check INI file for corruption!");
        update_success = false;
        break;
      }
      String new_tool;
      String ttype;
      // find mapping to new tool (might be the same name)
      if (p.exists(sec_inst + "type")) ttype = p.getValue(sec_inst + "type").toString();
      if (!updater.getNewToolName(sections[s], ttype, new_tool))
      {
        String type_text = ((ttype.empty()) ? "" : " with type '" + ttype + "' ");
        writeLogWarn_("Update for file '" + infile + "' failed because the tool '" + sections[s] + "'" + type_text + "is unknown. TOPPAS file seems to be corrupted!");
        update_success = false;
        break;
      }
      // get defaults of new tool by calling it
      QProcess pr;
      QStringList arguments;
      arguments << "-write_ini";
      arguments << tmp_ini_file.toQString();
      arguments << "-instance";
      arguments << String(this_instance).toQString();
      pr.start((path + "/" + new_tool).toQString(), arguments);
      if (!pr.waitForFinished(-1))
      {
        writeLogWarn_("Update for file '" + infile + "' failed because the tool '" + new_tool + "' returned with an error! Check if the tool works properly.");
        update_success = false;
        break;
      }

      // update defaults with old values
      Param new_param;
      paramFile.load(tmp_ini_file, new_param);
      new_param = new_param.copy(new_tool, true);
      Param old_param = p.copy(sections[s], true);
      new_param.update(old_param);
      // push back changes
      p.remove(sections[s] + ":");
      p.insert(new_tool, new_param);
    }

    if (!update_success)
    {
      failed_.push_back(infile);
      return;
    }

    // STORE
    if (outfile.empty()) // create a backup
    {
      QFileInfo fi(infile.toQString());
      String backup_filename = String(fi.path()) + "/" + fi.completeBaseName() + "_v" + version_old + ".ini";
      QFile::rename(infile.toQString(), backup_filename.toQString());
      std::cout << "Backup of input file created: " << backup_filename << std::endl;
      // write updated/new file
      paramFile.store(infile, p);
    }
    else
    {
      paramFile.store(outfile, p);
    }
  }

  ExitCodes main_(int, const char**) override
  {
    StringList in  = getStringList_("in");
    StringList out = getStringList_("out");
    bool inplace = getFlag_("i");

    // consistency checks
    if (out.empty() && !inplace)
    {
      writeLogError_("Cannot write output files, as neither -out nor -i are given. Use either of them, but not both!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    if (!out.empty() && inplace)
    {
      writeLogError_("Two incompatible arguments given (-out and -i). Use either of them, but not both!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!inplace && out.size() != in.size())
    {
      writeLogError_("Output and input file list length must be equal!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // do the conversion!
    FileHandler fh;
    for (Size i = 0; i < in.size(); ++i)
    {
      FileTypes::Type f_type = fh.getType(in[i]);
      if (f_type == FileTypes::INI) updateINI(in[i], inplace ? "" : out[i]);
      else if (f_type == FileTypes::TOPPAS) updateTOPPAS(in[i], inplace ? "" : out[i]);
    }

    for (Size i = 0; i < tmp_files_.size(); ++i)
    {
      // clean up
      File::remove(tmp_files_[i]);
    }


    if (!failed_.empty())
    {
      writeLogError_("The following INI/TOPPAS files could not be updated:\n  " + ListUtils::concatenate(failed_, "\n  "));
      return INPUT_FILE_CORRUPT;
    }

    return EXECUTION_OK;

  }

  StringList failed_; // list of failed INI/TOPPAS files

  StringList tmp_files_;

};

int main(int argc, const char** argv)
{
  TOPPINIUpdater tool;
  return tool.main(argc, argv);
}

/// @endcond
