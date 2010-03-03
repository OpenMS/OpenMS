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
// $Maintainer: Andreas Bertsch$
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <QFileInfo>
#include <QDir>
#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_GenericWrapper GenericWrapper
	
	@brief Allows generically the wrapping of external tools.
	

  This tool is solely a wrapper to call external (non-OpenMS) executables/scripts.
  The input is forwarded to the external tool, and the resulting output is
  generated as specified by '-out'.
  You must provide the command line to call the tool, which should contain the placeholders
  <b>$in</b> and <b>$out</b>, which will be substituted by the given input and output name given to this tool.

  Example:
    ProteinProphet $in $out

  Some external tools do not offer an output parameter (e.g. msConvert from the ProteoWizard suite).
  Thus you cannot specify $out in the command line.
  Instead, you can specify a set of rules via the '-output_forwarding' parameter that tell this wrapper what the name of output-file 
  will be once the external tool is finished. This is used to derive the generated filename and copy it 
  to the location specified by '-out'.<br>
  Supported commands are:<br>
  path:&lt;value&gt;       The path<br>
  prefix:&lt;value&gt;     The prefix of the filename (without path)<br>
  suffix:&lt;value&gt;     The suffix of the filename (without path)<br>
  &lt;value&gt;            Leaves the value as it is<br>

  Valid &lt;value&gt;'s:<br>
  - any string which is a valid file or path
  - any of: $in, $out or $cwd (current working directory)

  e.g. <tt>path:/home/user/somedir/somefile.txt</tt> generates <tt>/home/user/somedir/</tt><br>
       <tt>path:$in</tt> will extract the path of the input file.

  All strings of '-output_forwarding' will be evaluated and concatenated, resulting in the expected output file from the external tool.

  Example:<br>
  <tt>GenericWrapper -in /home/user/myfile.raw -out /network/converted/myfile.mzML -call "msConvert $in --mzML -o /network/tmp/" -output_forwarding "/network/tmp/" "prefix:$in" ".mzML"</tt><br>
  This tells GenericWrapper to expect an output file which has just a changed suffix named mzML. It will thus expect 
  a file named '/network/tmp/myfile.mzML' which it will move to '/network/converted/myfile.mzML'.
  
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_GenericWrapper.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

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
			registerInputFile_("in", "<file>", "", "input file ");
			//setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out", "<file>", "", "output file ");
	  	//setValidFormats_("out",StringList::create("mzML"));
			registerStringOption_("call", "<call>", "", "Command line which calls the external tool, e.g. 'ProteinProphet $in $out'");
			registerStringList_("output_forwarding", "<expression>", StringList(), "mapping which allows to bind the callee's output to the '-out' parameter, in case the outfile cannot be explicitly created (msConvert for example), e.g. 'base:$in','suffix:mzML'. The callee's output will be renamed to the '-out' param.", false);

			addEmptyLine_();
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
	
			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String call(getStringOption_("call"));
      StringList rename_rules(getStringList_("output_forwarding"));
			String logfile(getStringOption_("log"));
		
      //-------------------------------------------------------------
      // call external program
      //-------------------------------------------------------------

      String cwd = String(QDir::currentPath());

      String call_hot = call;
			writeDebug_("Original call: '" + call + "'", 1);
			call_hot.substitute("$in", in).substitute("$out", out).substitute("$cwd", cwd);
			writeDebug_("Final call: '" + call_hot + "'", 1);

      // what is the expected output?
      if (call.hasSubstring("$out") && rename_rules.size()>0)
      {
        writeLog_("Warning: The call to the external program already contains an $out parameter. Thus '-output_forwarding' should not be neccessary. Please check!");
      }
      String out_internal;
      for (Size i=0;i<rename_rules.size();++i)
      {
        if (rename_rules[i].hasPrefix("path:"))
        {
          out_internal += File::path(rename_rules[i].substr(5).substitute("$in",in).substitute("$out",out).substitute("$cwd",cwd)).ensureLastChar('/');
        }
        else if (rename_rules[i].hasPrefix("prefix:"))
        {
          out_internal += File::removeExtension(File::basename(rename_rules[i].substr(7).substitute("$in",in).substitute("$out",out).substitute("$cwd",cwd)));
        }
        else if (rename_rules[i].hasPrefix("suffix:"))
        {
          QFileInfo fi(rename_rules[i].substr(7).substitute("$in",in).substitute("$out",out).substitute("$cwd",cwd).c_str());
          out_internal += String(fi.suffix());
        }
        else // simply append
        {
          out_internal += rename_rules[i];
        }

      }
      writeDebug_("Expected internal output file: " + out_internal, 1);

			Int status(system(call_hot.c_str()));
      if (status != 0)
      {
        writeLog_("Error: External program problem! (Details can be seen in the logfile: \"" + logfile + "\")");
				writeLog_("Call was '" + call_hot + "'");
        return EXTERNAL_PROGRAM_ERROR;
      }


      // rename interal tools output file (if applicable)
      if (out_internal.trim().size()>0)
      {
        if (File::exists(out_internal))
        { // copy
          QFile qf(out_internal.toQString());
          qf.rename(out.toQString());
        }
        else
        { // file not found
          writeLog_("Expected internal program's output file '" + out_internal + "' not found, despite internal tool returning successfully. Please check your if the file is named differently and adjust '-output_forwarding'.");
  				writeLog_("Call was '" + call_hot + "'");
          return EXTERNAL_PROGRAM_ERROR;
        }

      }

			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPGenericWrapper tool;
	return tool.main(argc,argv);
}

