// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>
#include <OpenMS/FORMAT/DB/DBAdapter.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

/**
    @page TOPP_DBImporter DBImporter

    @brief Imports an mzML file to an %OpenMS database.

  @deprecated Deprecated in OpenMS 1.9

  [Tool is no longer supported and will be removed in OpenMS 2.0]

    Besides the file to import, only the connection data has to be given.
    The data can then be retrieved by the @ref TOPP_DBExporter.

    The @em init flag can be used to create a new %OpenMS database.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_DBImporter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_DBImporter.html
*/

// We do not want this class to show up in the docu -> cond
/// @cond TOPPCLASSES

class TOPPDBImporter :
  public TOPPBase
{
public:
  TOPPDBImporter() :
    TOPPBase("DBImporter", "Imports data to an OpenMS database.")
  {

  }

protected:
  void registerOptionsAndFlags_()
  {
    registerStringOption_("user", "<user>", "", "user/login of the DB");
    registerStringOption_("host", "<host>", "localhost", "host name of the DB server", false);
    registerStringOption_("password", "<password>", "", "password for the user");
    registerIntOption_("port", "<port>", 3306, "port the DB server is running on", false);
    registerStringOption_("db", "<name>", "", "DB name");
    registerInputFile_("in", "<file>", "", "input file ", false);
    setValidFormats_("in", StringList::create("mzML"));
    registerFlag_("init", "Deletes all tables and sets up a new OpenMS database.\n"
                          "The data of 'in' is not imported!");
  }

  ExitCodes main_(int, const char **)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //varaibles
    String db, user, password, host, in;
    Int port;

    bool init = getFlag_("init");
    if (!init)
    {
      in = getStringOption_("in");
    }

    db = getStringOption_("db");
    user = getStringOption_("user");
    password = getStringOption_("password");
    host = getStringOption_("host");
    port = getIntOption_("port");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    DBConnection con;
    con.connect(db, user, password, host, port);
    DBAdapter a(con);

    if (init)
    {
      a.createDB();
    }
    else
    {
      //load input file data
      MSExperiment<Peak1D> exp;
      MzMLFile f;
      f.setLogType(log_type_);
      f.load(in, exp);

      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FORMAT_CONVERSION));

      //store data
      a.storeExperiment(exp);

      writeLog_(String(" written file to DB (id: ") + (double)(exp.getPersistenceId()) + ")");
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPDBImporter tool;
  return tool.main(argc, argv);
}

/// @endcond

