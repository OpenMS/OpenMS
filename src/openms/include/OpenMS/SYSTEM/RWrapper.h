// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <QtCore/QString>

#include <QtCore/qcontainerfwd.h> // for QStringList

namespace OpenMS
{

  class String;

  /**
      @brief R-Wrapper Class.

      Call R scripts from OpenMS, mainly to produce
      plots.

      @ingroup System
  */
  class OPENMS_DLLAPI RWrapper
  {
public:

 
    /**
      @brief Look for an R script in the share/OpenMS/SCRIPT folder
      
      The script filename can be an absolute filename (in which case the filename is returned unchanged),
      or a relative filename or just a filename. In the latter cases, the script will be searched
      in the OpenMS 'SCRIPTS' path (share/OpenMS/SCRIPTS) and the full filename will be returned.

      An exception will be thrown if the file cannot be found.

      @param script_file Name of the R script
      @param verbose Print error message to OPENMS_LOG_ERROR upon FileNotFound
      @return Full filename with absolute path
      @throw Exception::FileNotFound
    */
    static String findScript(const String& script_file, bool verbose = true);

    /**
      @brief Check for presence of 'Rscript'.

      @param executable Name of the R interpreter
      @param verbose Print failure information?
      @return Success status
    */
    static bool findR(const QString& executable = QString("Rscript"), bool verbose = true);


    /**
      @brief Run an R script with certain arguments on the command line

      The following checks are done before running the script:
         1) [optional] 'Rscript' executable is searched (see findR() -- set @p find_R to true)
         2) The script_file is searched in 'OpenMS/share/SCRIPTS' (see findScript()).
         3) The script is run as $ Rscript &lt;path/to/script&gt; &lt;arg1&gt; &lt;arg2&gt; ...
      
      If any of the above steps fail, an error message is printed and false is returned.

      The 'cmd_args' are passed via commandline and should be read by the R script using R' commandArgs() function.
      Usually, the args are input and output filenames.

      @param script_file Filename of the R script
      @param cmd_args Command line arguments to the script
      @param executable Name of the R interpreter
      @param find_R Run findR()? May be skipped if runScript() is run repeatedly
      @param verbose Print status information; also passed internally to findR() and findScript().
      @return Success status
    */
    static bool runScript(const String& script_file, const QStringList& cmd_args, const QString& executable = QString("Rscript"), bool find_R = false, bool verbose = true);

  };

}



