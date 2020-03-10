// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <QString>

class QStringList;

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

      @param verbose Print failure information?
      @return Success status
    */
    static bool findR(const QString& executable = QString("Rscript"), bool verbose = true);


    /**
      @brief Run an R script with certain arguments on the command line

      The following checks are done before running the script:
         1) [optional] 'Rscript' executable is searched (see findR() -- set @p find_R to true)
         2) The script_file is searched in 'OpenMS/share/SCRIPTS' (see findScript()).
         3) The script is run as $ Rscript <path/to/script> <arg1> <arg2> ...
      
      If any of the above steps fail, an error message is printed and false is returned.

      The 'cmd_args' are passed via commandline and should be read by the R script using R' commandArgs() function.
      Usually, the args are input and output filenames.

      @param script_file Filename of the R script
      @param cmd_args Command line arguments to the script
      @param find_R Run findR()? May be skipped if runScript() is run repeatedly
      @param verbose Print status information; also passed internally to findR() and findScript().
      @return Success status

    */
    static bool runScript(const String& script_file, const QStringList& cmd_args, const QString& executable = QString("Rscript"), bool find_R = false, bool verbose = true);

  };

}



