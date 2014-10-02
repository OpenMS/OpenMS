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
// $Maintainer:$
// $Author: Adam Tenderholt $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_PSCPSECUREFILETRANSFER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_PSCPSECUREFILETRANSFER_H

#include <QtCore/QProcess>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/PeakInvestigatorImplConfig.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/SYSTEM/AbstractSecureFileTransfer.h>

namespace OpenMS
{
  /** @brief Wrapper around PuTTY's scp program (Windows).
   *
   * This class is used for transferring files using the secure file transfer protocol.
   */
  class PEAKINVESTIGATORIMPL_DLLAPI PSCPSecureFileTransfer 
      : public AbstractSecureFileTransfer
  {
    public:
      /** @brief Default constructor.
       *
       * This sets the QProcess channel mode to forwarding so that output is directly
       * echod to the terminal.
       */
      PSCPSecureFileTransfer();

      /** @brief Constructor that sets SFTP server hostname, SFTP username, and SFTP password.
       *
       * This sets the QProcess channel mode to forwarding so that output is directly
       * echod to the terminal.
       */
      PSCPSecureFileTransfer(String hostname, String username, String password);

      /// @brief Deconstructor.
      ~PSCPSecureFileTransfer();

      /** @brief Download a file from a SFTP server.
       *
       * @param fromFilename  The remote filename, including path.
       * @param toFilename    The local filename, including path.
       * @returns Bool indicating success.
       */
      bool downloadFile(String fromFilename, String toFilename);

      /** @brief Upload a file to a remote SFTP server.
       *
       * @param fromFilename  The local filename, including path.
       * @param toFilename    The remote filename, including path.
       * @returns Bool indicating success.
       */
      bool uploadFile(String fromFilename, String toFileName);

    protected:
      /// Display an error returned from the QProcess call to pscp.
      void displayPSCPError_(int error);

      QProcess process_; ///< @brief Class used to make system call to pscp.
  };

} //OpenMS

#endif // OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_PSCPSECUREFILETRANSFER_H
