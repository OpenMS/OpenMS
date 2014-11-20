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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_ABSTRACTSECUREFILETRANSFER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_ABSTRACTSECUREFILETRANSFER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/PeakInvestigatorImplConfig.h>

namespace OpenMS
{
  /** @brief Abstract base class for all classes that peform file transfers using a SSH-based
   * method such as SFTP or SCP.
   *
   * It should be subclassed for seperate implmentations such as those that make system
   * calls (e.g. @ref PSCPSecureFileTransfer) or use other libraries
   * (e.g. LibSSH2SecureFileTransfer).
   *
   * The downloadFile() and uploadFile() functions must be implemented in subclasses.
   */
  class PEAKINVESTIGATORIMPL_DLLAPI AbstractSecureFileTransfer
  {
    public:
      /// Constructor
      AbstractSecureFileTransfer();

      /// Constructor that sets hostname, username, and password
      AbstractSecureFileTransfer(String hostname, String username, String password);

      /// Destructor
      virtual ~AbstractSecureFileTransfer();

      String getHostname() { return hostname_; }
      void setHostname(String hostname) { hostname_ = hostname; }

      String getUsername() { return username_; }
      void setUsername(String username) { username_ = username; }

      String getPassword() { return password_; }
      void setPassword(String password) { password_ = password; }

      /** @brief Stub fuction for setting expected hash.
       *
       * This function should be overridden in any subclass that does not
       * make an external call to a system SFTP/SCP program, and performs its own
       * host authentication step.
       *
       */
      void setExpectedServerHash(String expected_hash);

      /// Function for downloading file that must be implemented in base class.
      virtual bool downloadFile(String fromFilename, String toFileName) = 0;

      /// Function for uploading file that must be implemented in base class.
      virtual bool uploadFile(String fromFilename, String toFileName) = 0;

    protected:
      String hostname_; ///< @brief SFTP server hostname
      String username_; ///< @brief SFTP username
      String password_; ///< @brief SFTP password
  };


} //OpenMS

#endif // OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_ABSTRACTSECUREFILETRANSFER_H
