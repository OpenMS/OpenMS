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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_LIBSSH2SECUREFILETRANSFER_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_LIBSSH2SECUREFILETRANSFER_H

#include <libssh2.h>
#include <libssh2_sftp.h>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/PeakInvestigatorImplConfig.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/SYSTEM/AbstractSecureFileTransfer.h>

namespace OpenMS
{
  class PEAKINVESTIGATORIMPL_DLLAPI LibSSH2SecureFileTransfer 
      : public AbstractSecureFileTransfer, public ProgressLogger
  {
    public:
      LibSSH2SecureFileTransfer();
      LibSSH2SecureFileTransfer(String hostname, String username, String password);
      ~LibSSH2SecureFileTransfer();

      bool downloadFile(String fromFilename, String toFilename);
      bool uploadFile(String fromFilename, String toFilename);

      void setExpectedServerHash(String expected_hash) { expected_hash_ = expected_hash; }

    protected:
      bool establishSSHSession_();
      bool confirmSSHServerIdentity_();
      bool authenticateUser_();
      bool establishSFTPSession_();
      void disconnect_();

      int state_;
      int socket_;
      struct addrinfo* host_info_;

      LIBSSH2_SESSION* ssh_session_;
      LIBSSH2_SFTP* sftp_session_;

      String expected_hash_;
  };
} // OpenMS

#endif // OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_LIBSSH2SECUREFILETRANSFER_H
