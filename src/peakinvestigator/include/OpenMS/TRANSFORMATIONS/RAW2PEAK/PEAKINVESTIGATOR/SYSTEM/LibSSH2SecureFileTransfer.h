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
  /** @brief Wrapper around libssh2.
   *
   * This class is used for transferring files via the secure file transfer protocol.
   *
   * For each call to uploadFile() and downloadFile(), the SSH/SFTP sessions are initialized
   * before file transfers can take place. Since there are multiple steps for connection that
   * involve creating various objects, and each step can fail, this class also keeps track
   * of the state of each session and should correctly disconnect sessions and free objects before
   * returning from downloadFile() or uploadFile().
   */
  class PEAKINVESTIGATORIMPL_DLLAPI LibSSH2SecureFileTransfer
      : public AbstractSecureFileTransfer, public ProgressLogger
  {
    public:
      /** @brief Default constructor.
       *
       * Sets logging type to ProgressLogger::CMD, initializes protected variables to sane values,
       * and initializes the libssh2 library.
       */
      LibSSH2SecureFileTransfer();

      /** @brief Constructor that sets SFTP server hostname, SFTP username, and SFTP password.
       *
       * Sets logging type to ProgressLogger::CMD, initializes protected variables to sane values,
       * and initializes the libssh2 library.
       */
      LibSSH2SecureFileTransfer(String hostname, String username, String password);

      /** @brief Deconstructor.
       *
       * De-initializes the libssh2 library.
       */
      ~LibSSH2SecureFileTransfer();

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
      bool uploadFile(String fromFilename, String toFilename);

      /** @brief Sets the expected server hash to bypass warning about the host in
       * confirmSSHServerIdentity_()
       */
      void setExpectedServerHash(String expected_hash) { expected_hash_ = expected_hash; }

    protected:
      /** @brief Takes care of steps for establishing the SSH session.
       *
       * This includes a hostname lookup, creating and connecting a network socket,
       * initializing the ssh session, and performing the SSH "handshake". It also calls
       * confirmSSHServerIdentity_() and authenticateUser_().
       * @returns Bool indicating success.
       */
      bool establishSSHSession_();

      /** @brief Confirm that the SSH server hash is correct, either automatically or prompting user.
       *
       * This function compares the MD5 hash against expected_hash_. If the two hashes do not match,
       * the user is presented the hashes and asked if they wish to continue.
       * @returns Bool indicating success.
       */
      bool confirmSSHServerIdentity_();

      /** @brief Authenticate the user using the username_ and password_ instance variables.
       *
       * @returns Bool indicating success.
       */
      bool authenticateUser_();

      /** @brief Establish an SFTP session on top of an existing SSH session.
       *
       * @returns Bool indicating success.
       */
      bool establishSFTPSession_();

      /// @brief Takes care of the appropriate disconnect/free steps depending on state.
      void disconnect_();

      int state_; ///< Variable to keep track of the state.
      int socket_; ///< The socket used for connection.
      struct addrinfo* host_info_; ///< Variable used for creating connecting the socket based on hostname

      LIBSSH2_SESSION* ssh_session_;
      LIBSSH2_SFTP* sftp_session_;

      String expected_hash_;
  };
} // OpenMS

#endif // OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_SYSTEM_LIBSSH2SECUREFILETRANSFER_H
