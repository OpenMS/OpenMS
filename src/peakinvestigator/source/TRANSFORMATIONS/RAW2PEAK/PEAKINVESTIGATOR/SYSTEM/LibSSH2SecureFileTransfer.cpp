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

#include <sys/socket.h>
#include <sys/types.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <fcntl.h>

#include <QtCore/QFileInfo>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/SYSTEM/LibSSH2SecureFileTransfer.h>

// states of connection for SFTP
#define LIBRARY_UNINITIALIZED     0
#define LIBRARY_INITIALIZED       1
#define SOCKET_CONNECTED          3
#define SSH_SESSION_INITIALIZED   7
#define SSH_SESSION_ESTABLISHED   15
#define SFTP_SESSION_ESTABLISHED  31

// state modification defines
#define INITIALIZE_LIBRARY        1
#define CONNECT_SOCKET            2
#define INITIALIZE_SSH_SESSION    4
#define ESTABLISH_SSH_SESSION     8
#define ESTABLISH_SFTP_SESSION    16

#define BUFFER_SIZE 131072

using namespace std;

namespace OpenMS
{
  LibSSH2SecureFileTransfer::LibSSH2SecureFileTransfer()
    : AbstractSecureFileTransfer()
  {
    setLogType(ProgressLogger::CMD);

    socket_ = 0;
    ssh_session_ = NULL;
    sftp_session_ = NULL;

    // initialize libssh2
    state_ = LIBRARY_UNINITIALIZED;
    if (libssh2_init(0) == 0)
    {
      state_ |= INITIALIZE_LIBRARY;
    }
  }

  LibSSH2SecureFileTransfer::LibSSH2SecureFileTransfer(String hostname, String username, String password)
    : AbstractSecureFileTransfer(hostname, username, password), ProgressLogger()
  {
    setLogType(ProgressLogger::CMD);

    socket_ = 0;
    ssh_session_ = NULL;
    sftp_session_ = NULL;

    // initialize libssh2
    state_ = LIBRARY_UNINITIALIZED;
    if (libssh2_init(0) == 0)
    {
      state_ |= INITIALIZE_LIBRARY;
    }

  }

  LibSSH2SecureFileTransfer::~LibSSH2SecureFileTransfer()
  {
    if (state_ != LIBRARY_INITIALIZED)
    {
      LOG_DEBUG << "Unwinding stack incomplete: " << state_ << endl;
    }

    libssh2_exit();
  }

  bool LibSSH2SecureFileTransfer::downloadFile(String fromFilename, String toFilename)
  {
    if (state_ != LIBRARY_INITIALIZED)
    {
      return false;
    }

    // opening path with ~/ seems broken, so remove if present
    if (fromFilename.hasPrefix("~/"))
    {
      Size size = fromFilename.length();
      fromFilename = fromFilename.suffix(size - 2);
    }

    FILE* fileObject = fopen(toFilename.c_str(), "w+");
    if (!fileObject)
    {
      LOG_ERROR << "Problem opening local file: " << toFilename << endl;
      return false;
    }

    bool retval = false;
    if(establishSSHSession_() && establishSFTPSession_())
    {
      int flags = LIBSSH2_FXF_READ;

      LIBSSH2_SFTP_HANDLE* sftp_handle = libssh2_sftp_open(sftp_session_, fromFilename.c_str(), flags, 0);
      if (sftp_handle)
      {
        char buffer[BUFFER_SIZE];
        int nwritten, length;

        do
        {
          length = libssh2_sftp_read(sftp_handle, buffer, sizeof(buffer));
          if (length == 0)
          {
            retval = true;
            break; // EOF
          }
          else if (length < 0)
          {
            LOG_ERROR << "Error while reading remote file: " << fromFilename << endl;
            break;
          }

          nwritten = fwrite(buffer, sizeof(char), length, fileObject);
          if (nwritten != length)
          {
            LOG_ERROR << "Can't write data to local file: " << toFilename << endl;
            break;
          }

        } while (true);

        fclose(fileObject);
        libssh2_sftp_close_handle(sftp_handle);
      }
      else
      {
        LOG_ERROR << "Problem creating SFTP handle for " << fromFilename << endl;
      }

    }

    disconnect_();
    return retval;
  }

  bool LibSSH2SecureFileTransfer::uploadFile(String fromFilename, String toFilename)
  {
    if (state_ != LIBRARY_INITIALIZED)
    {
      LOG_ERROR << "uploadFile: library not properly initialized.\n";
      return false;
    }

    // opening path with ~/ seems broken, so remove if present
    if (toFilename.hasPrefix("~/"))
    {
      Size size = toFilename.length();
      toFilename = toFilename.suffix(size - 2);
    }

    FILE* fileObject = fopen(fromFilename.c_str(), "r+");
    if (!fileObject)
    {
      LOG_ERROR << "Problem opening local file: " << fromFilename << endl;
      return false;
    }

    bool retval = false;
    if(establishSSHSession_() && establishSFTPSession_())
    {
      int flags = LIBSSH2_FXF_WRITE | LIBSSH2_FXF_CREAT | LIBSSH2_FXF_TRUNC;
      int mode = LIBSSH2_SFTP_S_IRUSR | LIBSSH2_SFTP_S_IWUSR | LIBSSH2_SFTP_S_IRGRP | LIBSSH2_SFTP_S_IROTH;

      LIBSSH2_SFTP_HANDLE* sftp_handle = libssh2_sftp_open(sftp_session_, toFilename.c_str(), flags, mode);
      if (sftp_handle)
      {
        QFileInfo info(fromFilename.toQString());
        qint64 size = info.size();
        qint64 progress = 0;
        char buffer[BUFFER_SIZE];
        int length, nwritten;

        startProgress(0, size, "Uploading file to " + hostname_);

        do
        {
          length = fread(buffer, sizeof(char), sizeof(buffer), fileObject);
          if (length <= 0)
          {
            break;
          }

          char *ptr = buffer;
          do
          {
            nwritten = libssh2_sftp_write(sftp_handle, ptr, length);
            if (nwritten < 0)
            {
              break;
            }

            ptr += nwritten;
            length -= nwritten;
            progress += nwritten;
            setProgress(progress);

          } while (length);

        } while(nwritten);

        if (progress == size)
        {
          retval = true;
        }

        endProgress();
        fclose(fileObject);
        libssh2_sftp_close_handle(sftp_handle);
      }
      else
      {
        LOG_ERROR << "Problem creating SFTP handle for " << toFilename << endl;
      }

    }

    disconnect_();
    return retval;

  }

  bool LibSSH2SecureFileTransfer::establishSSHSession_()
  {
    LOG_DEBUG << "Trying to establish session." << endl;

    int retval;

    // setup hints to for host lookup
    struct addrinfo hints;
    memset(&hints, 0, sizeof(hints));
    hints.ai_socktype = SOCK_STREAM;

    // get information about host (i.e. IP address) and setup socket
    retval = getaddrinfo(hostname_.c_str(), "22", &hints, &host_info_);
    socket_ = socket(host_info_->ai_family, host_info_->ai_socktype, host_info_->ai_protocol);

    // connect socket!
    retval = connect(socket_, host_info_->ai_addr, host_info_->ai_addrlen);
    if (retval != 0)
    {
      LOG_ERROR << "Unable to connect socket." << endl;
      return false;
    }
    state_ |= CONNECT_SOCKET;

    ssh_session_ = libssh2_session_init();
    if (!ssh_session_)
    {
      LOG_ERROR << "Problem with libssh2_session_init()." << endl;
      return false;
    }
    state_ |= INITIALIZE_SSH_SESSION;

    // Sounds like bocking is what is desired, but not positive. Documentation says:
    // If a read is performed on a session with no data currently available, a blocking session will
    // wait for data to arrive and return what it receives. A non-blocking session will return
    // immediately with an empty buffer. If a write is performed on a session with no room for more data,
    // a blocking session will wait for room. A non-blocking session will return immediately without
    // writing anything.
    libssh2_session_set_blocking(ssh_session_, 1);

    retval = libssh2_session_handshake(ssh_session_, socket_);
    if (retval != 0)
    {
      LOG_ERROR << "Problem establishing SSH session: " << retval << endl;
      return false;
    }
    state_ |= ESTABLISH_SSH_SESSION;

    if(!confirmSSHServerIdentity_())
    {
      return false;
    }

    if(!authenticateUser_())
    {
      return false;
    }

    return true;
  }

  bool LibSSH2SecureFileTransfer::confirmSSHServerIdentity_()
  {
    const char* fingerprint = libssh2_hostkey_hash(ssh_session_, LIBSSH2_HOSTKEY_HASH_MD5);
    int size = 16; // MD5 is 16 bytes long
    char hash[16 * 3 - 1]; // seperate pairs of bytes with colons, the first number should match size

    // clear hash with NULL, and then copy fingerprint to human readable form
    memset(&hash, 0, size * 3 - 1);
    for (int i = 0; i < size - 1; i++){
      sprintf(&(hash[i*3]), "%02X:", (unsigned char) fingerprint[i]);
    }
    sprintf(&(hash[(size - 1) * 3]), "%02X", (unsigned char) fingerprint[size - 1]);

    if (expected_hash_.toUpper() != hash)
    {
      LOG_WARN << "\n**************************************************************************\n";
      LOG_WARN << "The host key for " << hostname_ << " does not match expected.\n\n";
      LOG_WARN << "Expected: '" << expected_hash_ << "'.\n";
      LOG_WARN << "Found:    '" << hash << "'.\n";
      LOG_WARN << "\n**************************************************************************\n";
      LOG_WARN << "Do you wish to proceed anyways (yes/no)?" << endl;

      String answer;
      cin >> answer;

      if (answer != "yes")
      {
        LOG_INFO << "Exiting." << endl;
        return false;
      }

      LOG_DEBUG << "Proceeding as requested.\n";

    }

    return true;

  }

  bool LibSSH2SecureFileTransfer::authenticateUser_()
  {
    if (libssh2_userauth_password(ssh_session_, username_.c_str(), password_.c_str()))
    {
      LOG_ERROR << "Authentication by password failed." << endl;
      return false;
    }

    return true;
  }

  bool LibSSH2SecureFileTransfer::establishSFTPSession_()
  {
    sftp_session_ = libssh2_sftp_init(ssh_session_);
    if (!sftp_session_)
    {
      LOG_ERROR << "Unable to start SFTP session." << endl;
      return false;
    }

    state_ |= ESTABLISH_SFTP_SESSION;
    return true;
  }

  void LibSSH2SecureFileTransfer::disconnect_()
  {
    LOG_DEBUG << "LibSSH2SecureTransfer::disconnect_(): The current state is " << state_ << endl;
    switch(state_)
    {
    case SFTP_SESSION_ESTABLISHED:
      libssh2_sftp_shutdown(sftp_session_);
      state_ ^= ESTABLISH_SFTP_SESSION;
    case SSH_SESSION_ESTABLISHED:
      libssh2_session_disconnect(ssh_session_, "Disconnecting.");
      state_ ^= ESTABLISH_SSH_SESSION;
    case SSH_SESSION_INITIALIZED:
      libssh2_session_free(ssh_session_);
      state_ ^= INITIALIZE_SSH_SESSION;
    case SOCKET_CONNECTED:
      close(socket_);
      freeaddrinfo(host_info_);
      state_ ^= CONNECT_SOCKET;
    case LIBRARY_INITIALIZED:
      break;
    default:
      LOG_ERROR << "Problem in disconnect. State: " << state_ << endl;
    }

    return;
  }

} //OpenMS
