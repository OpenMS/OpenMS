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

#include <fcntl.h> // used for SFTP transfers
#include <zlib.h> // used for g'zipping tar files

#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakInvestigator.h>

#include <QtCore/QBuffer>
#include <QtCore/QCoreApplication>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QEventLoop>
#include <QtNetwork/QNetworkAccessManager>
#include <QtNetwork/QNetworkReply>
#include <QtNetwork/QNetworkRequest>
#include <QtCore/QProcess>
#include <QtCore/QStringList>
#include <QtCore/QTextStream>
#include <QtCore/QtDebug>

#define VI_API_SUFFIX "/interface/API.php"
#define VI_SSH_HASH QString("c4:fc:76:94:f8:7a:dc:8a:64:f2:71:bc:46:20:85:99")
#define BUFFER_SIZE 65536

// From https://github.com/lindenb/cclindenb/blob/master/src/core/lindenb/io/tarball.cpp
struct PosixTarHeader
{
    char name[100];
    char mode[8];
    char uid[8];
    char gid[8];
    char size[12];
    char mtime[12];
    char checksum[8];
    char typeflag[1];
    char linkname[100];
    char magic[6];
    char version[2];
    char uname[32];
    char gname[32];
    char devmajor[8];
    char devminor[8];
    char prefix[155];
    char pad[12];
};

// Based on code from https://github.com/lindenb/cclindenb/blob/master/src/core/lindenb/io/tarball.cpp
uint headerChecksum(void* header)
{
  unsigned int sum = 0;
  char *p = (char *) header;
  char *q = p + sizeof(PosixTarHeader);

  while (p < static_cast<PosixTarHeader*>(header)->checksum)
  {
    sum += *p++ & 0xff;
  }

  for (int i = 0; i < 8; ++i)
  {
    sum += ' ';
    ++p;
  }

  while (p < q)
  {
    sum += *p++ & 0xff;
  }

  return sum;
}

using namespace std;

namespace OpenMS
{
  PeakInvestigator::PeakInvestigator(QObject* parent) :
    QObject(parent),
    DefaultParamHandler("PeakInvestigator"),
    ProgressLogger()
  {
    // set default parameter values
    defaults_.setValue("server", "secure.veritomyx.com", "Server address for PeakInvestigator (without https://)");
    defaults_.setValue("username", "USERNAME", "Username for account registered with Veritomyx");
    defaults_.setValue("password", "PASSWORD", "Password for account registered with Veritomyx");
    defaults_.setValue("account", "0", "Account number");

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();

    // set application-wide settings (e.g. for saving jobs)
    QCoreApplication::setOrganizationName("Veritomyx");
    QCoreApplication::setOrganizationDomain("veritomyx.com");
    QCoreApplication::setApplicationName("PeakInvestigator");

  }

  PeakInvestigator::~PeakInvestigator()
  {
  }

  void PeakInvestigator::run()
  {

    switch(mode_)
    {

    case SUBMIT:
      if (!initializeJob_())
      {
        break;
      }

      bundleScans_();
      if(uploadBundle_() && submitJob_())
      {
        experiment_.setMetaValue("veritomyx:server", server_);
        experiment_.setMetaValue("veritomyx:job", job_);
        experiment_.setMetaValue("veritomyx:sftp_username", sftp_username_);
        experiment_.setMetaValue("veritomyx:sftp_password", sftp_password_);
        QString filename;
        if (out_filename_ == String::EMPTY)
        {
          filename = in_filename_.toQString().section(".", 0, -2) + "." + job_ + ".mzML"; // insert job number into filename
        }
        else
        {
          filename = out_filename_.toQString();
        }
        file_.store(filename, experiment_);
      }
      break;

    case CHECK:
      checkJob_();
      break;

    case FETCH:
      server_ = experiment_.getMetaValue("veritomyx:server");
      job_ = experiment_.getMetaValue("veritomyx:job").toQString();
      sftp_username_ = experiment_.getMetaValue("veritomyx:sftp_username").toQString();
      sftp_password_ = experiment_.getMetaValue("veritomyx:sftp_password").toQString();

      if(!checkJob_()) // Seems we need to check STATUS before file is moved to SFTP drop after completion
      {
        break;
      }

      if(!downloadBundle_())
      {
        break;
      }

      if (extractScans_() == 0)
      {
        LOG_INFO << "No results were found for " << job_.toAscii().data() << endl;
      }

      // remove SFTP username/password from file
      experiment_.removeMetaValue("veritomyx:sftp_username");
      experiment_.removeMetaValue("veritomyx:sftp_password");
      QString filename;
      if (out_filename_ == String::EMPTY)
      {
        filename = in_filename_.toQString().section('.', 0, -3) + ".peaks.mzML";
        file_.store(job_ + ".peaks.mzML", experiment_);
      }
      else
      {
        filename = out_filename_.toQString();
      }
      file_.store(filename, experiment_);

      removeJob_();
      break;

    } //end switch

    shutdown();

  }

  bool PeakInvestigator::loadFromInputFilename(String input_filename)
  {
    in_filename_ = input_filename;
    file_.load(input_filename, experiment_);
    if (experiment_.empty())
    {
      LOG_ERROR << "The given file appears to not contain any m/z-intensity data points.";
      return false;
    }

    //check for peak type (profile data required)
    if (PeakTypeEstimator().estimateType(experiment_[0].begin(), experiment_[0].end()) == SpectrumSettings::PEAKS)
    {
      LOG_ERROR << "OpenMS peak type estimation indicates that this is not profile data!";
      return false;
    }

    return true;
  }

  QBuffer* PeakInvestigator::saveDataToBuffer_(MSSpectrum<Peak1D> &spectrum)
  {
    QBuffer* buffer = new QBuffer;
    buffer->open(QIODevice::WriteOnly | QIODevice::Text);

    QTextStream out(buffer);

    for (Size i = 0; i < spectrum.size(); i++)
    {
      out << QString("%1").arg(spectrum[i].getMZ(), 0, 'f', 6) << '\t';
      out << QString("%1").arg(spectrum[i].getIntensity(), 0, 'f', 6) << endl;
    }

    buffer->close();

    return buffer;
  }

  void PeakInvestigator::loadDataFromBuffer_(QBuffer *buffer, MSSpectrum<Peak1D> &peaklist)
  {
    buffer->open(QIODevice::ReadOnly);

    double mz, counts;
    QString line;
    bool ok1, ok2;

    QTextStream in(buffer);
    while(!in.atEnd())
    {
      line = in.readLine();

      // skip comments
      if (line[0] == '#')
      {
        continue;
      }

      mz = line.section('\t', 0, 0).toDouble(&ok1);
      counts = line.section('\t', 1, 1).toDouble(&ok2);

      if (!ok1 || !ok2)
      {
        LOG_WARN << "Problem decoding line to buffer: '" << line.toAscii().data() << "'.\n";
      }

      Peak1D peak;
      peak.setMZ(mz);
      peak.setIntensity(counts);
      peaklist.push_back(peak);
    }

    peaklist.updateRanges();
    buffer->close();

  }

#if defined(__APPLE__) || defined(linux)
  ssh_session PeakInvestigator::establishSSHSession_(QString hostname, QString username)
  {
    int retval;
    ssh_session session = ssh_new();

    if (session == NULL)
    {
      LOG_FATAL_ERROR << "Unable to create SSH session.\n";
    }
    else
    {
      ssh_options_set(session, SSH_OPTIONS_HOST, hostname.toAscii().data());
      ssh_options_set(session, SSH_OPTIONS_USER, username.toAscii().data());

      retval = ssh_connect(session);
      if (retval != SSH_OK)
      {
        LOG_FATAL_ERROR << "Unabled to connect to SSH session: " << ssh_get_error(session) << endl;
        ssh_free(session);
        session = NULL;
      }
    }

    return session;
  }

  bool PeakInvestigator::confirmSSHServerIdentity_(ssh_session session, QString exptected_hash)
  {
    ssh_key key;
    int retval;
    bool shouldProceed = false;

    retval = ssh_get_publickey(session, &key);

    if (retval != SSH_OK)
    {
      LOG_FATAL_ERROR << "Unable to get server public key: " << ssh_get_error(session) << endl;
      ssh_disconnect(session);
      ssh_free(session);
      return false;
    }

    unsigned char* key_hash;
    size_t key_len;
    int state;

    retval = ssh_get_publickey_hash(key, SSH_PUBLICKEY_HASH_MD5, &key_hash, &key_len);
    state = ssh_is_server_known(session);

    switch (state)
    {
    case SSH_SERVER_KNOWN_OK:
      shouldProceed = true;
      break;

    case SSH_SERVER_KNOWN_CHANGED:
      LOG_WARN << "\n***********************************************************************\n";
      LOG_WARN << "The identity of " << server_ << " has changed, which is potentially \n";
      LOG_WARN << "a security problem. Please contact Veritomyx support.";
      LOG_WARN << "\n***********************************************************************" << endl;
      break;

    case SSH_SERVER_FOUND_OTHER:
      LOG_WARN << "\n***********************************************************************\n";
      LOG_WARN << "The host key for " << server_ << " was not found, but another";
      LOG_WARN << "type of key exits. This is potentially a security problem.\n";
      LOG_WARN << "Please contact Veritomyx support.";
      LOG_WARN << "\n***********************************************************************" << endl;
      break;

    case SSH_SERVER_FILE_NOT_FOUND:
    case SSH_SERVER_NOT_KNOWN:
      char* found_hash = ssh_get_hexa(key_hash, key_len);
      LOG_WARN << "\n**************************************************************************\n";
      if (exptected_hash == found_hash)
      {

        LOG_WARN << "Adding the following key to the known hosts file for " << server_ << ":\n";
        LOG_WARN << exptected_hash.toAscii().data() << endl;
        ssh_write_knownhost(session);
        shouldProceed = true;

      }
      else
      {
        LOG_WARN << "The host key for " << server_ << " does not match expected.\n";
        LOG_WARN << "Expected: " << exptected_hash.toAscii().data() << ".\n";
        LOG_WARN << "Found: " << found_hash << ".\n";
        LOG_WARN << "\n";
        LOG_WARN << "Do you wish to proceed anyways (yes/no)?" << endl;

        String answer;
        cin >> answer;

        if (answer == "yes")
        {
          ssh_write_knownhost(session);
          shouldProceed = true;
        }

        else
        {
          LOG_WARN << "Exiting." << endl;
        }
      }

      LOG_WARN << "\n**************************************************************************" << endl;
      ssh_string_free_char(found_hash);
      break;
    }

    ssh_clean_pubkey_hash(&key_hash);
    return shouldProceed;
  }

  bool PeakInvestigator::authenticateUser_(ssh_session session, QString password)
  {
    int retval = ssh_userauth_password(session, NULL, password.toAscii().data());
    if (retval != SSH_AUTH_SUCCESS)
    {
      ssh_disconnect(session);
      ssh_free(session);
      qDebug() << "Problem authenticating user.\n";
      return false;
    }

    return true;
  }

  sftp_session PeakInvestigator::establishSFTPSession_(ssh_session session)
  {
    int retval;

    sftp_session sftp = sftp_new(session);
    if (sftp == NULL)
    {
      ssh_free(session);
      qDebug() << "Unable to create SFTP session.\n";
    }
    else
    {
      retval = sftp_init(sftp);
      if (retval != SSH_OK)
      {
        qDebug() << "Unable to inialize SFTP session: " << sftp_get_error(sftp) << endl;
        sftp_free(sftp);
        ssh_disconnect(session);
        ssh_free(session);
        sftp = NULL;
      }
    }

    return sftp;
  }

  bool PeakInvestigator::uploadFile_(sftp_session sftp, QString localFileName, QString remoteFileName)
  {
    int access = O_WRONLY | O_CREAT | O_TRUNC;
    sftp_file file;
    uint length, nwritten;
    char buffer[BUFFER_SIZE];

    file = sftp_open(sftp, remoteFileName.toAscii().data(),
                     access, S_IRWXU | S_IRWXG | S_IRWXO);
    if (file == NULL)
    {
      LOG_FATAL_ERROR << "Can't open remote file for writing: " << sftp_get_error(sftp) << endl;
      return false;
    }

    QFileInfo info(localFileName);
    uint size = info.size();
    uint progress = 0;

    startProgress(0, size, "Uploading scans to " + server_);
    FILE* fileObject = fopen(localFileName.toAscii().data(), "r+");

    do
    {
      length = fread(buffer, sizeof(char), BUFFER_SIZE, fileObject);

      if (length <= 0)
      {
        break;
      }

      nwritten = sftp_write(file, buffer, length);
      if (nwritten != length)
      {
        LOG_ERROR << "Can't write data to file: " << sftp_get_error(sftp) << endl;
        break;
      }

      progress += nwritten;
      setProgress(progress);

    } while(true);

    endProgress();

    fclose(fileObject);

    int retval = sftp_close(file);
    if (retval != SSH_OK)
    {
      LOG_FATAL_ERROR << "Can't close the written file: " << sftp_get_error(sftp) << endl;
      return false;
    }

    return true;
  }

  bool PeakInvestigator::downloadFile_(sftp_session sftp, QString localFileName, QString remoteFileName)
  {
    int access = O_RDONLY;
    sftp_file file;
    int length, nwritten;
    char buffer[BUFFER_SIZE];


    file = sftp_open(sftp, remoteFileName.toAscii().data(), access, 0);
    if (file == NULL)
    {
      LOG_FATAL_ERROR << "Can't open remote file for reading: " << sftp_get_error(sftp) << endl;
      return false;
    }

    FILE* fileObject = fopen(localFileName.toAscii().data(), "w");

    do
    {
      length = sftp_read(file, buffer, sizeof(buffer));
      //length = stream.readRawData(buffer, BUFFER_SIZE);
      //length = fread(buffer, sizeof(char), BUFFER_SIZE, fileObject);

      if (length == 0)
      {
        break; //EOF
      }
      else if (length < 0)
      {
        LOG_FATAL_ERROR << "Error while reading remote file: " << sftp_get_error(sftp);
        sftp_close(file);
        fclose(fileObject);
        return false;
      }

      //nwritten = sftp_write(file, buffer, length);
      nwritten = fwrite(buffer, sizeof(char), length, fileObject);
      if (nwritten != length)
      {
        qDebug() << "Can't write data to file: " << sftp_get_error(sftp) << endl;
        break;
      }

    } while(true);

    fclose(fileObject);

    int retval = sftp_close(file);
    if (retval != SSH_OK)
    {
      LOG_FATAL_ERROR << "Can't close the remote file: " << sftp_get_error(sftp) << endl;
      return false;
    }

    return true;
  }

#elif _WIN32
  void PeakInvestigator::displayPSCPError_(int error)
  {
    LOG_ERROR << "\n***********************************************************************\n";
    LOG_ERROR << "Problem with PSCP process.\n\n";
    switch(error)
    {
    case QProcess::FailedToStart:
      LOG_ERROR << "The PSCP executable is missing from your path. Please download it from\n";
      LOG_ERROR << "http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html, and either\n";
      LOG_ERROR << "add its directory to your PATH or copy it to the same directory from which\n";
      LOG_ERROR << "you call the PeakInvestigator tool.";
      break;
    case QProcess::Crashed:
      LOG_ERROR << "The PSCP program crashed. Please consult Veritomyx for support.";
      break;
    case QProcess::Timedout:
    case QProcess::UnknownError:
      LOG_ERROR << "There is an unknown or timeout error for starting the PSCP program.\n";
      LOG_ERROR << "Please contact Veritomyx for support.";
      break;
    case QProcess::ReadError:
    case QProcess::WriteError:
      LOG_ERROR << "There was a problem reading/writing to the PSCP program. Please contact\n";
      LOG_ERROR << "Veritomyx for support.";
      break;
    }
    LOG_ERROR << "\n***********************************************************************\n" << endl;
  }

#endif //APPLE or Linux

  void PeakInvestigator::bundleScans_(QString zipfilename)
  {
    int numBytes;
    char temp[8192];

    if(zipfilename.isEmpty())
    {
      zipfilename = QDir::tempPath() + "/" + job_ + ".scans.tar";
    }

    gzFile file = gzopen(zipfilename.toAscii().data(), "wb");
    PosixTarHeader header;

    startProgress(0, experiment_.size(), "Bundling scans for upload");
    for (Size i = 0; i < experiment_.size(); ++i)
    {
      QBuffer* buffer = saveDataToBuffer_(experiment_[i]);
      experiment_[i].clear(false);  // remove data, but keep metainfo for saving

      // initialize archive header
      std::memset(&header, 0, sizeof(PosixTarHeader)); // fill with NULL chars

      // set archive filename
      std::sprintf(header.name, "scan_%06i.txt", (int) i);

      // set tar format
      std::sprintf(header.magic, "ustar");
      std::sprintf(header.version, "  ");

      // set modification time, mode, and filetype
      std::sprintf(header.mtime, "%011lo", time(NULL));
      std::sprintf(header.mode, "%07o", 0644);
      header.typeflag[0] = 0;

      // set size of file
      qint64 size = buffer->size();
      std::sprintf(header.size, "%011llo", (long long unsigned int) size);

      // set header checksum
      std::sprintf(header.checksum, "%06o", headerChecksum(&header));

      // now write archive file header (required before writing data)
      numBytes = gzwrite(file, &header, sizeof(PosixTarHeader));
      if(numBytes != sizeof(PosixTarHeader))
      {
        LOG_ERROR << "Not all of the header was written for scan " << i << "!\n";
        delete buffer;
        continue;
      }

      // copy buffer to file in archive
      buffer->open(QBuffer::ReadOnly);
      numBytes = buffer->read(temp, sizeof(temp));
      while (numBytes > 0)
      {
        gzwrite(file, temp, numBytes);
        numBytes = buffer->read(temp, sizeof(temp));
      }
      buffer->close();

      // fill remaining 512-byte block with NULL
      while (size % sizeof(PosixTarHeader) != 0)
      {
        gzputc(file, 0);
        ++size;
      }

      delete buffer;

      // occasionally update progress
      if (i % 10 == 0)
      {
        setProgress(i);
      }
    }

    // now close out tar format by writing with two NULL header entries
    std::memset(&header, 0, sizeof(PosixTarHeader)); // fill with NULL chars
    gzwrite(file, &header, sizeof(PosixTarHeader));
    gzwrite(file, &header, sizeof(PosixTarHeader));

    gzclose(file);

    endProgress();

  }

  int PeakInvestigator::extractScans_(QString zipfilename)
  {
    char temp[8192];
    int numBytes, retval;

    if (zipfilename.isEmpty())
    {
      zipfilename = QDir::tempPath() + "/" + job_ + ".vcent.tar";
    }

    gzFile file = gzopen(zipfilename.toAscii().data(), "rb");

    if (file == NULL)
    {
      LOG_FATAL_ERROR << "Problem opening " << zipfilename.toAscii().data() << endl;
      return 0;
    }

    // set-up DataProcessing metadata to add to each spectrum
    DataProcessing dp;
    std::set<DataProcessing::ProcessingAction> actions;
    actions.insert(DataProcessing::PEAK_PICKING);
    dp.setProcessingActions(actions);
    dp.getSoftware().setName("PeakInvestigator");
    dp.setCompletionTime(DateTime::now());
    dp.setMetaValue("parameter: veritomyx:server", server_);
    dp.setMetaValue("parameter: veritomyx:username", username_);
    dp.setMetaValue("parameter: veritomyx:account", account_number_);
    dp.setMetaValue("veritomyx:job", job_);

    int i = 0;
    do
    {
      // read header of entry i
      PosixTarHeader header;
      numBytes = gzread(file, &header, sizeof(PosixTarHeader));

      if (numBytes == 0) // end of file
      {
        break;
      }

      else if (numBytes != sizeof(PosixTarHeader))
      {
        LOG_ERROR << "Problem parsing header for entry " << i << endl;
        return i;
      }

      bool ok;
      qlonglong fileSize = QString(header.size).toLongLong(&ok, 8);
      if (!ok) // ugly hack to handle end of tarfile, as well as bad-formed header
      {
        continue;
      }
      int remainder = sizeof(PosixTarHeader) - (fileSize % sizeof(PosixTarHeader)); // entries are NULL terminated

      QBuffer* buffer = new QBuffer;
      buffer->open(QIODevice::WriteOnly);

      while (fileSize >= sizeof(temp))
      {
        numBytes = gzread(file, temp, sizeof(temp));
        buffer->write(temp, numBytes);
        fileSize -= numBytes;
      }
      // get leftover
      gzread(file, temp, fileSize);
      buffer->write(temp, fileSize);

      // read rest of entry (NULLs)
      gzread(file, temp, remainder);

      buffer->close();

      // now load the buffer into a Spectrum, and place it in the MSExperiment
      loadDataFromBuffer_(buffer, experiment_[i]);
      experiment_[i].getDataProcessing().push_back(dp);
      experiment_[i].setType(SpectrumSettings::PEAKS);

      delete buffer;
      i++;


    } while (true);

    qDebug() << "Processed" << i << "files";

    gzclose(file);

    return i;
  }

  bool PeakInvestigator::uploadBundle_()
  {

#if defined(__APPLE__) || defined(linux)
    ssh_session session = establishSSHSession_(server_.toQString(), sftp_username_);
    if (session == NULL)
    {
      return false;
    }

    if (! confirmSSHServerIdentity_(session, VI_SSH_HASH))
    {
      return false;
    }

    if (! authenticateUser_(session, sftp_password_))
    {
      return false;
    }

    sftp_session sftp = establishSFTPSession_(session);
    if (sftp == NULL)
    {
      return false;
    }

    QString tarFile = job_ + ".scans.tar";
    QString localFile = QDir::tempPath() + "/" + tarFile;
    QString remoteFile = "accounts/" + account_number_.toQString() + "/batches/" + tarFile;

    qDebug() << "Uploading" << localFile << "to" << remoteFile;
    uploadFile_(sftp, localFile, remoteFile);

    ssh_disconnect(session);
    ssh_free(session);
    return true;

#elif _WIN32
    QString tarFile = job_ + ".scans.tar";
    QString localFile = QDir::tempPath() + "/" + tarFile;
    QString remoteFile = "accounts/" + account_number_.toQString() + "/batches/" + tarFile;
    QString remoteURL = server_.toQString() + ":" + remoteFile;

    QProcess process;
    process.setProcessChannelMode(QProcess::ForwardedChannels);

    QStringList arguments;
    arguments << "-l" << sftp_username_ << "-pw" << sftp_password_;
    arguments << localFile << remoteURL;

    process.start("pscp", arguments);

    if (!process.waitForStarted())
    {
      displayPSCPError_((int) process.error());
      return false;
    }

    if (!process.waitForFinished(-1))
    {
      displayPSCPError_((int) process.error());
      return false;
    }

    return true;

#else
    LOG_ERROR << "Platform not correctly identified or supported. Exiting.";
    return false;

#endif // APPLE or Linux

  }

  bool PeakInvestigator::downloadBundle_()
  {

#if defined(__APPLE__) || defined(linux)
    ssh_session session = establishSSHSession_(server_.toQString(), sftp_username_);
    if (session == NULL)
    {
      return false;
    }

    if (! confirmSSHServerIdentity_(session, VI_SSH_HASH))
    {
      return false;
    }

    if (! authenticateUser_(session, sftp_password_))
    {
      return false;
    }

    sftp_session sftp = establishSFTPSession_(session);
    if (sftp == NULL)
    {
      return false;
    }

    QString tarFile = job_ + ".vcent.tar";
    QString localFile = QDir::tempPath() + "/" + tarFile;
    QString remoteFile = "accounts/" + account_number_.toQString() + "/results/" + tarFile;

    qDebug() << "Downloading" << remoteFile << "to" << localFile;
    downloadFile_(sftp, localFile, remoteFile);

    ssh_disconnect(session);
    ssh_free(session);
    return true;

#elif _WIN32
    QString tarFile = job_ + ".vcent.tar";
    QString localFile = QDir::tempPath() + "/" + tarFile;
    QString remoteFile = "accounts/" + account_number_.toQString() + "/results/" + tarFile;
    QString remoteURL = server_.toQString() + ":" + remoteFile;

    QProcess process;
    process.setProcessChannelMode(QProcess::ForwardedChannels);

    QStringList arguments;
    arguments << "-l" << sftp_username_ << "-pw" << sftp_password_;
    arguments << remoteURL << localFile;

    process.start("pscp", arguments);

    if (!process.waitForStarted())
    {
      displayPSCPError_((int) process.error());
      return false;
    }

    if (!process.waitForFinished(-1))
    {
      displayPSCPError_((int) process.error());
      return false;
    }

    return true;


#else
    LOG_ERROR << "Platform not correctly identified or supported. Exiting.";
    return false;

#endif // Apple or Linux

  }

  bool PeakInvestigator::initializeJob_()
  {
    LOG_DEBUG << "Requsting credentials for " + username_ + "..." << endl;

    url_.setUrl("https://" + server_.toQString() + VI_API_SUFFIX);
    url_.addQueryItem("Version", "1.25");
    url_.addQueryItem("User", username_.toQString());
    url_.addQueryItem("Code", password_.toQString());
    url_.addQueryItem("Action", "INIT");
    url_.addQueryItem("Account", account_number_.toQString());
    url_.addQueryItem("Command", "ckm");
    url_.addQueryItem("Count", QString::number(experiment_.size()));

    QNetworkRequest request(url_);
    reply_ = manager_.get(request);

    QEventLoop loop;
    QObject::connect(reply_, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();

    if (reply_->error() != QNetworkReply::NoError)
    {
      LOG_ERROR << "There was an error making a network request:\n";
      LOG_ERROR << reply_->errorString().toAscii().data() << endl;
      reply_->deleteLater();
      return false;
    }

    QString contents(reply_->readAll());
    reply_->deleteLater();

    if (contents.startsWith("Error"))
    {
      QStringList list = contents.split(":");
      LOG_ERROR << "Error occurred:" << list[1].toAscii().data() << endl;
      return false;
    }
    else if (contents.startsWith("<!DOCTYPE HTML"))
    {
      LOG_ERROR << "There is a problem with the specified server address." << endl;
      return false;
    }

    QStringList list = contents.split(" ");
    job_ = list[2];
    sftp_username_ = list[3];
    sftp_password_ = list[4];

    return true;
  }

  bool PeakInvestigator::submitJob_()
  {
    url_.setUrl("https://" + server_.toQString() + VI_API_SUFFIX);
    url_.addQueryItem("Version", "1.25");
    url_.addQueryItem("User", username_.toQString());
    url_.addQueryItem("Code", password_.toQString());
    url_.addQueryItem("Action", "RUN");
    url_.addQueryItem("Job", job_);

    QNetworkRequest request(url_);
    reply_ = manager_.get(request);

    QEventLoop loop;
    QObject::connect(reply_, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();

    if (reply_->error() != QNetworkReply::NoError)
    {
      LOG_ERROR << "There was an error making a network request:\n";
      LOG_ERROR << reply_->errorString().toAscii().data() << endl;
      reply_->deleteLater();
      return false;
    }

    QString contents(reply_->readAll());
    reply_->deleteLater();

    if (contents.startsWith("Error")) {
      QStringList list = contents.split(":");
      cout << "Error occurred:" << list[1].toAscii().data() << endl;
      return false;
    }

    cout << contents.toAscii().data() << endl;
    return true;

  }

  bool PeakInvestigator::checkJob_()
  {
    bool retval = false;

    server_ = experiment_.getMetaValue("veritomyx:server");
    job_ = experiment_.getMetaValue("veritomyx:job").toQString();

    if (job_.isEmpty())
    {
      LOG_WARN << "Problem getting job ID from meta data.\n";
      return retval;
    }

    url_.setUrl("https://" + server_.toQString() + VI_API_SUFFIX);
    url_.addQueryItem("Version", "1.25");
    url_.addQueryItem("User", username_.toQString());
    url_.addQueryItem("Code", password_.toQString());
    url_.addQueryItem("Action", "STATUS");
    url_.addQueryItem("Job", job_);

    QNetworkRequest request(url_);
    reply_ = manager_.get(request);

    QEventLoop loop;
    QObject::connect(reply_, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();

    if (reply_->error() != QNetworkReply::NoError)
    {
      LOG_ERROR << "There was an error making a network request:\n";
      LOG_ERROR << reply_->errorString().toAscii().data() << endl;
      reply_->deleteLater();
      return false;
    }

    QString contents(reply_->readAll());
    reply_->deleteLater();

    if (contents.startsWith("Error"))
    {
      QStringList list = contents.split(":");
      cout << "Error occurred:" << list[1].toAscii().data() << endl;
      retval = false;
    }
    else if (contents.startsWith("Running"))
    {
      LOG_INFO << job_.toAscii().data() << " is still running.\n";
      retval = false;
    }
    else if (contents.startsWith("Done"))
    {
      LOG_INFO << job_.toAscii().data() << " has finished.\n";
      retval = true;
    }

    return retval;
  }

  bool PeakInvestigator::removeJob_()
  {
    url_.setUrl("https://" + server_.toQString() + VI_API_SUFFIX);
    url_.addQueryItem("Version", "1.25");
    url_.addQueryItem("User", username_.toQString());
    url_.addQueryItem("Code", password_.toQString());
    url_.addQueryItem("Action", "DONE");
    url_.addQueryItem("Job", job_);

    QNetworkRequest request(url_);
    reply_ = manager_.get(request);

    QEventLoop loop;
    QObject::connect(reply_, SIGNAL(finished()), &loop, SLOT(quit()));
    loop.exec();

    if (reply_->error() != QNetworkReply::NoError)
    {
      LOG_ERROR << "There was an error making a network request:\n";
      LOG_ERROR << reply_->errorString().toAscii().data() << endl;
      reply_->deleteLater();
      return false;
    }

    QString contents(reply_->readAll());
    reply_->deleteLater();

    if (contents.startsWith("Error")) {
      QStringList list = contents.split(":");
      cout << "Error occurred:" << list[1].toAscii().data() << endl;
      return false;
    }

    cout << contents.toAscii().data() << endl;
    return true;

  }

  void PeakInvestigator::updateMembers_()
  {
    server_ = param_.getValue("server");
    username_ = param_.getValue("username");
    password_ = param_.getValue("password");
    account_number_ = param_.getValue("account");
  }

}

