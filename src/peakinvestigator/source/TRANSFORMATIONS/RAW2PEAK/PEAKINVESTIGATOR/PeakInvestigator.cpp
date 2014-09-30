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
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/PeakInvestigator.h>

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
#define VI_SSH_HASH String("7E:6D:03:89:68:38:0B:9F:C7:E5:13:26:56:46:08:FF")

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

  }

  PeakInvestigator::~PeakInvestigator()
  {
  }

  void PeakInvestigator::run()
  {

    // filenames for the tar'd scans/results
    QString zipfilename;
    QString localFilename;
    QString remoteFilename;

    switch(mode_)
    {

    case SUBMIT:
      if (!initializeJob_())
      {
        break;
      }

      // Generate local and remote filenames of tar'd scans
      zipfilename = job_ + ".scans.tar";
      localFilename = QDir::tempPath() + "/" + zipfilename;
      remoteFilename = "accounts/" + account_number_.toQString() + "/batches/" + zipfilename;
      tar.store(localFilename, experiment_);

      // Remove data values from scans in exp now that they have been bundled
      for (Size i = 0; i < experiment_.size(); i++)
      {
        experiment_[i].clear(false);
      }

      // Set SFTP host paramters and upload file
      sftp.setHostname(server_);
      sftp.setUsername(sftp_username_);
      sftp.setPassword(sftp_password_);
      sftp.setExpectedServerHash(VI_SSH_HASH);

      if(sftp.uploadFile(localFilename, remoteFilename) && submitJob_())
      {
        experiment_.setMetaValue("veritomyx:server", server_);
        experiment_.setMetaValue("veritomyx:job", job_);
        experiment_.setMetaValue("veritomyx:sftp_username", sftp_username_);
        experiment_.setMetaValue("veritomyx:sftp_password", sftp_password_);
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

      // Set SFTP host paramters and upload file
      sftp.setHostname(server_);
      sftp.setUsername(sftp_username_);
      sftp.setPassword(sftp_password_);
      sftp.setExpectedServerHash(VI_SSH_HASH);

      // Generate local and remote filenames of tar'd scans
      zipfilename = job_ + ".vcent.tar";
      localFilename = QDir::tempPath() + "/" + zipfilename;
      remoteFilename = "accounts/" + account_number_.toQString() + "/results/" + zipfilename;

      if (!sftp.downloadFile(remoteFilename, localFilename))
      {
        break;
      }

      tar.load(localFilename, experiment_);

      // Set-up data processing meta data to add to each scan
      DataProcessing dp;
      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(DataProcessing::PEAK_PICKING);
      dp.setProcessingActions(actions);
      dp.getSoftware().setName("PeakInvestigator");
      dp.setCompletionTime(DateTime::now());
      dp.setMetaValue("paramter: veritomyx:server", server_);
      dp.setMetaValue("paramter: veritomyx:username", username_);
      dp.setMetaValue("parameter: veritomyx:account", account_number_);
      dp.setMetaValue("veritomyx:job", job_);

      // Now add meta data to the scans
      for (Size i = 0; i < experiment_.size(); i++)
      {
        experiment_[i].getDataProcessing().push_back(dp);
        experiment_[i].setType(SpectrumSettings::PEAKS);
      }

      // remove SFTP username/password from file
      experiment_.removeMetaValue("veritomyx:sftp_username");
      experiment_.removeMetaValue("veritomyx:sftp_password");
      removeJob_();
      break;

    } //end switch

    shutdown();

  }

  bool PeakInvestigator::setExperiment(MSExperiment<Peak1D>& experiment)
  {
    if (experiment.empty())
    {
      LOG_ERROR << "The given file appears to not contain any m/z-intensity data points.";
      return false;
    }

    //check for peak type (profile data required)
    if (PeakTypeEstimator().estimateType(experiment[0].begin(), experiment[0].end()) == SpectrumSettings::PEAKS)
    {
      LOG_ERROR << "OpenMS peak type estimation indicates that this is not profile data!";
      return false;
    }

    experiment_ = experiment;
    return true;
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

