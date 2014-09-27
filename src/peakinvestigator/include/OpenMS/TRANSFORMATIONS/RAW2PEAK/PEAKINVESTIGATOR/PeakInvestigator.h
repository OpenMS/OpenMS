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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_H

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/FORMAT/TarFile.h>

#include <QtNetwork/QNetworkAccessManager>
#include <QtCore/QUrl>

#if defined(Q_OS_MAC) || defined(Q_OS_LINUX)
  #include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/SYSTEM/LibSSH2SecureFileTransfer.h>
#elif defined(Q_OS_WIN32)
  #include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/SYSTEM/PSCPSecureFileTransfer.h>
#endif

class QBuffer;
class QNetworkReply;

namespace OpenMS
{
  /**
    @brief This class implements the <a href=https://secure.veritomyx.com/interface/API.php>
    PeakInvestigator public API</a> provided by <a href=http://www.veritomyx.com>Veritomyx</a>.

    This class has three modes of operation (submit, check, and fetch) which are
    specified by the PeakInvestigator::setMode() function. Since this class depends
    on a Qt event loop for some network operations, a QCoreApplication (or QApplication) needs to be
    instatiated and its exec() function called before PeakInvestigator::run() is executed. For example:

    @code
    QCoreApplication app(argc, argv2);

    PeakInvestigator pp(&app);
    pp.setLogType(log_type_);
    pp.setParameters(pepi_param); // set parameters containing username, password, and account number

    if(!pp.loadFromInputFilename(in)) // make sure the file has correct type/format
    {
        return TOPPBase::INCOMPATIBLE_INPUT_DATA;
    }

    pp.setOutputFilename(out); // set output name, which is used in PeakInvestigator::run()

    if(mode == "submit")
    {
        pp.setMode(PeakInvestigator::SUBMIT);
    }

    else if(mode == "check")
    {
        pp.setMode(PeakInvestigator::CHECK);
    }

    else if(mode == "fetch")
    {
        pp.setMode(PeakInvestigator::FETCH);
    }


    QObject::connect(&pp, SIGNAL(finishedRequest()), &app, SLOT(quit()));

    // delay calling PeakInvestigator::run() until the event loop is established
    QTimer::singleShot(100, &pp, SLOT(run()));

    app.exec();

    @endcode

    @ingroup PeakPicking
  */

  class OPENMS_DLLAPI PeakInvestigator :
    public QObject,
    public DefaultParamHandler,
    public ProgressLogger
  {
    Q_OBJECT

public:

    enum PIMode
    {
        SUBMIT,
        CHECK,
        FETCH
    };

    /// Constructor
    PeakInvestigator(QObject* parent = 0);

    /// Destructor
    virtual ~PeakInvestigator();
    
    /// Set the mode to one of the following: SUBMIT, CHECK, FETCH
    void setMode(PIMode mode) { mode_ = mode; }

    /** @brief Set the experiment for processing
     *
     * This function also makes sure the experiment being set has the correct data (i.e.
     * contains mass spectra and isn't already centroided).
     * @param experiment  Self-explanatory.
     * @return Bool indicating whether the experiment has the correct attributes.
     */
    bool setExperiment(MSExperiment<Peak1D>& experiment);

    /// Get the experiment after processing
    MSExperiment<Peak1D>& getExperiment() { return experiment_; }

    /// Get the jobID
    QString getJobID() { return job_; }

    /// Function that should be called to exit the Qt event loop.
    void shutdown() { emit finishedRequest(); }


public slots:
    /// Main function that should be called once a Qt event loop has been created and is executing.
    /// @todo When changing to GAMMA, no longer need to store SFTP credentials to file.
    void run();

signals:
    /// Signal that should be connected to the quit() slot of a QApplication or QCoreApplication.
    void finishedRequest();

protected:

//--------------------------------------------------------------------------------------------------------
// Bundling/Extracting and SFTP Upload/Download functions
//--------------------------------------------------------------------------------------------------------
    /** @name Packaging functions
     * Used for bundling/unbundling files, and trasmission to/from SFTP server. Makes calls to the
     * SFTP-related and Spectrum <--> QBuffer functions.
     */
///@{

    /// Upload the bundle containing scans to the Veritomyx SFTP directory.
    bool uploadBundle_();

    /// Download the bundle containing results from the Veritomyx SFTP directory.
    bool downloadBundle_();

///@}
//--------------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------
// PeakInvestigator public API functions
//--------------------------------------------------------------------------------------------------------
    /** @name PeakInvestigator public API functions
     * Used for calling the <a href="https://secure.veritomyx.com/interface/API.php">PeakInvestigator
     * public API</a> at various stages of SUBMIT, CHECK, or FETCH modes.
     */
///@{

    /** @brief Initialize a job using the <a href="https://secure.veritomyx.com/interface/API.php">
     * PeakInvestigator API</a>.
     *
     * This requires the Veritomyx username, password, and account parameters to be
     * correctly specified with setParameters(). It sets the sftp_username_, sftp_password_, and
     * job_ variables.
     */
    bool initializeJob_();

    /** @brief Submit the job using the <a href="https://secure.veritomyx.com/interface/API.php">
     * PeakInvestigator API</a>.
     *
     * This assumes that a job ID has been assigned and that SFTP credentials have been obtained
     * (i.e. initializeJob() should have been called), and that the scans have been bundled and
     * uploaded to the Veritomyx SFTP drop.
     */
    bool submitJob_();

    /** @brief Check the job status using the <a href="https://secure.veritomyx.com/interface/API.php">
     * PeakInvestigator API</a>.
     *
     * This assumes that the Veritomyx username, password, and account variables have been correctly
     * specified with setParameters(). It obtains the job meta data value from the inputfile.
     */
    bool checkJob_();

    /** @brief Remove the job using the <a href="https://secure.veritomyx.com/interface/API.php">
     * PeakInvestigator API</a>.
     *
     * This removes the job from Veritomyx's servers, so should only be called once a job has finished
     * and scans have been downloaded (unless the job *should* be deleted without obtaining results).
     * It assumes that the Veritomyx username, password, and account variables have been correctly
     * specified with setParameters().
     */
    bool removeJob_();

///@}
//--------------------------------------------------------------------------------------------------------


    // Veritomyx account info
    String server_; ///< @brief Veritomyx server address. Should be provided using the TOPP interface.
    String username_; ///< @brief Veritomyx account username. Should be provided using the TOPP interface.
    String password_; ///< @brief Veritomyx account password. Should be provided using the TOPP interface.
    String account_number_; ///< @brief Veritomyx account number. Should be provided using the TOPP interface.
    QString job_; ///< @brief Job number obtained from public API during INIT request.
    QString sftp_username_; ///< @brief Username for Veritomyx SFTP server, obtained from public API.
    QString sftp_password_; ///< @brief Password for Veritomyx SFTP server, obtained from public API.

    // docu in base class
    void updateMembers_();

    // Network Stuff
    QNetworkAccessManager manager_; ///< @brief Class used for making request to public API.
    QNetworkReply* reply_; ///< @brief Pointer to response from public API.
    QUrl url_; ///< @brief Url of the public API.

    // Utility classes
    TarFile tar;
#if defined(Q_OS_MAC) || defined(Q_OS_LINUX)
    LibSSH2SecureFileTransfer sftp;
#elif defined(Q_OS_WIN32)
    PSCPSecureFileTransfer sftp;
#endif

    // Misc variables
    MSExperiment<Peak1D> experiment_; ///< @brief Class used to hold spectra (raw or peak data) in memory.
    PIMode mode_;

  }; // end PeakInvestigator

} // namespace OpenMS

#endif
