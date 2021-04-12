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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtCore/QVector>

namespace OpenMS
{
  class TOPPASScene;

  /**
      @brief A vertex representing a TOPP tool

      Besides TOPPASScene, this class contains most of the remaining functionality of
      TOPPAS regarding the execution of pipelines. Once a pipeline run is started
      from TOPPASScene, the execution is propagated from tool to tool and the
      TOPP tools are actually called from here.

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASToolVertex :
    public TOPPASVertex
  {
    Q_OBJECT

public:
    /// current status of the vertex
    enum TOOLSTATUS {TOOL_READY, TOOL_SCHEDULED, TOOL_RUNNING, TOOL_SUCCESS, TOOL_CRASH, TOOLSTATUS_SIZE};

    /// Stores the information for input/output files/lists
    struct IOInfo
    {
      /// Standard constructor
      IOInfo() :
        type(IOT_FILE),
        param_name(),
        valid_types()
      {
      }

      /// Copy constructor
      IOInfo(const IOInfo& rhs) :
        type(rhs.type),
        param_name(rhs.param_name),
        valid_types(rhs.valid_types)
      {
      }

      /// The type
      enum IOType
      {
        IOT_FILE,
        IOT_LIST
      };

      /// Comparison operator
      bool operator<(const IOInfo& rhs) const
      {
        if (type != rhs.type)
        {
          return type == IOT_FILE;
        }
        else
        {
          return param_name.compare(rhs.param_name) < 0;
        }
      }

      /// Assignment operator
      IOInfo& operator=(const IOInfo& rhs)
      {
        type = rhs.type;
        param_name = rhs.param_name;
        valid_types = rhs.valid_types;

        return *this;
      }

      /// Is any of the input/output parameters a list?
      static bool isAnyList(const QVector<IOInfo>& params)
      {
        for (QVector<IOInfo>::const_iterator it = params.begin();
             it != params.end(); ++it)
        {
          if (it->type == IOT_LIST) return true;
        }
        return false;
      }

      ///The type of the parameter
      IOType type;
      ///The name of the parameter
      String param_name;
      ///The valid file types for this parameter
      StringList valid_types;
    };

    /// Default constructor
    TOPPASToolVertex();
    /// Constructor
    TOPPASToolVertex(const String& name, const String& type = "");
    /// Copy constructor
    TOPPASToolVertex(const TOPPASToolVertex& rhs);
    /// Destructor
    ~TOPPASToolVertex() override = default;
    /// Assignment operator
    TOPPASToolVertex& operator=(const TOPPASToolVertex& rhs);

    /// returns the name of the TOPP tool
    String getName() const override;
    /// Returns the type of the tool
    const String& getType() const;
    /// Fills @p input_infos with the required input file/list parameters together with their valid types.
    void getInputParameters(QVector<IOInfo>& input_infos) const;
    /// Fills @p output_infos with the required output file/list parameters together with their valid types.
    void getOutputParameters(QVector<IOInfo>& output_infos) const;
    // documented in base class
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) override;
    // documented in base class
    QRectF boundingRect() const override;
    // documented in base class
    void setTopoNr(UInt nr) override;
    // documented in base class
    void reset(bool reset_all_files = false) override;
    /// Sets the Param object of this tool
    void setParam(const Param& param);
    /// Returns the Param object of this tool
    const Param& getParam();
    /// Checks if all parent nodes have finished the tool execution and, if so, runs the tool
    void run() override;
    /// Updates the vector containing the lists of current output files for all output parameters
    /// using the input files as guidance
    /// Returns true on success, on failure the error_message is filled
    bool updateCurrentOutputFileNames(const RoundPackages& pkg, String& error_message);
    /// return if tool failed or is ready etc.
    TOOLSTATUS getStatus() const;
    /// Lets the user edit the parameters of the tool
    void editParam();
    /// Returns the number of iterations this tool has to perform
    int numIterations();
    /// Returns the full directory (including preceding tmp path)
    String getFullOutputDirectory() const;
    /// Returns the directory where this tool stores its output files
    String getOutputDir() const;
    /// Creates all necessary directories
    void createDirs();
    /// Opens the folder where the file is contained
    void openContainingFolder();
    /// Opens the files in TOPPView
    void openInTOPPView();
    /// Refreshes the parameters of this tool, returns if their has been a change
    bool refreshParameters();
    /// underlying TOPP tool found and parameters fetched?! (done in C'Tor)
    bool isToolReady() const;
    /// Toggle breakpoint
    void toggleBreakpoint();
    /// Called when the QProcess in the queue is called: emits 'toolStarted()'
    virtual void emitToolStarted();
    /// invert status of recycling (overriding base class)
    bool invertRecylingMode() override;

public slots:

    /// Called when the execution of this tool has finished
    void executionFinished(int ec, QProcess::ExitStatus es);
    /// Called when the running TOPP tool produces output
    void forwardTOPPOutput();
    /// Called when the tool is started
    void toolStartedSlot();
    /// Called when the tool has finished
    void toolFinishedSlot();
    /// Called when the tool has crashed
    void toolCrashedSlot();
    /// Called when the tool has failed
    void toolFailedSlot();
    /// Called when the tool was scheduled for running
    virtual void toolScheduledSlot();
    /// Called by an incoming edge when it has changed
    void inEdgeHasChanged() override;
    /// Called by an outgoing edge when it has changed
    void outEdgeHasChanged() override;

signals:

    /// Emitted when the tool is started
    void toolStarted();
    /// Emitted when the tool is finished
    void toolFinished();
    /// Emitted when the tool crashes
    void toolCrashed();
    /// Emitted when the tool execution fails
    void toolFailed(const QString& message = "");
    /// Emitted from forwardTOPPOutput() to forward the signal outside
    void toppOutputReady(const QString& out);

protected:

    ///@name reimplemented Qt events
    //@{
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e) override;
    //@}
	

    /// get parent Scene
    TOPPASScene* getScene_() const;

    /// determines if according to current status_, a parameter change would invalidate the pipeline status (e.g., because this node was already processed)
    bool doesParamChangeInvalidate_();
    /// renames SUFFICES of the output files created by the TOPP tool by inspecting file content
    bool renameOutput_();
    /// Initializes the parameters with standard values (from -write_ini), uses the parameters from the old_ini_file if given, returns if parameters have changed (if old_ini_file was given)
    bool initParam_(const QString& old_ini_file = "");
    /// Fills @p io_infos with the required input/output file/list parameters. If @p input_params is true, input params are returned, otherwise output params.
    void getParameters_(QVector<IOInfo>& io_infos, bool input_params) const;
    /// Writes @p param to the @p ini_file
    void writeParam_(const Param& param, const QString& ini_file);
    /// Helper method for finding good boundaries for wrapping the tool name. Returns a string with whitespaces at the preferred boundaries.
    QString toolnameWithWhitespacesForFancyWordWrapping_(QPainter* painter, const QString& str);
    /// smart naming of round-based filenames
    /// when basename is not unique we take the preceding directory name
    void smartFileNames_(std::vector<QStringList>& filenames);

    /// The name of the tool
    String name_;
    /// The type of the tool, or "" if it does not have a type
    String type_;
    /// The temporary path
    String tmp_path_;
    /// The parameters of the tool
    Param param_;
    /// current status of the tool
    TOOLSTATUS status_{TOOL_READY};
    /// tool initialization status: if C'tor was successful in finding the TOPP tool, this is set to 'true'
    bool tool_ready_{true};
    /// Breakpoint set?
    bool breakpoint_set_{false};
  };
}

