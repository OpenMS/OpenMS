// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASSCENE_H
#define OPENMS_VISUAL_TOPPASSCENE_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QGraphicsScene>
#include <QtCore/QProcess>

namespace OpenMS
{
	class TOPPASVertex;
	class TOPPASToolVertex;
  class TOPPASMergerVertex;
	class TOPPASOutputFileListVertex;
	class TOPPASEdge;
	class TOPPASResources;
	
  /**
    @brief A FakeProcess class.
  */
  class FakeProcess
    : public QProcess
  {
    Q_OBJECT

    public:
      virtual void 	start ( const QString & program, const QStringList & arguments, OpenMode mode = ReadWrite );
  };

	/**
		@brief A container for all visual items of a TOPPAS workflow
		
		TOPPASScene is a subclass of QGraphicsScene and acts as a container for all visual items
		(i.e. all vertices and edges). It is visualized by a TOPPASWidget (a subclass of QGraphicsView).
		This class also provides large parts of the functionality of TOPPAS, e.g., the methods for loading,
		saving, running, and aborting pipelines are located here.
		
		TOPPASScene can also be used without a visualizing TOPPASWidget (i.e., without a gui) which can
		be indicated via the constructor. In this case, the signals for log message output are connected
		to standard out. This is utilized for the ExecutePipeline tool.
	
    Temporary files of the pipeline are stored in the member tmp_path_. Update it when loading a pipeline which has
    tmp data from an old run. TOPPASToolVertex will ask its parent scene() whenever it wants to know the tmp directory.

		@ingroup TOPPAS_elements
	*/
	class OPENMS_GUI_DLLAPI TOPPASScene
		:	public QGraphicsScene
	{
		Q_OBJECT
		
		public:
			
			/// Stores the information for a TOPP process
			struct TOPPProcess
			{
				///Constructor
				TOPPProcess(QProcess* p, const QString& cmd, const QStringList& arg)
					:	proc(p),
						command(cmd),
						args(arg)
				{
				}
				
				///The process
				QProcess* proc;
				///The command
				QString command;
				///The arguments
				QStringList args;
			};
			
			/// The current action mode (creation of a new edge, or panning of the widget)
			enum ActionMode
      {
      	AM_NEW_EDGE,
      	AM_MOVE
      };
      
      /// The container for edges
      typedef QList<TOPPASEdge*> EdgeContainer;
      /// A mutable iterator for edges
			typedef EdgeContainer::iterator EdgeIterator;
			/// A const iterator for edges
			typedef EdgeContainer::const_iterator ConstEdgeIterator;
			/// The container for vertices
			typedef QList<TOPPASVertex*> VertexContainer;
			/// A mutable iterator for vertices
			typedef VertexContainer::iterator VertexIterator;
			/// A const iterator for vertices
			typedef VertexContainer::const_iterator ConstVertexIterator;
			
			/// Constructor
			TOPPASScene(QObject* parent, const QString& tmp_path, bool gui = true);
			
			/// Destructor
			virtual ~TOPPASScene();
			
			/// Adds a vertex
			void addVertex(TOPPASVertex* tv);
			/// Adds an edge
			void addEdge(TOPPASEdge* te);
			/// Sets the action mode
			void setActionMode(ActionMode mode);
			/// Returns the action mode
			ActionMode getActionMode();
			/// Returns begin() iterator of all vertices
			VertexIterator verticesBegin();
			/// Returns end() iterator of all vertices
			VertexIterator verticesEnd();
			/// Returns begin() iterator of all edges
			EdgeIterator edgesBegin();
			/// Returns end() iterator of all edges
			EdgeIterator edgesEnd();
			/// Copies all currently selected edges and vertices
			void copySelected();
			/// Pastes the copied items
			void paste(QPointF pos = QPointF());
			/// Removes all currently selected edges and vertices
			void removeSelected();
			/// Unselects all items
			void unselectAll();
			/// Updates all edge colors (color of green and yellow edges can change when edges are added/removed)
			void updateEdgeColors();
      /// Called when user fires "Resume" action, to clear downstream nodes from previous results
      void resetDownstream(TOPPASVertex* vertex);
			/// Runs the pipeline
			void runPipeline();
			/// Stores the pipeline to @p file
			void store(const String& file);
			/// Loads the pipeline from @p file
			void load(const String& file);
			/// Includes the pipeline @p scene
			void include(TOPPASScene* new_scene, QPointF pos = QPointF());
			/// Returns the file name
			const String& getSaveFileName();
			/// Sets the file name
			void setSaveFileName(const String& name);
			/// Performs a topological sort of all vertices
			void topoSort();
			/// Returns the name of the directory for output files
			const QString& getOutDir();
      /// Returns the name of the directory for temporary files
      const QString& getTempDir();
			/// Sets the name of the directory for output files
			void setOutDir(const QString& dir);
			/// Saves the pipeline if it has been changed since the last save.
			bool saveIfChanged();
			/// Sets the changed flag
			void setChanged(bool b);
			/// Returns if a pipeline is currently running
			bool isPipelineRunning();
			/// Shows a dialog that allows to specify the output directory. If @p always_ask == false, the dialog won't be shown if a directory has been set, already.
			bool askForOutputDir(bool always_ask = true);
			/// Enqueues the process, it will be run when the currently pending processes have finished
			void enqueueProcess(QProcess* p, const QString& command, const QStringList& args);
			/// Runs the next process in the queue, if any
			void runNextProcess();
			/// Resets the processes queue
			void resetProcessesQueue();
			/// Sets the clipboard content
			void setClipboard(TOPPASScene* clipboard);
			///Connects the signals to slots
			void connectVertexSignals(TOPPASVertex* tv);
			///Connects the signals to slots
			void connectToolVertexSignals(TOPPASToolVertex* ttv);
			///Connects the signals to slots
      void connectOutputVertexSignals(TOPPASOutputFileListVertex* oflv);
			///Connects the signals to slots
      void connectMergerVertexSignals(TOPPASMergerVertex* tmv);
			///Connects the signals to slots
			void connectEdgeSignals(TOPPASEdge* e);
			///Loads the @p resources into the input nodes of this workflow
			void loadResources(const TOPPASResources& resources);
			///Create @p resources from the current workflow
			void createResources(TOPPASResources& resources);
			///Returns whether the workflow has been changed since the latest "save"
			bool wasChanged();
			/// Refreshes the parameters of the TOPP tools in this workflow
			bool refreshParameters();
      /// determine dry run status (are tools actually called?)
      bool isDryRun() const;
			/// workflow description (to be displayed in TOPPAS window)
      QString getDescription() const;
      /// when description is updated by user, use this to update the description for later storage in file
      void setDescription(const QString& desc);
      /// sets the maximum number of jobs
      void setAllowedThreads(int num_threads);


		public slots:

			/// Terminates the currently running pipeline
			void abortPipeline();
			/// Called when an item is clicked
			void itemClicked();
			/// Called when an item is released
			void itemReleased();
			/// Called when the position of the hovering edge changes
			void updateHoveringEdgePos(const QPointF& new_pos);
			/// Called when a new out edge is supposed to be created
			void addHoveringEdge(const QPointF& pos);
			/// Called when the new edge is being "released"
			void finishHoveringEdge();
			/// Checks whether all output vertices are finished, and if yes, emits entirePipelineFinished() (called by finished output vertices)
			void checkIfWeAreDone();
			/// Called by vertices at which an error occurred during pipeline execution
			void pipelineErrorSlot(const QString msg = "");
			/// Moves all selected items by dx, dy
			void moveSelectedItems(qreal dx, qreal dy);
			/// Makes all vertices snap to the grid
			void snapToGrid();
			/// Sets if the running_ flag to true
			void setPipelineRunning(bool b = true);
			/// Called by a finished QProcess to indicate that we are free to start a new one
      void processFinished();
      /// dirty solution: when using ExecutePipeline this slot is called when the pipeline crashes. This will quit the app
      void quitWithError();


			///@name Slots for printing log/error output when no GUI is available
      //@{
      /// Writes the TOPP tool output to the logfile (and to stdout if no gui available)
      void logTOPPOutput(const QString& out);
      /// Writes the "tool started" message to the logfile (and to stdout if no gui available)
      void logToolStarted();
      /// Writes the "tool finished" message to the logfile (and to stdout if no gui available)
      void logToolFinished();
      /// Writes the "tool failed" message to the logfile (and to stdout if no gui available)
      void logToolFailed();
      /// Writes the "tool crashed" message to the logfile (and to stdout if no gui available)
      void logToolCrashed();
      /// Writes the "output file written" message to the logfile (and to stdout if no gui available)
      void logOutputFileWritten(const String& file);
			//@}
			
		signals:
			
			/// Emitted when the entire pipeline execution is finished
			void entirePipelineFinished();
			/// Emitted when the pipeline execution has failed
			void pipelineExecutionFailed();
			/// Emitted when the pipeline should be saved (showing a save as file dialog and so on)
			void saveMe();
			/// Kills all connected TOPP processes
			void terminateCurrentPipeline();
			/// Emitted when a selection is copied to the clipboard
			void selectionCopied(TOPPASScene* ts);
			/// Requests the clipboard content from TOPPASBase, will be stored in clipboard_
			void requestClipboardContent();
			/// Emitted when the main window needs to be updated
			void mainWindowNeedsUpdate();
      /// Emitted when files are triggered for opening in TOPPView
      void openInTOPPView(QStringList all_files);
      /// Emitted when in dry run mode and asked to run a TOPP tool (to fake success)
      void dryRunFinished(int, QProcess::ExitStatus);
      /// Emitted when there is an important message that needs to be printed in TOPPAS
      void messageReady(const QString& msg);


		protected:
			
			/// The current action mode
			ActionMode action_mode_;
			/// The list of all vertices
			VertexContainer vertices_;
			/// The list of all edges
			EdgeContainer edges_;
			/// The hovering edge which is currently being created
			TOPPASEdge* hover_edge_;
			/// The current potential target vertex of the hovering edge
			TOPPASVertex* potential_target_;
			/// The file name of this pipeline
			String file_name_;
			/// The path for temporary files
			QString tmp_path_;
			/// Are we in a GUI or is the scene used by ExecutePipeline (at the command line)?
			bool gui_;
			/// The directory where the output files will be written
			QString out_dir_;
			/// Flag that indicates if the pipeline has been changed since the last save
			bool changed_;
			/// Indicates if a pipeline is currently running
			bool running_;
      /// true if an error occurred during pipeline execution
      bool error_occured_;
			/// Indicates if the output directory has been specified by the user already
			bool user_specified_out_dir_;
			/// The queue of pending TOPP processes
			QList<TOPPProcess> topp_processes_queue_;
			/// Stores the clipboard content when requested from TOPPASBase
			TOPPASScene* clipboard_;
      /// dry run mode (no tools are actually called)
      bool dry_run_;
			/// currently running processes...
      int threads_active_;
      /// description text
      QString description_text_;
      /// maximum number of allowed threads
      int allowed_threads_;


			/// Returns the vertex in the foreground at position @p pos , if existent, otherwise 0.
			TOPPASVertex* getVertexAt_(const QPointF& pos);
			/// Returns whether an edge between node u and v would be allowed
			bool isEdgeAllowed_(TOPPASVertex* u, TOPPASVertex* v);
			/// DFS helper method. Returns true, if a back edge has been discovered
			bool dfsVisit_(TOPPASVertex* vertex);
			/// Performs a sanity check of the pipeline and notifies user when it finds something strange. Returns if pipeline OK.
      bool sanityCheck_();
			
			///@name reimplemented Qt events
      //@{
      void contextMenuEvent(QGraphicsSceneContextMenuEvent* event);
			//@}
			
			///Writes the @p text to the logfile
			void writeToLogFile_(const QString& text);
	};

}

#endif
