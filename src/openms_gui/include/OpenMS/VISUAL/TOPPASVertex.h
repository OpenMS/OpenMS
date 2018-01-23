// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_VISUAL_TOPPASVERTEX_H
#define OPENMS_VISUAL_TOPPASVERTEX_H

// ------------- DEBUGGING ----------------

// ---- Uncomment to enable debug mode ----
//#define TOPPAS_DEBUG
// ----------------------------------------

#ifdef TOPPAS_DEBUG
#define __DEBUG_BEGIN_METHOD__ \
  {                                                                                                                                                                                   \
    for (int dbg_indnt_cntr = 0; dbg_indnt_cntr < global_debug_indent_; ++dbg_indnt_cntr)   \
    {                                                                                                                                                                           \
      std::cout << "  ";                                                                                                                                  \
    }                                                                                                                                                                           \
    std::cout << "BEGIN [" << topo_nr_ << "] " << OPENMS_PRETTY_FUNCTION << std::endl;             \
    ++global_debug_indent_;                                                                                                                             \
  }

#define __DEBUG_END_METHOD__ \
  {                                                                                                                                                                                       \
    --global_debug_indent_;                                                                                                                             \
    if (global_debug_indent_ < 0) global_debug_indent_ = 0;                                                             \
    for (int dbg_indnt_cntr = 0; dbg_indnt_cntr < global_debug_indent_; ++dbg_indnt_cntr)   \
    {                                                                                                                                                                           \
      std::cout << "  ";                                                                                                                                  \
    }                                                                                                                                                                           \
    std::cout << "END [" << topo_nr_ << "] " << OPENMS_PRETTY_FUNCTION << std::endl;                   \
  }
#else
#define __DEBUG_BEGIN_METHOD__ {}
#define __DEBUG_END_METHOD__ {}
#endif

// ----------------------------------------

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtGui/QGraphicsSceneMouseEvent>
#include <QtGui/QGraphicsSceneContextMenuEvent>
#include <QtGui/QGraphicsItem>
#include <QtCore/QProcess>
#include <QtGui/QMenu>

namespace OpenMS
{
  class TOPPASEdge;

  /**
      @brief The base class of the different vertex classes.

      This class contains the common functionality (such as
      event handling for mouse clicks and drags) and holds the common
      members of all different kinds of vertices (e.g., containers
      for all in and out edges, the vertex ID, the number of a
      topological sort of the whole graph, etc.)

      @ingroup TOPPAS_elements
  */
  class OPENMS_GUI_DLLAPI TOPPASVertex :
    public QObject,
    public QGraphicsItem
  {
    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

public:

    /// The container for in/out edges
    typedef QList<TOPPASEdge *> EdgeContainer;
    /// A mutable iterator for in/out edges
    typedef EdgeContainer::iterator EdgeIterator;
    /// A const iterator for in/out edges
    typedef EdgeContainer::const_iterator ConstEdgeIterator;
	/// A class which interfaces with QStringList for holding filenames
	/// Incoming filenames are checked, and an exception is thrown if they are too long
	/// to avoid issues with common filesystems (due to filesystem limits).
	class TOPPASFilenames
	{
	  public:
		TOPPASFilenames()
		{
		}
	  
		int size() const;
		const QStringList& get() const;
		const QString& operator[](int i) const;

		///@name Setters; their all use check_() and can throw!
		//@{
		void set(const QStringList& filenames);
		void set(const QString& filename, int i);
		void push_back(const QString& filename);
		void append(const QStringList& filenames);
		//@}

	  private:
		/*
		@brief Check length of filename and throw Exception::FileNotWritable() if too long
		
		@param filename Full path to file (using relative paths will circumvent the effectiveness)
		@throw Exception::FileNotWritable() if too long (>=255 chars)
		*/
		void check_(const QString& filename);
		QStringList filenames_;   ///< filenames passed from upstream node in this round
	};
	/// Info for one edge and round, to be passed to next node
    struct VertexRoundPackage
    {
      VertexRoundPackage() :
        filenames(),
        edge(nullptr)
      {
      }

	  TOPPASFilenames filenames; ///< filenames passed from upstream node in this round
      TOPPASEdge* edge;  ///< edge that connects the upstream node to the current one
    };

	


    /// all infos to process one round for a vertex (from all incoming vertices)
    /// indexing via "parameter_index" of adjacent edge (could later be param_name) -> filenames
    /// Index for input and output edges is (-1) implicitly, thus we need signed type
    /// warning: the index refers to either input OR output (depending on if this structure is used for input files storage or output files storage)
    typedef std::map<Int, VertexRoundPackage> RoundPackage;
    typedef RoundPackage::const_iterator RoundPackageConstIt;
    typedef RoundPackage::iterator RoundPackageIt;

    /// all information a node needs to process all rounds
    typedef std::vector<RoundPackage> RoundPackages;

    /// The color of a vertex during depth-first search
    enum DFS_COLOR
    {
      DFS_WHITE,
      DFS_GRAY,
      DFS_BLACK
    };

    /// The color of a vertex during depth-first search
    enum SUBSTREESTATUS
    {
      TV_ALLFINISHED,  // all downstream nodes are done (including the ones which are feed by a parallel subtree)
      TV_UNFINISHED,   // some direct downstream node is not done
      TV_UNFINISHED_INBRANCH // a parallel subtree which merged with some downstream node A was not done (which prevented processing of the node A)
    };

    /// Default Constructor
    TOPPASVertex();
    /// Copy constructor
    TOPPASVertex(const TOPPASVertex & rhs);
    /// Destructor
    ~TOPPASVertex() override;
    /// Assignment operator
    TOPPASVertex & operator=(const TOPPASVertex & rhs);

    /// get the round package for this node from upstream
    /// -- indices in 'RoundPackage' mapping are thus referring to incoming edges of this node
    /// returns false on failure
    bool buildRoundPackages(RoundPackages & pkg, String & error_msg);

    /// check if all upstream nodes are ready to go ( 'finished_' is true)
    bool isUpstreamFinished() const;

    /// Returns the bounding rectangle of this item
    QRectF boundingRect() const override = 0;
    /// Returns a more precise shape
    QPainterPath shape() const override = 0;
    /// Paints the item
    void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget) override = 0;
    /// Returns begin() iterator of outgoing edges
    ConstEdgeIterator outEdgesBegin() const;
    /// Returns end() iterator of outgoing edges
    ConstEdgeIterator outEdgesEnd() const;
    /// Returns begin() iterator of incoming edges
    ConstEdgeIterator inEdgesBegin() const;
    /// Returns end() iterator of incoming edges
    ConstEdgeIterator inEdgesEnd() const;
    /// Returns the number of incoming edges
    Size incomingEdgesCount();
    /// Returns the number of outgoing edges
    Size outgoingEdgesCount();
    /// Adds an incoming edge
    void addInEdge(TOPPASEdge * edge);
    /// Adds an outgoing edge
    void addOutEdge(TOPPASEdge * edge);
    /// Removes an incoming edge
    void removeInEdge(TOPPASEdge * edge);
    /// Removes an outgoing edge
    void removeOutEdge(TOPPASEdge * edge);
    /// Returns the DFS color of this node
    DFS_COLOR getDFSColor();
    /// Sets the DFS color of this node
    void setDFSColor(DFS_COLOR color);
    /// Returns the DFS parent of this node
    TOPPASVertex * getDFSParent();
    /// Sets the DFS parent of this node
    void setDFSParent(TOPPASVertex * parent);
    /// Checks if all tools in the subtree below this node are finished
    TOPPASVertex::SUBSTREESTATUS getSubtreeStatus() const;
    /// Returns whether the vertex has been marked already (during topological sort)
    bool isTopoSortMarked();
    /// (Un)marks the vertex (during topological sort)
    void setTopoSortMarked(bool b);
    /// Returns the topological sort number
    UInt getTopoNr();
    /// Sets the topological sort number (overridden in tool and output vertices)
    virtual void setTopoNr(UInt nr);
    /// Resets the status
    /// @param reset_all_files Not used in this implementation, but in derived classes
    virtual void reset(bool reset_all_files = false);
    /// Marks this node (and everything further downstream) as unreachable. Overridden behavior in mergers.
    virtual void markUnreachable();
    /// Returns whether this node is reachable
    bool isReachable();
    /// Returns whether this node has already been processed during the current pipeline execution
    bool isFinished() const;
    /// run the tool (either ToolVertex, Merger, or OutputNode)
    /// @exception NotImplemented
    virtual void run();
    /// invert status of recycling
    virtual bool invertRecylingMode();
    /// get status of recycling
    bool isRecyclingEnabled() const;
    /// set status of recycling
    void setRecycling(const bool is_enabled);

    // get the name of the vertex (to be overridden by derived classes)
    virtual String getName() const = 0;

    /**
      @brief gets filenames for a certain output parameter (from this vertex), for a certain TOPPAS round


    */
    QStringList getFileNames(int param_index, int round) const;

    /// get all output files for all parameters for all rounds
    QStringList getFileNames() const;

    // get the output structure directly
    const RoundPackages & getOutputFiles() const;


    /// check if all upstream nodes are finished
    bool allInputsReady();


public slots:

    /// Called by an incoming edge when it has changed
    virtual void inEdgeHasChanged();
    /// Called by an outgoing edge when it has changed
    virtual void outEdgeHasChanged();

signals:

    /// Emitted when this item is clicked
    void clicked();
    /// Emitted when this item is released
    void released();
    /// Emitted when the position of the hovering edge changes
    void hoveringEdgePosChanged(const QPointF & new_pos);
    /// Emitted when a new out edge is supposed to be created
    void newHoveringEdge(const QPointF & pos);
    /// Emitted when the mouse is released after having dragged a new edge somewhere
    void finishHoveringEdge();
    /// Emitted when something has changed
    void somethingHasChanged();
    /// Emitted when the item is dragged
    void itemDragged(qreal dx, qreal dy);
    /// Emitted if an INI parameter or recycling mode or whatever was edited by the user - depending on the current state of the tool an action is taken
    /// each node type decides itself if this will invalidate the pipeline, depending on its internal status
    void parameterChanged(const bool invalidates_running_pipeline);

protected:

    /// The list of incoming edges
    EdgeContainer in_edges_;
    /// The list of outgoing edges
    EdgeContainer out_edges_;
    /// Indicates whether a new out edge is currently being created
    bool edge_being_created_;
    /// The color of the pen
    QColor pen_color_;
    /// The color of the brush
    QColor brush_color_;
    /// The DFS color of this node
    DFS_COLOR dfs_color_;
    /// The DFS parent of this node
    TOPPASVertex * dfs_parent_;
    /// "marked" flag for topological sort
    bool topo_sort_marked_;
    /// The number in a topological sort of the entire graph
    UInt topo_nr_;
    /// Stores the current output file names for each output parameter
    RoundPackages output_files_;
    /// number of rounds this node will do ('Merge All' nodes will pass everything, thus do only one round)
    int round_total_;
    /// currently finished number of rounds (TODO: do we need that?)
    int round_counter_;
    /// Stores whether this node has already been processed during the current pipeline execution
    bool finished_;
    /// Indicates whether this node is reachable (i.e. there is an input node somewhere further upstream)
    bool reachable_;
    /// shall subsequent tools be allowed to recycle the output of this node to match the number of rounds imposed by other parent nodes?
    bool allow_output_recycling_;


#ifdef TOPPAS_DEBUG
    // Indentation level for nicer debug output
    static int global_debug_indent_;
#endif

    ///@name reimplemented Qt events
    //@{
    void mouseReleaseEvent(QGraphicsSceneMouseEvent * e) override;
    void mousePressEvent(QGraphicsSceneMouseEvent * e) override;
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent * e) override;
    void mouseMoveEvent(QGraphicsSceneMouseEvent * e) override;
    void contextMenuEvent(QGraphicsSceneContextMenuEvent * event) override;
    //@}

    /// Moves the target pos of the edge which is just being created to @p pos
    virtual void moveNewEdgeTo_(const QPointF & pos);
    /// Returns a three character string (i.e. 001 instead of 1) for the given @p number
    String get3CharsNumber_(UInt number) const;

    /// Displays the debug output @p message, if TOPPAS_DEBUG is defined
    void debugOut_(const String &
#ifdef TOPPAS_DEBUG
                   message
#endif
                   ) const
    {
#ifdef TOPPAS_DEBUG
      for (int i = 0; i < global_debug_indent_; ++i)
      {
        std::cout << "  ";
      }
      std::cout << "[" << topo_nr_ << "] " << message << std::endl;
#endif
    }

  };
}

#endif
