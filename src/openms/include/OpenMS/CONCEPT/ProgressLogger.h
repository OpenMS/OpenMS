// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
  class String;

  /**
    @brief Base class for all classes that want to report their progress.

    Per default the progress log is disabled. Use setLogType to enable it.

    Use startProgress, setProgress and endProgress for the actual logging.

    @note All methods are const, so it can be used through a const reference or in const methods as well!
  */
  class OPENMS_DLLAPI ProgressLogger
  {
public:
    /// Constructor
    ProgressLogger();

    /// Destructor
    virtual ~ProgressLogger();

    /// Copy constructor
    ProgressLogger(const ProgressLogger& other);

    /// Assignment Operator
    ProgressLogger& operator=(const ProgressLogger& other);

    /// Possible log types
    enum LogType
    {
      CMD, ///< Command line progress
      GUI, ///< Progress dialog
      NONE ///< No progress logging
    };

    /**
      @brief This class represents an actual implementation of a logger.
    */
    class OPENMS_DLLAPI ProgressLoggerImpl
    {
public:
      virtual void startProgress(const SignedSize begin, const SignedSize end, const String& label, const int current_recursion_depth) const = 0;
      virtual void setProgress(const SignedSize value, const int current_recursion_depth) const = 0;
      virtual SignedSize nextProgress() const = 0; //< does not print/show anything; returns current progress
      /// finalize; usually stops the clock and prints a summary; 
      /// You may optionally pass the number of bytes read/written to get a MB/sec estimate
      virtual void endProgress(const int current_recursion_depth, UInt64 bytes_processed = 0) const = 0;

      virtual ~ProgressLoggerImpl() {}

      /// Factory requirements
  

    };

    /// Sets the progress log that should be used. The default type is NONE!
    void setLogType(LogType type) const;

    /// Returns the type of progress log being used.
    LogType getLogType() const;

    /// @brief  Sets the logger to be used for progress logging
    /// @param logger 
    void setLogger(ProgressLoggerImpl* logger);

    /**
      @brief Initializes the progress display

      Sets the progress range from @p begin to @p end.
      If @p begin equals @p end, setProgress only indicates that
      the program is still running, but without showing any absolute progress value.

      Sets the label to @p label.

      @note Make sure to call setLogType first!
    */
    void startProgress(SignedSize begin, SignedSize end, const String& label) const;

    /// Sets the current progress
    void setProgress(SignedSize value) const;

    /// Ends the progress display
    /// You may optionally pass the number of bytes read/written to get a MB/sec estimate (if -1 or 0, it will be ignored)
    void endProgress(UInt64 bytes_processed = 0) const;

    /// increment progress by 1 (according to range begin-end)
    void nextProgress() const;

protected:
    mutable LogType type_;
    mutable time_t last_invoke_;
    static int recursion_depth_;

    mutable ProgressLoggerImpl* current_logger_;

  };

  // Function pointer for injecting the GUI progress logger implementation
  typedef ProgressLogger::ProgressLoggerImpl* (*MakeGUIProgressLoggerFunc)();
  extern OPENMS_DLLAPI MakeGUIProgressLoggerFunc make_gui_progress_logger;

} // namespace OpenMS

