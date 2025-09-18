// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <map>
#include <iosfwd>

namespace OpenMS
{
  /**
   * @brief Provides a central class to register globally used output streams. Currently supported streams are
   *
   * <ul>
   *  <li>std::ofstream</li>
   *  <li>std::ostringstream</li>
   * </ul>
   *
   * You can simply request a stream from StreamHandler and it will manage the process of construction and
   * destruction (if the stream is not needed any more).
   *
   * A normal scenario would be
   *
   * <code>
   * STREAM_HANDLER.registerStream(StreamHandler::FILE, "name_of_the_output_file");<br>
   * STREAM_HANDLER.getStream(StreamHandler::FILE, "name_of_the_output_file") << "This will be send to the file" << std::endl;<br>
   * // ...<br>
   * <br>
   * // you do not need the file any more<br>
   * STREAM_HANDLER.unregisterStream(StreamHandler::FILE, "name_of_the_output_file");<br>
   * // StreamHandler will close the file stream if no other part of your program requested the same stream
   * </code>
   *
   * @ingroup Concept
   */
  class OPENMS_DLLAPI StreamHandler
  {
public:
    /**
     * @brief Defines the type of the stream that should be handled
     */
    enum StreamType
    {
      FILE,
      STRING
    };

    /// Default constructor
    StreamHandler();

    /// destructor
    virtual ~StreamHandler();

    /**
     * @brief Creates a stream of type @p type and with name @p stream_name
     *
     * If the stream is already registered the reference counter is increased.
     *
     * @note The stream name must be unique. You cannot register the same stream name with two different types.
     *
     * @param type  Type of the stream (e.g. FILE)
     * @param stream_name Name of the stream (e.g. the file name for a file stream).
     *
     * @return An integer indicating if the operation was completed successfully (@p value != 1 means a failure occurred).
     */
    Int registerStream(StreamType const type, const String & stream_name);

    /**
     * @brief De-registers a stream of type @p type and with name @p stream_name from the handler.
     *
     * It also decreases the reference counter for the named stream. If the counter
     * reaches 0. The stream will be closed.
     *
     * @param type  Type of the stream (e.g. FILE)
     * @param stream_name Name of the stream (e.g. the file name for a file stream).
     *
     */
    void unregisterStream(StreamType const type, const String & stream_name);

    /**
     * @brief Returns a reference to the stream of type @p type and with name @p stream_name.
     *
     * If the stream was not registered before an ElementNotFoundException will be thrown.
     *
     * @param type  Type of the stream (e.g. FILE)
     * @param stream_name Name of the stream (e.g. the file name for a file stream).
     *
     * @throw ElementNotFoundException
     *
     * @return A reference to the requested stream.
     */
    std::ostream & getStream(StreamType const type, const String & stream_name);

    /**
     * @brief Returns true if the stream @p stream_name with type @p type is
     * registered.
     *
     * @param type  Type of the stream (e.g. FILE)
     * @param stream_name Name of the stream (e.g. the file name for a file stream).
     *
     * @return bool indication if the stream is known.
     */
    bool hasStream(const StreamType type, const String & stream_name);

protected:

    std::map<String, std::ostream *>  name_to_stream_map_;  ///< Maps all registered stream names to the corresponding std::ostream.
    std::map<String, StreamType> name_to_type_map_;  ///< Maps all registered stream names to the corresponding StreamHandler::StreamType
    std::map<String, Size>      name_to_counter_map_;   ///< Maps all registered stream names to the number of times it was registered. If the counter goes to zero, the stream will be closed and removed.

    /**
     * @brief Creates a stream with the given type and the given name.
     *
     * @param type  Type of the stream (e.g. FILE)
     * @param stream_name Name of the stream (e.g. the file name for a file stream).
     *
     * @return A pointer to the created stream.
     */
    std::ostream * createStream_(const StreamType type, const String & stream_name);

private:
    // copy constructor and assignment operator are hidden to avoid
    // creating multiple pointers to a single filestream instance

    /// copy constructor
    StreamHandler(const StreamHandler & source);

    /// assignment operator
    virtual StreamHandler & operator=(const StreamHandler & source);

    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, StreamHandler const & stream_handler);
  };

  /// Overload for the \a insertion \a operator (operator<<) to have a formatted output of the StreamHandler
  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, StreamHandler const & stream_handler);

  /// Global StreamHandler instance.
  OPENMS_DLLAPI extern StreamHandler STREAM_HANDLER;
} // end namespace OpenMS

