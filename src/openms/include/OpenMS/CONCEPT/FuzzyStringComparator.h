// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h> 
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>


#include <map>
#include <sstream>

namespace OpenMS
{
  namespace Internal
  {
    namespace ClassTest
    {
      void OPENMS_DLLAPI testStringSimilar(const char * file,
                                           int line,
                                           const std::string & string_1,
                                           const char * string_1_stringified,
                                           const std::string & string_2,
                                           const char * string_2_stringified);

      bool OPENMS_DLLAPI isFileSimilar(const std::string &,
                                       const std::string &);
    }
  }

  /**
    @brief Fuzzy comparison of strings, tolerates numeric differences.
  */
  class OPENMS_DLLAPI FuzzyStringComparator
  {

    friend void OPENMS_DLLAPI
    Internal::ClassTest::testStringSimilar(
      const char * file,
      int line,
      const std::string & string_1,
      const char * string_1_stringified,
      const std::string & string_2,
      const char * string_2_stringified);

    friend bool OPENMS_DLLAPI
    Internal::ClassTest::isFileSimilar(const std::string &,
                                       const std::string &);

    /// %Internal exception class.
    struct AbortComparison
    {
    };

public:

    ///@name the fabulous four
    //@{

    /// Constructor
    FuzzyStringComparator();

    /// Destructor
    virtual
    ~FuzzyStringComparator();

    /// Copy constructor intentionally not implemented
    FuzzyStringComparator(const FuzzyStringComparator & rhs);

    /// Assignment operator intentionally not implemented
    FuzzyStringComparator & operator=(const FuzzyStringComparator & rhs);

    //@}

    /// Acceptable relative error (a number >= 1.0)
    const double & getAcceptableRelative() const;

    /// Acceptable relative error (a number >= 1.0)
    void setAcceptableRelative(const double rhs);

    /// Acceptable absolute difference (a number >= 0.0)
    const double & getAcceptableAbsolute() const;

    /// Acceptable absolute difference (a number >= 0.0)
    void setAcceptableAbsolute(const double rhs);

    /// White list.  If both lines contain the same element from this list, they are skipped over.
    const StringList & getWhitelist() const;

    /// White list.  If both lines contain the same element from this list, they are skipped over.
    StringList & getWhitelist();

    /// White list.  If both lines contain the same element from this list, they are skipped over.
    void setWhitelist(const StringList & rhs);

    /// Matched white list. If file 1 contains element 1 and file 2 contains element 2, they are skipped over.
    void setMatchedWhitelist(const std::vector< std::pair<std::string, std::string> >& rhs); 

    /// Matched white list. If file 1 contains element 1 and file 2 contains element 2, they are skipped over.
    const std::vector< std::pair<std::string, std::string> >& getMatchedWhitelist() const; 

    /**
      @brief verbose level

      - 0 = very quiet mode (absolutely no output)
      - 1 = quiet mode (no output unless differences detected)
      - 2 = default (include summary at end)
      - 3 = continue after errors
    */
    const int & getVerboseLevel() const;

    /**
      @brief verbose level

      - 0 = very quiet mode (absolutely no output)
      - 1 = quiet mode (no output unless differences detected)
      - 2 = default (include summary at end)
      - 3 = continue after errors
     */
    void setVerboseLevel(const int rhs);

    /**
      @brief get tab width (for column numbers)
    */
    const int & getTabWidth() const;

    /**
      @brief set tab width (for column numbers)
    */
    void setTabWidth(const int rhs);

    /**
      @brief get first column (for column numbers)
    */
    const int & getFirstColumn() const;

    /**
      @brief set first column (for column numbers)
    */
    void setFirstColumn(const int rhs);

    /**
      @brief Log output is written to this destination.

      The default is std::cout. Use std::ostringstream etc. to save the output
      in a string.
    */
    std::ostream & getLogDestination() const;

    /**
      @brief Log output is written to this destination.

      The default is std::cout. Use std::ostringstream etc. to save the output
      in a string.

      @internal There seems to be an issue with this under Windows, see comment
      in FuzzyStringComparator_test.cpp

    */
    void setLogDestination(std::ostream & rhs);

    /**
      @brief Compare two strings.

      This compares all lines of the input.

      @return true in case of no differences found
    */
    bool compareStrings(std::string const & lhs, std::string const & rhs);

    /**
      @brief Compare two streams of input.

      This compares all lines of the input.  Intended to be used for file
      streams.

      @return true in case of no differences found
    */
    bool compareStreams(std::istream & input_1, std::istream & input_2);

    /**
      @brief Simple diff-like application to compare two input files.
      Numeric differences are tolerated up to a certain ratio or absolute
      difference.

      where
      @param filename_1 first input file
      @param filename_2 second input file
      @return true in case of no differences found

      @sa ratio_max_allowed_
      @sa absdiff_max_allowed_
      @sa verbose_level_
    */
    bool compareFiles(const std::string & filename_1,
                      const std::string & filename_2);

protected:

    /**
      @brief Compare two lines of input.

      This implements the core functionality.  Intended to be used for a single
      line of input.

      @return true (non-zero) in case of success
     */
    bool compareLines_(std::string const & line_str_1, std::string const & line_str_2);

    /// Report good news.
    void reportSuccess_() const;

    /// Report bad news.
    /// @exception AbortComparison
    void reportFailure_(char const * const message) const;

    /// Write info about hits in the whitelist
    void writeWhitelistCases_(const std::string & prefix) const;

    /// read the next line in input stream, skipping over empty lines
    /// and lines consisting of whitespace only
    void readNextLine_(std::istream & input_stream, std::string & line_string, int & line_number) const;

    /// opens and checks an input file stream std::ifstream
    bool openInputFileStream_(const std::string & filename, std::ifstream & input_stream) const;

    /// Log and results output goes here
    std::ostream * log_dest_;

    /// Name of first input e.g., filename
    std::string input_1_name_;
    /// Name of second input e.g., filename
    std::string input_2_name_;

    /// Stores information about the current input line (i.e., stream for the line and the current position in the stream)
    struct InputLine
    {
      std::stringstream line_;
      std::ios::pos_type line_position_;

      InputLine();

      /// Initialize the input line to the passed string
      void setToString(const std::string & s);

      /// Save current position of the stream
      void updatePosition();

      /// Resets the stream to the last saved position
      void seekGToSavedPosition();

      /**
        @brief Convert to bool

        The function indicates success when none of the error flags (either failbit or badbit of the nested std::stringstream) are set.

        @return False on error, true otherwise.
      */
      bool ok() const;
    };

    InputLine input_line_1_;
    InputLine input_line_2_;

    int line_num_1_;
    int line_num_2_;

    int line_num_1_max_;
    int line_num_2_max_;

    std::string line_str_1_max_;
    std::string line_str_2_max_;

    /// Maximum ratio of numbers allowed, see @em ratio_max_.
    double ratio_max_allowed_;

    /// Maximum ratio of numbers observed so far, see @em ratio_max_allowed_.
    double ratio_max_;

    /// Maximum absolute difference of numbers allowed, see @em absdiff_max_.
    double absdiff_max_allowed_;

    /// Maximum difference of numbers observed so far, see @em absdiff_max_allowed_.
    double absdiff_max_;

    /// Stores information about characters, numbers, and white spaces loaded from the InputStream
    struct StreamElement_
    {
      double number;
      unsigned char letter;
      bool is_number;
      bool is_space;

      StreamElement_();

      /// reset all elements of the element to default value
      void reset();

      /// Read the next element from an InputLine and update the InputLine accordingly
      /// The @p str_line contains the same data as the stream, since it saves some forth-and-back conversion internally
      /// TODO: avoid streams all together (slow, and no random access, required by boost::qi) at some point
      void fillFromInputLine(InputLine& input_line, const std::string& str_line);
    };

    /// Stores information about characters, numbers, and white spaces loaded from the first input stream
    StreamElement_ element_1_;
    /// Stores information about characters, numbers, and white spaces loaded from the second input stream
    StreamElement_ element_2_;

    /// Wrapper for the prefix information computed for the failure report
    struct PrefixInfo_
    {
      OpenMS::String prefix;
      OpenMS::String prefix_whitespaces;
      int line_column;

      PrefixInfo_(const InputLine & input_line, const int tab_width_, const int first_column_);
    };

    bool is_absdiff_small_;

    int verbose_level_;
    int tab_width_;
    int first_column_;

    /**
      @brief Has comparison been successful so far?  Note: this flag is
             changed in reportFailure_();
     */
    bool is_status_success_;

    /// use a prefix when reporting
    bool use_prefix_;

    /// Whitelist
    StringList whitelist_;
    /// Occurrences of whitelist entries
    std::map<String, UInt> whitelist_cases_;

    /// Alternative Whitelist
    std::vector< std::pair<std::string, std::string> > matched_whitelist_; 
  }; // class FuzzyStringComparator

} //namespace OpenMS

