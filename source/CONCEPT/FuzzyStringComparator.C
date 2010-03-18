// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>
#include <sstream>

namespace OpenMS
{

  FuzzyStringComparator::FuzzyStringComparator()
    :
    log_dest_(&std::cout),
    input_1_name_("input_1"),
    input_2_name_("input_2"),
    line_1_(),
    line_2_(),
    ratio_max_allowed_(1.0),
    ratio_max_(1.0),
    absdiff_max_allowed_(0.0),
    absdiff_max_(0.0),
    number_1_(0),
    letter_1_(0),
    is_number_1_(false),
    is_space_1_(false),
    number_2_(0),
    letter_2_(0),
    is_number_2_(false),
    is_space_2_(false),
    is_absdiff_small_(false),
    line_num_1_(0),
    line_num_2_(0),
    line_num_1_max_(-1),
    line_num_2_max_(-1),
    verbose_level_(2),
		tab_width_(8),
		first_column_(1),
    is_status_success_(true),
    line_str_1_max_(),
    line_str_2_max_(),
		use_prefix_(false),
		whitelist_(),
		whitelist_cases_()
  {
	}

  FuzzyStringComparator::~FuzzyStringComparator(){}

  void FuzzyStringComparator::reportFailure_( char const * const message ) const
  {
    // We neither want this entire method be non-const nor make
    // is_status_success_ a mutable.  So lets hack around it.  (Documented in
    // class.)
    const_cast<bool&>(is_status_success_) = false;

    if ( verbose_level_ >= 1 )
    {
 			int line_1_col = 0;
			OpenMS::String pre1(line_1_.str());
			pre1 = pre1.prefix(size_t(line_1_pos_));
			OpenMS::String pre1_white(pre1);
			for ( String::iterator iter = pre1_white.begin(); iter != pre1_white.end(); ++iter )
			{
				if ( *iter != '\t' )
				{
					*iter = ' ';
					++line_1_col;
				}
				else
				{
					line_1_col = (line_1_col/tab_width_+1)*tab_width_;
				}
			}
			line_1_col += first_column_;

			int line_2_col = 0;
			OpenMS::String pre2(line_2_.str());
			pre2 = pre2.prefix(size_t(line_2_pos_));
			OpenMS::String pre2_white(pre2);
			for ( String::iterator iter = pre2_white.begin(); iter != pre2_white.end(); ++iter )
			{
				if ( *iter != '\t' )
				{
					*iter = ' ';
					++line_2_col;
				}
				else
				{
					line_2_col = (line_2_col/tab_width_+1)*tab_width_;
				}
			}
			line_2_col += first_column_;

			std::string prefix;
			if ( use_prefix_ )
			{
				prefix = "   :|:  ";
			}

			*log_dest_ << std::boolalpha <<
				prefix << "FAILED: '" << message << "'\n" <<
				prefix << "\n" <<
				prefix << "  input:\tin1\tin2\n" <<
				prefix << "  line:\t" << line_num_1_ << '\t' << line_num_2_ << "\n" <<
				prefix << "  pos/col:\t" << line_1_pos_ << '/' << line_1_col << '\t' << line_2_pos_ << '/' << line_2_col << "\n" <<
				prefix << " --------------------------------\n" <<
				prefix << "  is_number:\t" << is_number_1_ << '\t' << is_number_2_ << "\n" <<
				prefix << "  numbers:\t" << number_1_ << '\t' << number_2_ << "\n" <<
				prefix << "  is_space:\t" << is_space_1_ << '\t' << is_space_2_ << "\n" <<
				prefix << "  is_letter:\t" << (!is_number_1_&&!is_space_1_) << '\t' << (!is_number_2_&&!is_space_2_) << "\n" <<
				prefix << "  letters:\t\"" << letter_1_ << "\"\t\"" << letter_2_ << "\"\n" <<
				prefix << "  char_codes:\t" << static_cast<UInt>(letter_1_) << "\t" << static_cast<UInt>(letter_2_) << "\n" <<
				prefix << " --------------------------------\n" <<
				prefix << "  relative_max:        " << ratio_max_ << "\n" <<
				prefix << "  relative_acceptable: " << ratio_max_allowed_ << "\n" <<
				prefix << " --------------------------------\n" <<
				prefix << "  absolute_max:        " << absdiff_max_ << "\n" <<
				prefix << "  absolute_acceptable: " << absdiff_max_allowed_ << std::endl;

			writeWhitelistCases_(prefix);

			*log_dest_
				<< prefix << "\n"
				<< prefix << "Offending lines:\t\t\t(tab_width = " << tab_width_ << ", first_column = " << first_column_ << ")\n"
				<< prefix << "\n"
				<< prefix << "in1:  " << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() << "   (line: " << line_num_1_ << ", position/column: " << line_1_pos_ << '/' << line_1_col << ")\n"
				<< prefix << pre1 << "!\n"
				<< prefix << pre1_white << OpenMS::String(line_1_.str()).suffix(line_1_.str().size()-pre1.size()) << "\n"
				<< prefix <<  "\n"
				<< prefix << "in2:  " << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() << "   (line: " << line_num_2_ << ", position/column: " << line_2_pos_ << '/' << line_2_col << ")\n"
				<< prefix << pre2 << "!\n"
				<< prefix << pre2_white << OpenMS::String(line_2_.str()).suffix(line_2_.str().size()-pre2.size()) << "\n"
				<< prefix << "\n\n"
        << "Easy Access:" << "\n"
				<< QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() << ':' << line_num_1_ << ":" << line_1_col << ":\n"
				<< QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() << ':' << line_num_2_ << ":" << line_2_col << ":\n"
				<< "\n"
#ifdef WIN32
        << "TortoiseMerge"
        << " /base:\"" << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() << "\""
        << " /mine:\"" << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() << "\"" 
#else
        << "diff"
        << " " << QDir::toNativeSeparators(File::absolutePath(input_1_name_).toQString()).toStdString() 
        << " " << QDir::toNativeSeparators(File::absolutePath(input_2_name_).toQString()).toStdString() 
#endif
        << std::endl;
		}

		// If verbose level is low, report only the first error.
    if ( verbose_level_ < 3 )
    {
      throw FuzzyStringComparator::AbortComparison();
    }

    return;
  } // reportFailure_()

  void FuzzyStringComparator::reportSuccess_() const
  {
    if (  is_status_success_ && verbose_level_ >= 2 )
    {
			std::string prefix;
			if ( use_prefix_ )
			{
				prefix = "   :|:  ";
			}

      *log_dest_ <<
				prefix << "PASSED.\n" <<
				prefix << '\n' <<
				prefix << "  relative_max:        " << ratio_max_ << '\n' <<
				prefix << "  relative_acceptable: " << ratio_max_allowed_ << '\n' <<
				prefix << '\n' <<
				prefix << "  absolute_max:        " << absdiff_max_ << '\n' <<
				prefix << "  absolute_acceptable: " << absdiff_max_allowed_ << std::endl;

			writeWhitelistCases_(prefix);

			*log_dest_ << prefix << std::endl;

      if ( line_num_1_max_ == -1 && line_num_2_max_ == -1 )
      {
				*log_dest_ <<
					prefix << "No numeric differences were found.\n" <<
					prefix << std::endl;
      }
      else
      {
				*log_dest_ <<
					prefix << "Maximum relative error was attained at these lines, enclosed in \"\":\n" <<
					prefix << '\n' <<
					QDir::toNativeSeparators(input_1_name_.c_str()).toStdString() << ':' << line_num_1_max_ << ":\n" <<
					"\""<< line_str_1_max_ << "\"\n" <<
					'\n' <<
					QDir::toNativeSeparators(input_2_name_.c_str()).toStdString() << ':' << line_num_2_max_ << ":\n" <<
					"\""<< line_str_2_max_ << "\"\n" <<
					std::endl;
      }
    }
    return;
  }

  bool FuzzyStringComparator::compareLines_( std::string const & line_str_1, std::string const & line_str_2 )
  {

		for ( StringList::const_iterator slit = whitelist_.begin();
					slit != whitelist_.end();
					++slit
				)
		{
			if ( line_str_1.find(*slit)!=std::string::npos &&
					 line_str_2.find(*slit)!=std::string::npos
				 )
			{
				++whitelist_cases_[*slit];
				// *log_dest_ << "whitelist_ case: " << *slit << '\n';
				return is_status_success_;
			}
		}

    line_1_.str(line_str_1);
    line_1_.seekp(0);
    line_1_.clear();
		line_1_.unsetf(std::ios::skipws);

    line_2_.str(line_str_2);
    line_2_.seekp(0);
    line_2_.clear();
		line_2_.unsetf(std::ios::skipws);

    try
    {
      while ( line_1_ && line_2_ )
      {
				is_number_1_ = false;
				is_number_2_ = false;
				is_space_1_ = false;
				is_space_2_ = false;
				letter_1_ = '\0';
				letter_2_ = '\0';
				number_1_ = std::numeric_limits<double>::quiet_NaN();
				number_2_ = std::numeric_limits<double>::quiet_NaN();

				line_1_pos_ = line_1_.tellg(); // save current reading position
				line_1_ >> letter_1_; // read letter
				// std::cout << ":::" << letter_1_ << line_1_pos_ << std::endl;
				if ( ( is_space_1_ = (isspace(letter_1_)!=0) ) ) // is whitespace?
				{
					line_1_ >> std::ws; // skip over further whitespace
				}
				else
				{
					line_1_.seekg(line_1_pos_); // rewind to saved position
					if ( ( is_number_1_ = ( ( line_1_ >> number_1_ )!=0) ) ) // is a number?
					{
						// letter_1_ = '\0';
						// std::cout << line_1_pos_ << std::endl;
					}
					else
					{
						line_1_.clear(); // reset status
						line_1_.seekg(line_1_pos_); // rewind to saved position
						line_1_ >> letter_1_; // read letter
					}
				}

				line_2_pos_ = line_2_.tellg(); // save current reading position
				line_2_ >> letter_2_; // read letter
				if ( ( is_space_2_ = ( isspace(letter_2_)!=0 ) ) ) // is whitespace?
				{
					line_2_ >> std::ws; // skip over further whitespace
				}
				else
				{
					line_2_.seekg(line_2_pos_); // rewind to saved position
					if ( ( is_number_2_ = ( ( line_2_ >> number_2_ )!=0) ) ) // is a number?
					{
						// letter_2_ = '\0';
					}
					else
					{
						line_2_.clear(); // reset status
						line_2_.seekg(line_2_pos_); // rewind to saved position
						line_2_ >> letter_2_; // read letter
					}
				}


				if ( is_number_1_ )
				{
					if ( is_number_2_ )
					{ // we are comparing numbers

						// check if absolute difference is small
						double absdiff = number_1_ - number_2_;
						if ( absdiff < 0 )
						{
							absdiff = -absdiff;
						}
						if ( absdiff > absdiff_max_ )
						{
							absdiff_max_ = absdiff;
						}
						// If absolute difference is small, large relative errors will be
						// tolerated in the cases below.  But a large absolute difference is
						// not an error, if relative error is small.  We do not jump out of
						// the case distinction here because we want to record the relative
						// error even in case of a successful comparison.
						is_absdiff_small_ = ( absdiff <= absdiff_max_allowed_ );

						if ( !number_1_ )
						{ // number_1_ is zero
							if (!number_2_ )
							{ // both numbers are zero
								continue;
							}
							else
							{
								if ( !is_absdiff_small_ )
								{
									reportFailure_("number_1_ is zero, but number_2_ is not");
									continue;
								}
							}
						}
						else
						{ // number_1_ is not zero
							if ( !number_2_ )
							{
								if ( !is_absdiff_small_ )
								{
									reportFailure_("number_1_ is not zero, but number_2_ is");
									continue;
								}
							}
							else
							{ // both numbers are not zero
								double ratio = number_1_ / number_2_;
								if ( ratio < 0 )
								{
									if ( !is_absdiff_small_ )
									{
										reportFailure_("numbers have different signs");
										continue;
									}
								}
								else
								{ // ok, numbers have same sign, but we still need to check their ratio
									if ( ratio < 1 )
									{ // take reciprocal value
										ratio = 1.0 / ratio;
									}
									// by now, we are sure that ratio >= 1
									if ( ratio > ratio_max_ )
									{ // update running max
										ratio_max_ = ratio;
										line_num_1_max_ = line_num_1_;
										line_num_2_max_ = line_num_2_;
										line_str_1_max_ = line_str_1;
										line_str_2_max_ = line_str_2;
										if ( ratio_max_ > ratio_max_allowed_ )
										{
											if ( !is_absdiff_small_ )
											{
												reportFailure_("ratio of numbers is too large");
												continue;
											}
										}
									}
								}
								// okay
								continue;
							}
						}
					}
					else
					{
						reportFailure_("input_1 is a number, but input_2 is not");
						continue;
					}
				}
				else
				{ // input_1 is not a number
					if ( is_number_2_ )
					{
						reportFailure_("input_1 is not a number, but input_2 is");
						continue;
					}
					else
					{ // ok, both inputs are not numbers, let us compare them as characters or whitespace
						if ( is_space_1_ )
						{
							if ( is_space_2_ )
							{ // ok, both inputs are whitespace
								continue;
							}
							else
							{
                if ( letter_1_ == ASCII__CARRIAGE_RETURN ) // should be 13 == ascii carriage return char
                {
									// we skip over '\r'
									line_2_.clear(); // reset status
									line_2_.seekg(line_2_pos_); // rewind to saved position
                  continue;
									//reportFailure_("input_1 is carriage return, but input_2_ is not whitespace");
                }
                else
                {
                  reportFailure_("input_1 is whitespace, but input_2 is not");
                }
								continue;
							}
						}
						else
						{ // input_1 is not whitespace
							if ( is_space_2_ )
							{
							  if ( letter_2_ == ASCII__CARRIAGE_RETURN ) // should be 13 == ascii carriage return char
							  {
									// we skip over '\r'
									line_1_.clear(); // reset status
									line_1_.seekg(line_1_pos_); // rewind to saved position
									continue;
							    //reportFailure_("input_1 is not whitespace, but input_2 is carriage return");
							  }
							  else
                {
							    reportFailure_("input_1 is not whitespace, but input_2 is");
                }
								continue;
							}
							else
							{ // both inputs are neither numbers nor whitespace, let us compare them as characters
								if ( letter_1_ == letter_2_ )
								{ // ok, same characters
									continue;
								}
								else
								{
									reportFailure_("different letters");
									continue;
								}
							}
						}
					}
				}

				if ( is_absdiff_small_ )
				{
					is_absdiff_small_ = false;
					continue;
				}

				verbose_level_ = 10000;
				reportFailure_
					("This cannot happen.  You should never get here ... "
					 "please report this bug along with the data that produced it."
					);

      } // while ( line_1_ || line_2_ )
			if ( line_1_ && !line_2_ )
			{
				reportFailure_("line from input_2 is shorter than line from input_1");
			}
			if ( !line_1_ && line_2_ )
			{
				reportFailure_("line from input_1 is shorter than line from input_2");
			}
    }
    catch ( FuzzyStringComparator::AbortComparison const& )
    {
      // *log_dest_ << "compareLines_(): Caught FuzzyStringComparator::AbortComparison\n";
    }

    return is_status_success_;
  } // compareLines_()

  bool FuzzyStringComparator::compareStrings( std::string const & lhs, std::string const & rhs )
  {
		std::istringstream input_1(lhs);
		std::istringstream input_2(rhs);

    std::string line_str_1;
    std::string line_str_2;

    while ( input_1 || input_2 )
    {

			// read the next line in both input streams, skipping over empty lines
			// and lines consisting of whitespace only

			for ( line_str_1.clear(); ++line_num_1_, std::getline(input_1,line_str_1); )
			{
				if ( line_str_1.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_1.begin(); // loop initialization
				for ( ; iter != line_str_1.end() && isspace((unsigned char)*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_1.end() ) break; // line is not empty or whitespace only
			}

			for ( line_str_2.clear(); ++line_num_2_, std::getline(input_2,line_str_2); )
			{
				if ( line_str_2.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_2.begin(); // loop initialization
				for ( ; iter != line_str_2.end() && isspace((unsigned char)*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_2.end() ) break; // line is not empty or whitespace only
			}

			// compare the two lines of input
			if ( !compareLines_(line_str_1, line_str_2) && verbose_level_ < 3 ) break;

    } // while ( input_1 || input_2 )

		reportSuccess_();

		return is_status_success_;

  } // compareStrings()

  bool FuzzyStringComparator::compareStreams( std::istream & input_1, std::istream & input_2 )
  {
    std::string line_str_1;
    std::string line_str_2;

    while ( input_1 || input_2 )
    {

			// read the next line in both input streams, skipping over
			// - empty lines
			// - lines consisting of whitespace only

			for ( line_str_1.clear(); ++line_num_1_, std::getline(input_1,line_str_1); )
			{
				if ( line_str_1.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_1.begin(); // loop initialization
				for ( ; iter != line_str_1.end() && isspace((unsigned char)*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_1.end() ) break; // line is not empty or whitespace only
			}
			//std::cout << "eof: " << input_1.eof() << " failbit: " << input_1.fail() << " badbit: " << input_1.bad() << " reading " << input_1.tellg () << "chars\n";

			for ( line_str_2.clear(); ++line_num_2_, std::getline(input_2,line_str_2); )
			{
				if ( line_str_2.empty() ) continue; // shortcut
				std::string::const_iterator iter = line_str_2.begin(); // loop initialization
				for ( ; iter != line_str_2.end() && isspace((unsigned char)*iter); ++iter ) ; // skip over whitespace
				if ( iter != line_str_2.end() ) break; // line is not empty or whitespace only
			}
			//std::cout << "eof: " << input_2.eof() << " failbit: " << input_2.fail() << " badbit: " << input_2.bad() << " reading " << input_2.tellg () << "chars\n";

			// compare the two lines of input
			if ( !compareLines_(line_str_1, line_str_2) && verbose_level_ < 3 ) break;

    } // while ( input_1 || input_2 )

		reportSuccess_();

    return is_status_success_;

  } // compareStreams()

	bool FuzzyStringComparator::compareFiles(const std::string & filename_1, const std::string & filename_2)
  {

		input_1_name_ = filename_1;
		input_2_name_ = filename_2;

		if ( input_1_name_ == input_2_name_ )
		{
			*log_dest_ << "Error: first and second input file have the same name.  That's cheating!\n";
			return false;
		}

		std::ifstream  input_1_f;
		input_1_f.open(input_1_name_.c_str(), std::ios::in | std::ios::binary);
		if ( !input_1_f )
		{
			*log_dest_ << "Error opening first input file '" << input_1_name_ <<"'.\n";
			return false;
		}
		input_1_f.unsetf(std::ios::skipws);

		std::ifstream  input_2_f;
		input_2_f.open(input_2_name_.c_str(), std::ios::in | std::ios::binary);
		if ( !input_2_f )
		{
			*log_dest_ << "Error opening second input file '" << input_2_name_ <<"'.\n";
			return false;
		}
		input_2_f.unsetf(std::ios::skipws);

		//------------------------------------------------------------
		// main loop

		compareStreams(input_1_f, input_2_f);

		return is_status_success_;

	} // compareFiles()


	void FuzzyStringComparator::writeWhitelistCases_(const std::string& prefix) const
	{
		if (!whitelist_cases_.empty())
		{
			*log_dest_ <<
				prefix << '\n' <<
				prefix << "  whitelist cases:\n";
			Size length = 0;
			for (std::map<String,UInt>::const_iterator wlcit = 
						 whitelist_cases_.begin(); wlcit != whitelist_cases_.end(); 
					 ++wlcit)
			{
				if (wlcit->first.size() > length) length = wlcit->first.size();
			}
			for (std::map<String, UInt>::const_iterator wlcit = 
						 whitelist_cases_.begin(); wlcit != whitelist_cases_.end();
					 ++wlcit)
			{
				*log_dest_ << 
					prefix << "    " << std::setw(length+3) << std::left << 
					( "\"" + wlcit->first + "\"" ) << std::setw(3) << std::right << 
					wlcit->second << "x\n";
			}
		}
	}


} //namespace OpenMS
