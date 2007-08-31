// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS 
{

	XTandemInfile::XTandemInfile()
		: peak_mass_tolerance_(0.4),
			precursor_mass_tolerance_plus_(100.0),
			precursor_mass_tolerance_minus_(100.0),
      precursor_monoisotopic_error_(XTandemInfile::MONOISOTOPIC),
      precursor_mass_error_unit_(XTandemInfile::PPM),
      peak_mass_error_unit_(XTandemInfile::DALTONS),
      dynamic_range_(100.0),
      total_number_peaks_(50),
      max_precursor_charge_(4),
      noise_supression_(true),
			precursor_lower_mz_(500.0),
			peak_lower_mz_(150.0),
			min_number_peaks_(15),
			number_of_threads_(1),
      batch_size_(1000),
			fixed_modifications_(""),
			variable_modifications_(""),
			variable_modification_motif_(""),
      input_filename_(""),
			output_filename_("")

	{
	  	
	}
	
	XTandemInfile::~XTandemInfile()
	{
	}
	
  void XTandemInfile::load(const String& filename) throw (Exception::FileNotFound, Exception::ParseError)
  {
  	//try to open file
		if (!File::exists(filename))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
		
		// initialize parser
		try 
		{
			XMLPlatformUtils::Initialize();
		}
		catch (const XMLException& toCatch) 
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + Internal::StringManager().convert(toCatch.getMessage()) );
	  }

		SAX2XMLReader* parser = XMLReaderFactory::createXMLReader();
		parser->setFeature(XMLUni::fgSAX2CoreNameSpaces,false);
		parser->setFeature(XMLUni::fgSAX2CoreNameSpacePrefixes,false);

		Internal::XTandemInfileXMLHandler handler(filename, this);


		parser->setContentHandler(&handler);
		parser->setErrorHandler(&handler);
		
		LocalFileInputSource source( Internal::StringManager().convert(filename.c_str()) );
		try 
    {
    	parser->parse(source);
    	delete(parser);
    }
    catch (const XMLException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + Internal::StringManager().convert(toCatch.getMessage()) );
    }
    catch (const SAXException& toCatch) 
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + Internal::StringManager().convert(toCatch.getMessage()) );
    }

	}
		
	void XTandemInfile::write(const String& filename) throw (Exception::UnableToCreateFile)
	{
		if (!File::writable(filename))
		{
			throw (Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename));
		}
		ofstream os(filename.c_str());
		writeTo_(os);
		return;
	}

	void XTandemInfile::writeTo_(ostream& os)
	{
		os << "<?xml version=\"1.0\"?>" << endl
			 << "<?xml-stylesheet type=\"text/xsl\" href=\"tandem-input-style.xsl\"?>" << endl
			 << "<bioml>" << endl;

		writeNote_(os, "input", "list path, default parameters", "/home/andreas/DATA/MSSoftware/tandem-linux-07-04-01-1/bin/default_input.xml");
		writeNote_(os, "input", "list path, taxonomy information", "/home/andreas/DATA/MSSoftware/tandem-linux-07-04-01-1/bin/taxonomy.xml");

		//writeNote_(os, "input", "spectrum, fragment monoisotopic mass error", String(peak_mass_tolerance_));
		//writeNote_(os, "input", "spectrum, parent monoisotopic mass error plus", String(precursor_mass_tolerance_plus_)); 
		//writeNote_(os, "input", "spectrum, parent monoisotopic mass error minus", String(precursor_mass_tolerance_minus_));
		writeNote_(os, "input", "output, maximum valid expectation value", String(1000));
		writeNote_(os, "input", "output, results", "all");
		writeNote_(os, "input", "output, sort results by", "spectrum");


		writeNote_(os, "input", "protein, taxon", "yeast");
		writeNote_(os, "input", "spectrum, path", input_filename_);
		writeNote_(os, "input", "output, path", output_filename_);

		os << "</bioml>" << endl;
	}

	void XTandemInfile::writeNote_(ostream& os, const String& type, const String& label, const String& value)
	{
		os << "\t<note type=\"" << type << "\" label=\"" << label  << "\">" << value << "</note>" << endl;
	}

	void XTandemInfile::setOutputFilename(const String& filename)
	{
		output_filename_ = filename;
	}

	const String& XTandemInfile::getOutputFilename() const
	{
		return output_filename_;
	}

	void XTandemInfile::setInputFilename(const String& filename)
	{
		input_filename_ = filename;
	}

	const String& XTandemInfile::getInputFilename() const
	{
		return input_filename_;
	}

  					 
} // namespace OpenMS
