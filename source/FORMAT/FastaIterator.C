// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Clemens Groepl, Andreas Bertsch $
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/FastaIterator.h>
#include <iostream>


namespace OpenMS
{

typedef std::pair <String, String> FASTAEntry;

FastaIterator::FastaIterator() :  
  PepIterator(),
  actual_seq_(),
  is_at_end_(false),
  fasta_file_(),
  input_file_()
{
}

FastaIterator::~FastaIterator()
{
	
}

// not implemented (since copying this stuff will invalidate the istream
//FastaIterator::FastaIterator(const FastaIterator & source) : PepIterator(source)
//{
//}

// not implemented (since copying this stuff will invalidate the istream
//FastaIterator& operator=(const FastaIterator &);


PepIterator * FastaIterator::operator++(int)
{ // this operator requires copying, which we cannot support!
	throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
}

FASTAEntry FastaIterator::operator*() 
{
	if (last_header_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	return FASTAEntry (last_header_,actual_seq_);
}

PepIterator & FastaIterator::operator++()
{
	if (last_header_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	actual_seq_ = next_();
	return *this;	
}

void FastaIterator::setFastaFile (const String & f)
{
	std::fstream fs;
	fs.open(f.c_str(), std::fstream::in);
	if (!fs.is_open())
	{
		throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, f);
	}
  fs.close();
	fasta_file_ = f;
}

String FastaIterator::getFastaFile()
{
	return (fasta_file_);
}

std::string FastaIterator::next_()
{
	if (input_file_.eof())
	{
		is_at_end_ = true;
    input_file_.close();
		return ("");
	}
  is_at_end_ = false;
	std::string line;
	std::getline(input_file_, line);
	if (line[0] == '>' || input_file_.eof())
	{
		last_header_ = header_;
		header_ = line;
		return ("");
	}
	return (std::string(line)+next_());
}
	
bool FastaIterator::begin()
{
	if (fasta_file_=="")
	{
		throw Exception::InvalidIterator(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	input_file_.open(fasta_file_.c_str(), std::fstream::in);
	
	if (input_file_)
	{
		std::string line;
		std::getline(input_file_,line);
		header_ = line;
		last_header_ = line;
		actual_seq_ = next_();
		return (true);
	}
	
	return (false);
	
}
	
bool FastaIterator::isAtEnd ()
{
	return (is_at_end_);	
}

} //namespace OpenMS
