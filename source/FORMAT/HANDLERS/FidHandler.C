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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FidHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <climits>

using namespace std;

#ifdef OPENMS_BIG_ENDIAN
  template<typename T> T ByteReverse(const T in)
  {
    T out;
    const char* pin = (const char*) &in;
    char* pout= (char*) (&out+1) - 1 ;

    int i;
    for(i= sizeof(T) ; i>0 ; --i)
    {
      *pout-- = *pin++ ;
    }
    return out ;
  }
#endif

namespace OpenMS
{
	namespace Internal
	{
      FidHandler::FidHandler(const String& filename)
        : ifstream(filename.c_str(), ios_base::binary | ios_base::in)
      {
        index_ = 0;
        seekg(0, ios::beg);
			}

      FidHandler::~FidHandler()
      {
      }
      
      Size FidHandler::getIndex()
      {
        return index_;
      }
      
      Size FidHandler::getIntensity()
      {
        // intensity is coded in 32 bits little-endian integer format
        Int32 result = 0;
        read( (char*) &result, 4);
        #ifdef OPENMS_BIG_ENDIAN
          result = ByteReverse<Int32>(result);
        #endif
        index_++;        
        return (result > 0) ? result : 0;
      }
	} // namespace Internal
} // namespace OpenMS
