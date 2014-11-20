// From https://github.com/lindenb/cclindenb/blob/master/src/core/lindenb/io/tarball.cpp
/*
 * tarball.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: Pierre Lindenbaum PhD
 *              plindenbaum@yahoo.fr
 *              http://plindenbaum.blogspot.com
 *     License: Apache License, Version 2.0 (as communicated to Adam Tenderholt by Pierre
 *               Lindenbaum PhD
 *
 * Modified by: Adam Tenderholt
 *
 *              + Moved PosixTarHeader to seperate header file
 *              + Renamed LOCALNS::Tar::_checksum() to OpenMS::Internal::headerChecksum()
 *              + Modify headerChecksum() to not automatically store the checksum in the
 *                header, but instead returns it as a uint.
 *              + Removed remaining code
 */

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/FORMAT/INTERNAL/tarball.h>

namespace OpenMS
{
  namespace Internal
  {
    unsigned int headerChecksum(void* header)
    {
      unsigned int sum = 0;
      char *p = (char *) header;
      char *q = p + sizeof(PosixTarHeader);

      while (p < static_cast<PosixTarHeader*>(header)->checksum)
      {
        sum += *p++ & 0xff;
      }

      for (int i = 0; i < 8; ++i)
      {
        sum += ' ';
        ++p;
      }

      while (p < q)
      {
        sum += *p++ & 0xff;
      }

      return sum;
    }
  }
}
