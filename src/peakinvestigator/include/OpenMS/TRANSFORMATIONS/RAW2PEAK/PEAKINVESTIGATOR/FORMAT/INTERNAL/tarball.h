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
 *              + Renamed to tarball.h
 *              + Place PosixTarHeader in OpenMS::Internal
 *              + Removed remaining code
 */

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_FORMAT_INTERNAL_TARBALL_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKINVESTIGATOR_FORMAT_INTERNAL_TARBALL_H

namespace OpenMS
{
  namespace Internal
  {
    struct PosixTarHeader
    {
        char name[100];
        char mode[8];
        char uid[8];
        char gid[8];
        char size[12];
        char mtime[12];
        char checksum[8];
        char typeflag[1];
        char linkname[100];
        char magic[6];
        char version[2];
        char uname[32];
        char gname[32];
        char devmajor[8];
        char devminor[8];
        char prefix[155];
        char pad[12];
    };

    unsigned int headerChecksum(void* header);

  }
}

#endif
