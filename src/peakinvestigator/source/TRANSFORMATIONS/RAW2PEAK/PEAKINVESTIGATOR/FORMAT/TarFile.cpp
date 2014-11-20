// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer:$
// $Author: Adam Tenderholt $
// --------------------------------------------------------------------------
//

#include <zlib.h>

#include <QtCore/QBuffer>
#include <QtCore/QString>
#include <QtCore/QTextStream>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/FORMAT/TarFile.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PEAKINVESTIGATOR/FORMAT/INTERNAL/tarball.h>

namespace OpenMS
{

  typedef Internal::PosixTarHeader TarHeader;

  TarFile::TarFile(): ProgressLogger()
  {
	  setLogType(ProgressLogger::CMD);
  }

  TarFile::~TarFile()
  {
  }

  void TarFile::load(const String &filename, MSExperiment<Peak1D> &experiment)
  {
    char temp[8192];
    int numBytes;

    gzFile file = gzopen(filename.c_str(), "rb");

    if (file == NULL)
    {
      LOG_FATAL_ERROR << "Problem opening " << filename << std::endl;
      return;
    }

    int i = 0;
    do
    {
      // read header of entry i
      TarHeader header;
      numBytes = gzread(file, &header, sizeof(TarHeader));

      if (numBytes == 0) // end of file
      {
        break;
      }

      else if (numBytes != sizeof(TarHeader))
      {
        LOG_ERROR << "Problem parsing header for entry " << i << std::endl;
        return;
      }

      bool ok;
      qulonglong fileSize = QString(header.size).left(11).toULongLong(&ok, 8);
      if (!ok) // ugly hack to handle end of tarfile, as well as bad-formed header
      {
        continue;
      }
      int remainder = sizeof(TarHeader) - (fileSize % sizeof(TarHeader)); // entries are mulitply NULL terminated

      QBuffer* buffer = new QBuffer;
      buffer->open(QIODevice::WriteOnly);

      while (fileSize >= sizeof(temp))
      {
        numBytes = gzread(file, temp, sizeof(temp));
        buffer->write(temp, numBytes);
        fileSize -= numBytes;
      }
      // get leftover
      gzread(file, temp, fileSize);
      buffer->write(temp, fileSize);

      // read rest of entry (NULLs)
      gzread(file, temp, remainder);

      buffer->close();

      // now load the buffer into a Spectrum, and place it in the MSExperiment
      int num;
      int count = sscanf(header.name, "scan_%6d", &num);
      if (count == 1 && (Size) num < experiment.size())
      {
        LOG_DEBUG << header.name << " loading scan #" << num << ".\n";
        loadDataFromBuffer_(buffer, experiment[num]);
        i++;
      }
      else if (count == 1)
      {
        LOG_WARN << "Loaded scan #" << num << ", but the experiment does not have that many scans.\n";
      }
      else
      {
        LOG_WARN << "File entry " << header.name << " is not of expected format: e.g. 'scan_000102.txt'.\n";
      }

      delete buffer;


    } while (true);

    LOG_DEBUG << "Processed " << i << " files." << std::endl;

    gzclose(file);

  }

  void TarFile::store(const String& filename, const MSExperiment<Peak1D>& experiment)
  {
    int numBytes;
    char temp[8192];

    gzFile file = gzopen(filename.c_str(), "wb");
    TarHeader header;

    startProgress(0, experiment.size(), "Bundling scans for upload");
    for (Size i = 0; i < experiment.size(); ++i)
    {
      QBuffer* buffer = saveDataToBuffer_(experiment[i]);

      // initialize archive header
      std::memset(&header, 0, sizeof(TarHeader)); // fill with NULL chars

      // set archive filename
      std::sprintf(header.name, "scan_%06i.txt", (int) i);

      // set tar format
      std::sprintf(header.magic, "ustar");
      std::memcpy(header.version, "  ", sizeof(header.version));

      // set modification time, mode, and filetype
      std::sprintf(header.mtime, "%011lo", time(NULL));
      std::sprintf(header.mode, "%07o", 0644);
      header.typeflag[0] = 0;

      // set size of file
      qint64 size = buffer->size();
      std::sprintf(header.size, "%011llo", (long long unsigned int) size);

      // set header checksum
      std::sprintf(header.checksum, "%06o", Internal::headerChecksum(&header));

      // now write archive file header (required before writing data)
      numBytes = gzwrite(file, &header, sizeof(TarHeader));
      if(numBytes != sizeof(TarHeader))
      {
        LOG_ERROR << "Not all of the header was written for scan " << i << "!\n";
        delete buffer;
        continue;
      }

      // copy buffer to file in archive
      buffer->open(QBuffer::ReadOnly);
      numBytes = buffer->read(temp, sizeof(temp));
      while (numBytes > 0)
      {
        gzwrite(file, temp, numBytes);
        numBytes = buffer->read(temp, sizeof(temp));
      }
      buffer->close();

      // fill remaining 512-byte block with NULL
      while (size % sizeof(TarHeader) != 0)
      {
        gzputc(file, 0);
        ++size;
      }

      delete buffer;

      // occasionally update progress
      if (i % 10 == 0)
      {
        setProgress(i);
      }
    }

    // now close out tar format by writing with two NULL header entries
    std::memset(&header, 0, sizeof(TarHeader)); // fill with NULL chars
    gzwrite(file, &header, sizeof(TarHeader));
    gzwrite(file, &header, sizeof(TarHeader));

    gzclose(file);

    endProgress();

    return;
  }

  void TarFile::loadDataFromBuffer_(QBuffer *buffer, MSSpectrum<Peak1D> &peaklist)
  {
    buffer->open(QIODevice::ReadOnly);

    double mz, counts;
    QString line;
    bool ok1, ok2;

    QTextStream in(buffer);
    while(!in.atEnd())
    {
      line = in.readLine();

      // skip comments
      if (line[0] == '#')
      {
        continue;
      }

      mz = line.section('\t', 0, 0).toDouble(&ok1);
      counts = line.section('\t', 1, 1).toDouble(&ok2);

      if (!ok1 || !ok2)
      {
        LOG_WARN << "Problem decoding line to buffer: '" << line.toAscii().data() << "'.\n";
      }

      Peak1D peak;
      peak.setMZ(mz);
      peak.setIntensity(counts);
      peaklist.push_back(peak);
    }

    peaklist.updateRanges();
    buffer->close();

  }

  QBuffer* TarFile::saveDataToBuffer_(const MSSpectrum<Peak1D> &spectrum)
  {
    QBuffer* buffer = new QBuffer;
    buffer->open(QIODevice::WriteOnly | QIODevice::Text);

    QTextStream out(buffer);

    for (Size i = 0; i < spectrum.size(); i++)
    {
      out << QString("%1").arg(spectrum[i].getMZ(), 0, 'f', 6) << '\t';
      out << QString("%1").arg(spectrum[i].getIntensity(), 0, 'f', 6) << endl;
    }

    buffer->close();

    return buffer;
  }

}
