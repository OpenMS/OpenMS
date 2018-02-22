#pragma warning(push, 0)
/* C:/dev/OpenMS_build/src/openms/Release/cotire/OpenMS_CXX_prefix.hxx */
/* cotire.cmake 1.7.10 generated file */
#ifdef __cplusplus


#include <algorithm>
#include <cassert>
#include <cctype>
#include <clocale>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <numeric>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <ostream>
#include <set>
#include <string>
#include <vector>


#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/Map.h>


#include <QtCore/QByteArray>
#include <QtCore/QCoreApplication>
#include <QtCore/QCryptographicHash>
#include <QtCore/QDate>
#include <QtCore/QDateTime>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QFileSystemWatcher>
#include <QtCore/QList>
#include <QtCore/QObject>
#include <QtCore/QProcess>
#include <QtCore/QRegExp>
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtCore/QTimer>
#include <QtCore/QUrl>
#include <QtGui/QTextDocument>
#include <QtNetwork/QHostInfo>
#include <QtNetwork/QHttpRequestHeader>
#include <QtNetwork/QNetworkReply>
#include <QtNetwork/QNetworkRequest>
#include <WildMagic/Wm5ApprLineFit2.h>
#include <WildMagic/Wm5IntpAkimaNonuniform1.h>
#include <WildMagic/Wm5LinearSystem.h>
#include <WildMagic/Wm5Vector2.h>
// required to avoid: 'error: 'make_array' is not a member of 'boost::serialization'
#include <boost/version.hpp>
#if BOOST_VERSION >= 106400
#  include <boost/serialization/array_wrapper.hpp>
#else
#  include <boost/serialization/array.hpp>
#endif
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/assign.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bind.hpp>
#include <boost/cast.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/function/function_base.hpp>
#include <boost/functional/hash.hpp>
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/tr1.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/unordered_map.hpp>
#include <bzlib.h>
/* 
// lots of warnings in Eigen code...
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Sparse>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
*/
//#include <glpk.h> // this is trouble, since it defines a lot of bad stuff
#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io/guess_stream_format.h>
#include <seqan/seq_io/read_fasta_fastq.h>
#include <seqan/stream.h>
// sqlite not currently required
//#include <sqlite/sqlite3.h>
#include <svm.h>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNode.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/framework/XMLFormatter.hpp>
#include <xercesc/framework/psvi/XSValue.hpp>
#include <xercesc/internal/MemoryManagerImpl.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/Locator.hpp>
#include <xercesc/sax2/Attributes.hpp>
#include <xercesc/sax2/DefaultHandler.hpp>
#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <xercesc/util/BinInputStream.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XMLUniDefs.hpp>
#include <xercesc/validators/common/Grammar.hpp>
#include <zlib.h>

#endif
#pragma warning(pop)
