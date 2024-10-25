/*
        MSNumpress.hpp
        johan.teleman@immun.lth.se
        
        This distribution goes under the BSD 3-clause license. If you prefer to use Apache
        version 2.0, that is also available at https://github.com/fickludd/ms-numpress
        Copyright (c) 2013, Johan Teleman
        All rights reserved.

        Redistribution and use in source and binary forms, with or without modification,
        are permitted provided that the following conditions are met:

*         Redistributions of source code must retain the above copyright notice, this list
        of conditions and the following disclaimer.
*        Redistributions in binary form must reproduce the above copyright notice, this
        list of conditions and the following disclaimer in the documentation and/or other
        materials provided with the distribution.
*        Neither the name of the Lund University nor the names of its contributors may be
        used to endorse or promote products derived from this software without specific
        prior written permission.

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
        EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
        OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
        SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
        SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
        OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
        HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
        OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/*
	==================== encodeInt ====================
	Some of the encodings described below use a integer compression referred to simply as 
	
	encodeInt()
 
	The algorithm is similar to other variable length integer encodings,
	such as the SQLite Variable-Length Integers encoding, but it uses half
	bytes in its encoding procedure.

	This encoding works on a 4 byte integer, by truncating initial zeros or ones.
	If the initial (most significant) half byte is 0x0 or 0xf, the number of such 
	halfbytes starting from the most significant is stored in a halfbyte. This initial 
	count is then followed by the rest of the ints halfbytes, in little-endian order. 
	A count halfbyte c of

		0 <= c <= 8 		is interpreted as an initial c 		0x0 halfbytes 
		9 <= c <= 15		is interpreted as an initial (c-8) 	0xf halfbytes

	Example:

	int		c		rest
	0 	=> 	0x8
	-1	=>	0xf		0xf
	 2	=>	0x7		0x2
	23	=>	0x6 	0x7	0x1
	2047	=>	0x5 	0xf 0xf	0xf

	Note that the algorithm returns a char array in which the half bytes are
	stored in the lower 4 bits of each element. Since the first element is a
	count half byte, the maximal length of the encoded data is 9 half bytes
	(1 count half byte + 8 half bytes for a 4-byte integer).

 */

#pragma once

#include <cstddef>
#include <vector>

// defines whether to throw an exception when a number cannot be encoded safely
// with the given parameters
#ifndef MS_NUMPRESS_THROW_ON_OVERFLOW
#define MS_NUMPRESS_THROW_ON_OVERFLOW true
#endif

namespace ms {
namespace numpress {

namespace MSNumpress {
	
	/**
	 * Compute the maximal linear fixed point that prevents integer overflow.
   *
	 * @param data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @param dataSize	number of doubles from *data to encode
   *
	 * @return		the linear fixed point safe to use
	 */
	double optimalLinearFixedPoint(
		const double *data,
		size_t dataSize);
	
	/**
	 * Compute the optimal linear fixed point with a desired m/z accuracy.
     *
     * @note If the desired accuracy cannot be reached without overflowing 64
     * bit integers, then a negative value is returned. You need to check for
     * this and in that case abandon numpress or use optimalLinearFixedPoint
     * which returns the largest safe value.
     *
	 * @param data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @param dataSize	number of doubles from *data to encode
	 * @param mass_acc	desired m/z accuracy in Th
     *
	 * @return		the linear fixed point that satisfies the accuracy requirement (or -1 in case of failure).
	 */
    double optimalLinearFixedPointMass(
            const double *data,
            size_t dataSize,
            double mass_acc);

	/**
	 * Encodes the doubles in data by first using a 
	 *   - lossy conversion to a 4 byte 5 decimal fixed point representation
	 *   - storing the residuals from a linear prediction after first two values
	 *   - encoding by encodeInt (see above) 
	 * 
	 * The resulting binary is maximally 8 + dataSize * 5 bytes, but much less if the 
	 * data is reasonably smooth on the first order.
	 *
	 * This encoding is suitable for typical m/z or retention time binary arrays. 
	 * On a test set, the encoding was empirically show to be accurate to at least 0.002 ppm.
	 *
	 * @param data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @param dataSize	number of doubles from *data to encode
	 * @param result		pointer to where resulting bytes should be stored
	 * @param fixedPoint	the scaling factor used for getting the fixed point repr. 
	 * 				This is stored in the binary and automatically extracted
	 * 				on decoding.
	 * @return		the number of encoded bytes
	 */
	size_t encodeLinear(
		const double *data, 
		const size_t dataSize, 
		unsigned char *result,
		double fixedPoint);
	
	/**
	 * Calls lower level encodeLinear while handling vector sizes appropriately
	 *
	 * @param data		vector of doubles to be encoded
	 * @param result		vector of resulting bytes (will be resized to the number of bytes)
   * @param fixedPoint	the scaling factor used for getting the fixed point repr. This is stored in the binary and automatically extracted on decoding.
	 */
	void encodeLinear(
		const std::vector<double> &data, 
		std::vector<unsigned char> &result,
		double fixedPoint);

	/**
   * Decodes data encoded by encodeLinear. 
	 *
	 * result vector guaranteed to be shorter or equal to (|data| - 8) * 2
	 *
	 * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e.
	 * that the last encoded int does not use the last byte in the data. In addition the last encoded 
	 * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte. 
	 *
	 * @param data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @param dataSize	number of bytes from *data to decode
	 * @param result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded doubles, or -1 if dataSize < 4 or 4 < dataSize < 8
	 */
	size_t decodeLinear(
		const unsigned char *data,
		const size_t dataSize,
		double *result);
	
	/**
	 * Calls lower level decodeLinear while handling vector sizes appropriately
	 *
	 * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e..
	 * that the last encoded int does not use the last byte in the data. In addition the last encoded 
	 * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte. 
	 *
	 * @param data		vector of bytes to be decoded
	 * @param result		vector of resulting double (will be resized to the number of doubles)
	 */
	void decodeLinear(
		const std::vector<unsigned char> &data,
		std::vector<double> &result);
		
/////////////////////////////////////////////////////////////
	
	
	/**
	 * Encodes the doubles in data by storing the residuals from a linear prediction after first two values.
	 * 
	 * The resulting binary is the same size as the input data.
	 *
	 * This encoding is suitable for typical m/z or retention time binary arrays, and is
	 * intended to be used before zlib compression to improve compression.
	 *
	 * @param data		pointer to array of doubles to be encoded (need memorycont. repr.)
	 * @param dataSize	number of doubles from *data to encode
	 * @param result		pointer to were resulting bytes should be stored
	 */
	size_t encodeSafe(
		const double *data, 
		const size_t dataSize, 
		unsigned char *result);
	
	
	/**
	 * Decodes data encoded by encodeSafe. 
	 *
	 * result vector is the same size as the input data.
	 *
	 * Might throw const char* is something goes wrong during decoding.
	 *
	 * @param data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @param dataSize	number of bytes from *data to decode
	 * @param result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded bytes
	 */
	size_t decodeSafe(
		const unsigned char *data,
		const size_t dataSize,
		double *result);
	
/////////////////////////////////////////////////////////////

	/**
	 * Encodes ion counts by simply rounding to the nearest 4 byte integer, 
	 * and compressing each integer with encodeInt. 
	 *
	 * The handleable range is therefore 0 -> 4294967294.
	 * The resulting binary is maximally dataSize * 5 bytes, but much less if the 
	 * data is close to 0 on average.
	 *
	 * @param data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @param dataSize	number of doubles from *data to encode
	 * @param result		pointer to where resulting bytes should be stored
	 * @return		the number of encoded bytes
	 */
	size_t encodePic(
		const double *data, 
		const size_t dataSize, 
		unsigned char *result);
		
	/**
	 * Calls lower level encodePic while handling vector sizes appropriately
	 *
	 * @param data		vector of doubles to be encoded
	 * @param result		vector of resulting bytes (will be resized to the number of bytes)
	 */
	void encodePic(
		const std::vector<double> &data,
		std::vector<unsigned char> &result);

	/**
	 * Decodes data encoded by encodePic
	 *
	 * result vector guaranteed to be shorter of equal to |data| * 2
	 *
	 * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e.
	 * that the last encoded int does not use the last byte in the data. In addition the last encoded 
	 * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte. 
	 *
	 * @param data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @param dataSize	number of bytes from *data to decode
	 * @param result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded doubles
	 */
	size_t decodePic(
		const unsigned char *data,
		const size_t dataSize,
		double *result);
	
	/**
	 * Calls lower level decodePic while handling vector sizes appropriately
	 *
	 * Note that this method may throw a const char* if it deems the input data to be corrupt, i.e.
	 * that the last encoded int does not use the last byte in the data. In addition the last encoded 
	 * int need to use either the last halfbyte, or the second last followed by a 0x0 halfbyte. 
	 *
	 * @param data		vector of bytes to be decoded
	 * @param result		vector of resulting double (will be resized to the number of doubles)
	 */
	void decodePic(
		const std::vector<unsigned char> &data,
		std::vector<double> &result);

/////////////////////////////////////////////////////////////


	double optimalSlofFixedPoint(
		const double *data, 
		size_t dataSize);

	/**
	 * Encodes ion counts by taking the natural logarithm, and storing a
	 * fixed point representation of this. This is calculated as
	 * 
	 * unsigned short fp = log(d + 1) * fixedPoint + 0.5
	 *
	 * the result vector is exactly |data| * 2 + 8 bytes long
	 *
	 * @param data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @param dataSize	number of doubles from *data to encode
	 * @param result		pointer to were resulting bytes should be stored
   * @param fixedPoint	the scaling factor used for getting the fixed point repr. This is stored in the binary and automatically extracted on decoding.
	 * @return		the number of encoded bytes
	 */
	size_t encodeSlof(
		const double *data, 
		const size_t dataSize, 
		unsigned char *result,
		double fixedPoint);
		
	/**
	 * Calls lower level encodeSlof while handling vector sizes appropriately
	 *
	 * @param data		vector of doubles to be encoded
	 * @param result		vector of resulting bytes (will be resized to the number of bytes)
   * @param fixedPoint	the scaling factor used for getting the fixed point repr. This is stored in the binary and automatically extracted on decoding.
	 */
	void encodeSlof(
		const std::vector<double> &data,
		std::vector<unsigned char> &result,
		double fixedPoint);

	/**
	 * Decodes data encoded by encodeSlof
	 *
	 * The return will include exactly (|data| - 8) / 2 doubles.
	 *
	 * Note that this method may throw a const char* if it deems the input data to be corrupt.
	 *
	 * @param data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @param dataSize	number of bytes from *data to decode
	 * @param result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded doubles
	 */
	size_t decodeSlof(
		const unsigned char *data, 
		const size_t dataSize, 
		double *result);
	
	/**
	 * Calls lower level decodeSlof while handling vector sizes appropriately
	 *
	 * Note that this method may throw a const char* if it deems the input data to be corrupt.
	 *
	 * @param data		vector of bytes to be decoded
	 * @param result		vector of resulting double (will be resized to the number of doubles)
	 */
	void decodeSlof(
		const std::vector<unsigned char> &data,
		std::vector<double> &result);

} // namespace MSNumpress
} // namespace msdata
} // namespace pwiz

