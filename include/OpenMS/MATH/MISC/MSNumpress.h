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
 
	This encoding works on a 4 byte integer, by truncating initial zeros or ones.
	If the initial (most significant) half byte is 0x0 or 0xf, the number of such 
	halfbytes starting from the most significant is stored in a halfbyte. This initial 
	count is then followed by the rest of the ints halfbytes, in little-endian order. 
	A count halfbyte c of
	
		0 <= c <= 8 		is interpreted as an initial c 		0x0 halfbytes 
		9 <= c <= 15		is interpreted as an initial (c-8) 	0xf halfbytes

	Ex:
	int		c		rest
	0 	=> 	0x8
	-1	=>	0xf		0xf
	23	=>	0x6 	0x7	0x1
 */

#ifndef OPENMS_MATH_MISC_MSNUMPRESS_H
#define OPENMS_MATH_MISC_MSNUMPRESS_H

#include <cstddef>
#include <vector>

namespace ms {
namespace numpress {

namespace MSNumpress {
	
	double optimalLinearFixedPoint(
		const double *data, 
		size_t dataSize);
	
	/**
	 * Encodes the doubles in data by first using a 
	 *   - lossy conversion to a 4 byte 5 decimal fixed point representation
	 *   - storing the residuals from a linear prediction after first to values
	 *   - encoding by encodeInt (see above) 
	 * 
	 * The resulting binary is maximally dataSize * 5 bytes, but much less if the 
	 * data is reasonably smooth on the first order.
	 *
	 * This encoding is suitable for typical m/z or retention time binary arrays. 
	 * For masses above 100 m/z the encoding is accurate to at least 0.1 ppm.
	 *
	 * @data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @dataSize	number of doubles from *data to encode
	 * @result		pointer to where resulting bytes should be stored
	 * @fixedPoint	the scaling factor used for getting the fixed point repr. 
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
	 * @data		vector of doubles to be encoded
	 * @result		vector of resulting bytes (will be resized to the number of bytes)
	 */
	void encodeLinear(
		const std::vector<double> &data, 
		std::vector<unsigned char> &result,
		double fixedPoint);

	/**
	 * Decodes data encoded by encodeLinear. Note that the compression 
	 * discard any information < 1e-5, so data is only guaranteed 
	 * to be within +- 5e-6 of the original value.
	 *
	 * Further, values > ~42000 will also be truncated because of the
	 * fixed point representation, so this scheme is strongly discouraged 
	 * if values above might be above this size.
	 *
	 * result vector guaranteed to be shorter than twice the data length (in nbr of values)
	 *
	 * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @dataSize	number of bytes from *data to decode
	 * @result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded doubles, or -1 if dataSize < 4 or 4 < dataSize < 8
	 */
	size_t decodeLinear(
		const unsigned char *data,
		const size_t dataSize,
		double *result);
	
	/**
	 * Calls lower level decodeLinear while handling vector sizes appropriately
	 *
	 * @data		vector of bytes to be decoded
	 * @result		vector of resulting double (will be resized to the number of doubles)
	 */
	void decodeLinear(
		const std::vector<unsigned char> &data,
		std::vector<double> &result);

/////////////////////////////////////////////////////////////

	/**
	 * Encodes ion counts by simply rounding to the nearest 4 byte integer, 
	 * and compressing each integer with encodeInt. 
	 *
	 * The handleable range is therefore 0 -> 4294967294.
	 * The resulting binary is maximally dataSize * 5 bytes, but much less if the 
	 * data is close to 0 on average.
	 *
	 * @data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @dataSize	number of doubles from *data to encode
	 * @result		pointer to were resulting bytes should be stored
	 * @return		the number of encoded bytes
	 */
	size_t encodePic(
		const double *data, 
		const size_t dataSize, 
		unsigned char *result);
		
	/**
	 * Calls lower level encodePic while handling vector sizes appropriately
	 *
	 * @data		vector of doubles to be encoded
	 * @result		vector of resulting bytes (will be resized to the number of bytes)
	 */
	void encodePic(
		const std::vector<double> &data,
		std::vector<unsigned char> &result);

	/**
	 * Decodes data encoded by encodePic
	 *
	 * result vector guaranteed to be shorter than twice the data length (in nbr of values)
	 *
	 * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @dataSize	number of bytes from *data to decode
	 * @result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded doubles
	 */
	void decodePic(
		const std::vector<unsigned char> &data,
		std::vector<double> &result);
	
	/**
	 * Calls lower level decodePic while handling vector sizes appropriately
	 *
	 * @data		vector of bytes to be decoded
	 * @result		vector of resulting double (will be resized to the number of doubles)
	 */
	size_t decodePic(
		const unsigned char *data,
		const size_t dataSize,
		double *result);

/////////////////////////////////////////////////////////////


	double optimalSlofFixedPoint(
		const double *data, 
		size_t dataSize);

	/**
	 * Encodes ion counts by taking the natural logarithm, and storing a
	 * fixed point representation of this. This is calculated as
	 * 
	 * unsigned short fp = log(d + 1) * 3000.0 + 0.5
	 *
	 * Note that this fixed point will mean any d < 0.00016667 will be 
	 * stored as a zero and mapped back to a zero. 
	 *
	 * result vector is exactly twice the data length (in nbr of values)
	 *
	 * @data		pointer to array of double to be encoded (need memorycont. repr.)
	 * @dataSize	number of doubles from *data to encode
	 * @result		pointer to were resulting bytes should be stored
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
	 * @data		vector of doubles to be encoded
	 * @result		vector of resulting bytes (will be resized to the number of bytes)
	 */
	void encodeSlof(
		const std::vector<double> &data,
		std::vector<unsigned char> &result,
		double fixedPoint);

	/**
	 * Decodes data encoded by encodeSlof
	 *
	 * @data		pointer to array of bytes to be decoded (need memorycont. repr.)
	 * @dataSize	number of bytes from *data to decode
	 * @result		pointer to were resulting doubles should be stored
	 * @return		the number of decoded doubles
	 */
	size_t decodeSlof(
		const unsigned char *data, 
		const size_t dataSize, 
		double *result);
	
	/**
	 * Calls lower level decodeSlof while handling vector sizes appropriately
	 *
	 * @data		vector of bytes to be decoded
	 * @result		vector of resulting double (will be resized to the number of doubles)
	 */
	void decodeSlof(
		const std::vector<unsigned char> &data,
		std::vector<double> &result);

} // namespace MSNumpress
} // namespace msdata
} // namespace pwiz

#endif // _MSNUMPRESS_HPP_
