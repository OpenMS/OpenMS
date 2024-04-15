/*
        MSNumpress.cpp
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

#include <algorithm>  // for min() and max() in VS2013
#include <climits>
#include <cmath>
#include <iostream>
#include <OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h>


namespace ms::numpress::MSNumpress {

using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

// This is only valid on systems were ints use more bytes than chars...

const int ONE = 1;
static bool is_big_endian() {
	return *((char*)&(ONE)) == 1;
}
bool IS_BIG_ENDIAN = is_big_endian();



/////////////////////////////////////////////////////////////

static void encodeFixedPoint(
		double fixedPoint, 
		unsigned char *result
) {
	int i;
	unsigned char *fp = (unsigned char*)&fixedPoint;
	for (i=0; i<8; i++) {
		result[i] = fp[IS_BIG_ENDIAN ? (7-i) : i];
	}
}



static double decodeFixedPoint(
		const unsigned char *data
) {
	int i;
	double fixedPoint;
	unsigned char *fp = (unsigned char*)&fixedPoint;
		
	for (i=0; i<8; i++)
	{
		fp[i] = data[IS_BIG_ENDIAN ? (7-i) : i];
	}
	
	return fixedPoint;
}

/////////////////////////////////////////////////////////////

/**
 * Encodes the int x as a number of halfbytes in res. 
 * res_length is incremented by the number of halfbytes, 
 * which will be 1 <= n <= 9
 *
 * see header file for a detailed description of the algorithm.
 */
static void encodeInt(
		const unsigned int x,
		unsigned char* res,
		size_t *res_length	
) {
    // get the bit pattern of a signed int x_inp
	unsigned int m;
	unsigned char i, l; // numbers between 0 and 9

    unsigned int mask = 0xf0000000;
    unsigned int init = x & mask;

	if (init == 0)
	{
		l = 8;
		for (i=0; i<8; i++) 
		{
			m = mask >> (4*i);
			if ((x & m) != 0)
			{
				l = i;
				break;
			}
		}
		res[0] = l;
		for (i=l; i<8; i++)
		{
			res[1+i-l] = static_cast<unsigned char>( x >> (4*(i-l)) );
		}
		*res_length += 1+8-l;

	} 
	else if (init == mask)
	{
		l = 7;
		for (i=0; i<8; i++)
		{
			m = mask >> (4*i);
			if ((x & m) != m)
			{
				l = i;
				break;
			}
		}
		res[0] = l + 8;
		for (i=l; i<8; i++)
		{
			res[1+i-l] = static_cast<unsigned char>( x >> (4*(i-l)) );
		}
		*res_length += 1+8-l;

	} 
	else
	{
		res[0] = 0;
		for (i=0; i<8; i++)
		{
			res[1+i] = static_cast<unsigned char>( x >> (4*i) );
		}
		*res_length += 9;

	}
}



/**
 * Decodes an int from the half bytes in bp. Lossless reverse of encodeInt 
 */
static void decodeInt(
		const unsigned char *data,
		size_t *di,
		size_t max_di,
		size_t *half,
		unsigned int *res
) {
    size_t n, i;
    unsigned int mask, m;
    unsigned char head;
    unsigned char hb;

	if (*half == 0)
	{
		head = data[*di] >> 4;
	} 
	else 
	{
		head = data[*di] & 0xf;
		(*di)++;
	}

	*half = 1-(*half);
	*res = 0;
	
	if (head <= 8)
	{
		n = head;
	} 
	else 
	{ // leading ones, fill n half bytes in res
		n = head - 8;
		mask = 0xf0000000;
		for (i=0; i<n; i++)
		{
			m = mask >> (4*i);
			*res = *res | m;
		}
	}
	
	if (n == 8)
	{
		return;
	}
	
	if (*di + ((8 - n) - (1 - *half)) / 2 >= max_di)
	{
		throw "[MSNumpress::decodeInt] Corrupt input data! ";
	}
	
	for (i=n; i<8; i++)
	{
		if (*half == 0)
		{
			hb = data[*di] >> 4;
		} 
		else
		{
			hb = data[*di] & 0xf;
			(*di)++;
		}
		*res = *res | ( static_cast<unsigned int>(hb) << ((i-n)*4));
		*half = 1 - (*half);
	}
}




/////////////////////////////////////////////////////////////

double optimalLinearFixedPointMass(
		const double *data, 
		size_t dataSize,
        double mass_acc
) {
	if (dataSize < 3)
	{
		return 0; // we just encode the first two points as floats
	}
    // We calculate the maximal fixedPoint we need to achieve a specific mass
    // accuracy. Note that the maximal error we will make by encoding as int is
    // 0.5 due to rounding errors.
    double maxFP = 0.5 / mass_acc;

    // There is a maximal value for the FP given by the int length (32bit)
    // which means we cannot choose a value higher than that. In case we cannot
    // achieve the desired accuracy, return failure (-1).
    double maxFP_overflow = optimalLinearFixedPoint(data, dataSize);
    if (maxFP > maxFP_overflow)
	{
		return -1;
	}
    return maxFP;
}

double optimalLinearFixedPoint(
		const double *data,
		size_t dataSize
) {
	/*
	 * safer impl - apparently not needed though
	 *
	if (dataSize == 0) return 0;
	
	double maxDouble = 0;
	double x;

	for (size_t i=0; i<dataSize; i++) {
		x = data[i];
		maxDouble = max(maxDouble, x);
	}

	return floor(0xFFFFFFFF / maxDouble);
	*/
	if (dataSize == 0)
	{
		return 0;
	}
	if (dataSize == 1)
	{
		return floor(0x7FFFFFFFl / data[0]);
	}
	double maxDouble = max(data[0], data[1]);
	double extrapol;
	double diff;

	for (size_t i=2; i<dataSize; i++) {
		extrapol = data[i-1] + (data[i-1] - data[i-2]);
		diff = data[i] - extrapol;
		maxDouble = max(maxDouble, ceil(abs(diff)+1));
	}

	return floor(0x7FFFFFFFl / maxDouble);
}



size_t encodeLinear(
		const double *data, 
		size_t dataSize, 
		unsigned char *result,
		double fixedPoint
) {
	long long ints[3];
	size_t i, ri;
	unsigned char halfBytes[10];
	size_t halfByteCount;
	size_t hbi;
	long long extrapol;
	int diff;

	//printf("Encoding %d doubles with fixed point %f\n", (int)dataSize, fixedPoint);
	encodeFixedPoint(fixedPoint, result);


	if (dataSize == 0)
	{
		return 8;
	}
	ints[1] = static_cast<long long>(data[0] * fixedPoint + 0.5);
	for (i=0; i<4; i++)
	{
		result[8+i] = (ints[1] >> (i*8)) & 0xff;
	}

	if (dataSize == 1)
	{
		return 12;
	}
	ints[2] = static_cast<long long>(data[1] * fixedPoint + 0.5);
	for (i=0; i<4; i++) {
		result[12+i] = (ints[2] >> (i*8)) & 0xff;
	}

	halfByteCount = 0;
	ri = 16;

	for (i=2; i<dataSize; i++) {
		ints[0] = ints[1];
		ints[1] = ints[2];
		if (MS_NUMPRESS_THROW_ON_OVERFLOW && data[i] * fixedPoint + 0.5 > LLONG_MAX	)
		{
			throw "[MSNumpress::encodeLinear] Next number overflows LLONG_MAX.";
		}

		ints[2] = static_cast<long long>(data[i] * fixedPoint + 0.5);
		extrapol = ints[1] + (ints[1] - ints[0]);

		if (MS_NUMPRESS_THROW_ON_OVERFLOW && 
				(		ints[2] - extrapol > INT_MAX 
					|| 	ints[2] - extrapol < INT_MIN	)) {
			throw "[MSNumpress::encodeLinear] Cannot encode a number that exceeds the bounds of [-INT_MAX, INT_MAX].";
		}

		diff = static_cast<int>(ints[2] - extrapol);
		//printf("%lu %lu %lu,   extrapol: %ld    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
		encodeInt(
				static_cast<unsigned int>(diff), 
				&halfBytes[halfByteCount], 
				&halfByteCount
			);
		/*
		printf("%d (%d):  ", diff, (int)halfByteCount);
		for (size_t j=0; j<halfByteCount; j++) {
			printf("%x ", halfBytes[j] & 0xf);
		}
		printf("\n");
		*/
		
		
		for (hbi=1; hbi < halfByteCount; hbi+=2) {
			result[ri] = static_cast<unsigned char>(
					(halfBytes[hbi-1] << 4) | (halfBytes[hbi] & 0xf)
				);
			//printf("%x \n", result[ri]);
			ri++;
		}
		if (halfByteCount % 2 != 0)
		{
			halfBytes[0] = halfBytes[halfByteCount-1];
			halfByteCount = 1;
		} 
		else
		{
			halfByteCount = 0;
		}
	}
	if (halfByteCount == 1) 
	{
		result[ri] = static_cast<unsigned char>(halfBytes[0] << 4);
		ri++;
	}
	return ri;
}



size_t decodeLinear(
		const unsigned char *data,
		const size_t dataSize,
		double *result
) {
	size_t i;
	size_t ri = 0;
	unsigned int init, buff;
	int diff;
	long long ints[3];
	//double d;
	size_t di;
	size_t half;
	long long extrapol;
	long long y;
	double fixedPoint;
	
	//printf("Decoding %d bytes with fixed point %f\n", (int)dataSize, fixedPoint);

	if (dataSize == 8)
	{
		return 0;
	}
	if (dataSize < 8) 
	{
		throw "[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read fixed point! ";
	}
	fixedPoint = decodeFixedPoint(data);


	if (dataSize < 12) 
	{
		throw "[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read first value! ";
	}
	ints[1] = 0;
	for (i=0; i<4; i++)
	{
		ints[1] = ints[1] | ((0xff & (init = data[8+i])) << (i*8));
	}
	result[0] = ints[1] / fixedPoint;

	if (dataSize == 12)
	{
		return 1;
	}
	if (dataSize < 16) 
	{
		throw "[MSNumpress::decodeLinear] Corrupt input data: not enough bytes to read second value! ";
	}
	ints[2] = 0;
	for (i=0; i<4; i++)
	{
		ints[2] = ints[2] | ((0xff & (init = data[12+i])) << (i*8));
	}
	result[1] = ints[2] / fixedPoint;
		
	half = 0;
	ri = 2;
	di = 16;
	
	//printf("   di     ri      half    int[0]    int[1]    extrapol   diff\n");
	
	while (di < dataSize)
	{
		if (di == (dataSize - 1) && half == 1)
		{
			if ((data[di] & 0xf) == 0x0)
			{
				break;
			}
		}
		//printf("%7d %7d %7d %lu %lu %ld", di, ri, half, ints[0], ints[1], extrapol);
		
		ints[0] = ints[1];
		ints[1] = ints[2];
		decodeInt(data, &di, dataSize, &half, &buff);
		diff = static_cast<int>(buff);

		extrapol = ints[1] + (ints[1] - ints[0]);
		y = extrapol + diff;
		//printf(" %d \n", diff);
		result[ri++] 	= y / fixedPoint;
		ints[2] 		= y;
	}

	return ri;
}



void encodeLinear(
		const std::vector<double> &data, 
		std::vector<unsigned char> &result,
		double fixedPoint
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 5 + 8);
	size_t encodedLength = encodeLinear(&data[0], dataSize, &result[0], fixedPoint);
	result.resize(encodedLength);
}



void decodeLinear(
		const std::vector<unsigned char> &data,
		std::vector<double> &result
) {
	size_t dataSize = data.size();
	result.resize((dataSize - 8) * 2);
	size_t decodedLength = decodeLinear(&data[0], dataSize, &result[0]);
	result.resize(decodedLength);
}

/////////////////////////////////////////////////////////////


size_t encodeSafe(
		const double *data, 
		const size_t dataSize, 
		unsigned char *result
) {
	size_t i, j, ri = 0;
	double latest[3];
	double extrapol, diff;
	const unsigned char *fp; 
	
	//printf("d0 d1 d2 extrapol diff\n");
		
	if (dataSize == 0)
	{
		return ri;
	}
	latest[1] = data[0];
	fp = (unsigned char*)data;
	for (i=0; i<8; i++)
	{
		result[ri++] = fp[IS_BIG_ENDIAN ? (7-i) : i];
	}
	
	if (dataSize == 1)
	{
		return ri;
	}
	latest[2] = data[1];
	fp = (unsigned char*)&(data[1]);
	for (i=0; i<8; i++)
	{
		result[ri++] = fp[IS_BIG_ENDIAN ? (7-i) : i];
	}

	fp = (unsigned char*)&diff;
	for (i=2; i<dataSize; i++) {
		latest[0] = latest[1];
		latest[1] = latest[2];
		latest[2] = data[i];
		extrapol = latest[1] + (latest[1] - latest[0]);
		diff = latest[2] - extrapol;
		//printf("%f %f %f %f %f\n", latest[0], latest[1], latest[2], extrapol, diff);
		for (j=0; j<8; j++)
		{
			result[ri++] = fp[IS_BIG_ENDIAN ? (7-j) : j];
		}
	}
	
	return ri;
}



size_t decodeSafe(
		const unsigned char *data,
		const size_t dataSize,
		double *result
) {
	size_t i, di, ri;
	double extrapol, diff;
	double latest[3];
	unsigned char *fp;
	
	if (dataSize % 8 != 0) 
	{
		throw "[MSNumpress::decodeSafe] Corrupt input data: number of bytes needs to be multiple of 8! ";
	}
	//printf("d0 d1 extrapol diff\td2\n");
	
	try {
		fp = (unsigned char*)&(latest[1]);
		for (i=0; i<8; i++)
		{
			fp[i] = data[IS_BIG_ENDIAN ? (7-i) : i];
		}
		result[0] = latest[1];

		if (dataSize == 8)
		{
			return 1;
		}
		fp = (unsigned char*)&(latest[2]);
		for (i=0; i<8; i++)
		{
			fp[i] = data[8 + (IS_BIG_ENDIAN ? (7-i) : i)];
		}
		result[1] = latest[2];
		
		ri = 2;
		
		fp = (unsigned char*)&diff;
		for (di = 16; di < dataSize; di += 8)
		{
			latest[0] = latest[1];
			latest[1] = latest[2];
			
			for (i=0; i<8; i++)
			{
				fp[i] = data[di + (IS_BIG_ENDIAN ? (7-i) : i)];
			}
			
			extrapol = latest[1] + (latest[1] - latest[0]);
			latest[2] = extrapol + diff;
			
			//printf("%f %f %f %f\t%f \n", latest[0], latest[1], extrapol, diff, latest[2]);
		
			result[ri++] = latest[2];
		}
	} catch (...) {
		throw "[MSNumpress::decodeSafe] Unknown error during decode! ";
	}
	
	return ri;
}

/////////////////////////////////////////////////////////////


size_t encodePic(
		const double *data, 
		size_t dataSize, 
		unsigned char *result
) {
	size_t i, ri;
	unsigned int x;
	unsigned char halfBytes[10];
	size_t halfByteCount;
	size_t hbi;

	//printf("Encoding %d doubles\n", (int)dataSize);

	halfByteCount = 0;
	ri = 0;

	for (i=0; i<dataSize; i++) {
		
		if (MS_NUMPRESS_THROW_ON_OVERFLOW && 
				(data[i] + 0.5 > INT_MAX || data[i] < -0.5)		){
			throw "[MSNumpress::encodePic] Cannot use Pic to encode a number larger than INT_MAX or smaller than 0.";
		}
		x = static_cast<unsigned int>(data[i] + 0.5);
		//printf("%d %d %d,   extrapol: %d    diff: %d \n", ints[0], ints[1], ints[2], extrapol, diff);
		encodeInt(x, &halfBytes[halfByteCount], &halfByteCount);
		
		for (hbi=1; hbi < halfByteCount; hbi+=2) {
			result[ri] = static_cast<unsigned char>(
					(halfBytes[hbi-1] << 4) | (halfBytes[hbi] & 0xf)
				);
			//printf("%x \n", result[ri]);
			ri++;
		}
		if (halfByteCount % 2 != 0) {
			halfBytes[0] = halfBytes[halfByteCount-1];
			halfByteCount = 1;
		} else {
			halfByteCount = 0;
		}
	}
	if (halfByteCount == 1) {
		result[ri] = static_cast<unsigned char>(halfBytes[0] << 4);
		ri++;
	}
	return ri;
}



size_t decodePic(
		const unsigned char *data,
		const size_t dataSize,
		double *result
) {
	size_t ri;
	unsigned int x;
	size_t di;
	size_t half;

	//printf("ri      di      half    dSize   count\n");
	
	half = 0;
	ri = 0;
	di = 0;
	
	while (di < dataSize)
	{
		if (di == (dataSize - 1) && half == 1)
		{
			if ((data[di] & 0xf) == 0x0)
			{
				break;
			}
		}
		
		decodeInt(&data[0], &di, dataSize, &half, &x);
		
		//printf("%7d %7d %7d %7d %7d\n", ri, di, half, dataSize, count);
		
		//printf("count: %d \n", count);
		result[ri++] = static_cast<double>(x);
	}

	return ri;
}



void encodePic(
		const std::vector<double> &data,  
		std::vector<unsigned char> &result
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 5);
	size_t encodedLength = encodePic(&data[0], dataSize, &result[0]);
	result.resize(encodedLength);
}



void decodePic(
		const std::vector<unsigned char> &data,  
		std::vector<double> &result
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 2);
	size_t decodedLength = decodePic(&data[0], dataSize, &result[0]);
	result.resize(decodedLength);
}


/////////////////////////////////////////////////////////////


double optimalSlofFixedPoint(
		const double *data, 
		size_t dataSize
) {
	if (dataSize == 0)
	{
		return 0;
	}
	double maxDouble = 1;
	double x;
	double fp;

	for (size_t i=0; i<dataSize; i++)
	{
		x = log(data[i]+1);
		maxDouble = max(maxDouble, x);
	}

	// here we use 0xFFFE as maximal value as we add 0.5 during encoding (see encodeSlof)
	fp = floor(0xFFFE / maxDouble);

	//cout << "    max val: " << maxDouble << endl;
	//cout << "fixed point: " << fp << endl;

	return fp;
}



size_t encodeSlof(
		const double *data, 
		size_t dataSize, 
		unsigned char *result,
		double fixedPoint
) {
	size_t i, ri;
	double temp;
	unsigned short x;
	encodeFixedPoint(fixedPoint, result);

	ri = 8;
	for (i=0; i<dataSize; i++)
	{
		temp = log(data[i]+1) * fixedPoint;

		if (MS_NUMPRESS_THROW_ON_OVERFLOW && 
				temp > USHRT_MAX		) {
			throw "[MSNumpress::encodeSlof] Cannot encode a number that overflows USHRT_MAX.";
		}

		x = static_cast<unsigned short>(temp + 0.5);
		result[ri++] = x & 0xff;
		result[ri++] = (x >> 8) & 0xff; 
	}
	return ri;
}



size_t decodeSlof(
		const unsigned char *data, 
		const size_t dataSize, 
		double *result
) {
	size_t i, ri;
	unsigned short x;
	double fixedPoint;

	if (dataSize < 8) 
	{
		throw "[MSNumpress::decodeSlof] Corrupt input data: not enough bytes to read fixed point! ";
	}
	ri = 0;
	fixedPoint = decodeFixedPoint(data);

	for (i=8; i<dataSize; i+=2) {
		x = static_cast<unsigned short>(data[i] | (data[i+1] << 8));
		result[ri++] = exp(x / fixedPoint) - 1;
	}
	return ri;
}



void encodeSlof(
		const std::vector<double> &data,  
		std::vector<unsigned char> &result,
		double fixedPoint
) {
	size_t dataSize = data.size();
	result.resize(dataSize * 2 + 8);
	size_t encodedLength = encodeSlof(&data[0], dataSize, &result[0], fixedPoint);
	result.resize(encodedLength);
}



void decodeSlof(
		const std::vector<unsigned char> &data,  
		std::vector<double> &result
) {
	size_t dataSize = data.size();
	result.resize((dataSize - 8) / 2);
	size_t decodedLength = decodeSlof(&data[0], dataSize, &result[0]);
	result.resize(decodedLength);
}

} // namespace ms // namespace numpress // namespace MSNumpress 
