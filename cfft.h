#pragma once

#ifndef CFFT_H
#define CFFT_H

#include <complex.h>


/** @brief Compute the fast Fourier transform of @p data, placing the result
 *      back in @p data
 *  @param len
 *      Length of @p data. This MUST be a power of two. Zero-pad your sample if
 *      needed
 *  @param data
 *      Input data buffer, also the output data buffer
 *  @returns Zero on success, -1 if @p len is not a power of two
 */
int cfft_compute(unsigned len, complex double data[]);


/** @brief Compute the inverse fast Fourier transform of @p data, placing the
 *      result back in @p data
 *  @param len
 *      Length of @p data. This also *must* be a power of two
 *  @param data
 *      Input data buffer, also the output data buffer
 *  @returns Zero on success, -1 if @p len is not a power of two
 */
int cfft_inverse(unsigned len, complex double data[]);


/** @brief Compute the two-dimensional fast Fourier transform of @p data
 *  @param dim
 *      Dimensions of @p data as the pair (width, height). These must both be
 *      powers of two---zero pad rows and columns to ensure this
 *  @param data
 *      The data in a single linear array. The output is written here
 *  @returns Zero on success, -1 if @p dim does not contain powers of two
 */
int cfft2_compute(const unsigned dim[], complex double data[]);


/** @brief Swap the left and right halves of @p data to move the zero-frequency
 *      component to the halfway point
 *  @param len
 *      Length of @p data. This does not need to be a power of two, but it
 *      should be even. If it is odd, then the middle element of @p data will be
 *      left unmoved
 *  @param data
 *      Raw FFT data as returned by cfft_compute
 *  @note This operation is an involution (its own inverse). You do NOT need to
 *      save your raw data before this call if you want to invert it later. Just
 *      apply this function again
 */
void cfft_shift(unsigned len, complex double data[]);


/** @brief Swap diagonal quadrants of @p data to move the zero-frequency
 *      component to the center
 *  @param dim
 *      Dimensions. These also do not need to be a power of two, but should be
 *      even. If a dimension is not even, the middle element will not be moved
 *  @param data
 *      Data buffer
 *  @note Again, this operation is an involution
 */
void cfft2_shift(const unsigned dim[], complex double data[]);


/** @brief Round up @p x to the nearest power of two. This does nothing if @p x
 *      is already a power of two
 *  @param x
 *      Test value
 *  @returns The first power of two not less than @p x
 */
static inline unsigned cfft_next_pow2(unsigned x)
{
    x--;
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return x + 1;
}


#endif /* CFFT_H */
