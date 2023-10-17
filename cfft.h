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


#endif /* CFFT_H */
