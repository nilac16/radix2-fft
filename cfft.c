#include <math.h>
#include <stdio.h>
#include "cfft.h"


/** @brief Count leading zeros in @p x
 *  @param x
 *      Bit string
 *  @returns The leading zero count in @p x
 */
static int cfft_clz(unsigned x)
{
    int res;

#if __GNUC__
    res = __builtin_clz(x);
#elif MSC_VER
    res = __lzcnt(x);
#endif
    return res;
}


/** @brief Compute the Hamming weight of @p x
 *  @param x
 *      Bit string
 *  @returns The number of set bits in @p x
 */
static int cfft_popcount(unsigned x)
{
    int res;

#if __GNUC__
    res = __builtin_popcount(x);
#elif MSC_VER
    res = __popcnt(x);
#endif
    return res;
}


/** @brief Reverse the bit string composing @p byte
 *  @param byte
 *      Byte to be reversed
 *  @returns The reverse of @p byte
 */
static unsigned char byte_reverse(unsigned char byte)
{
    /* yoink'd from SO
    thanks Matt J */
    static const unsigned char table[] = {
        0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0,
        0x08, 0x88, 0x48, 0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8,
        0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4, 0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4,
        0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C, 0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC,
        0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2, 0x32, 0xB2, 0x72, 0xF2,
        0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A, 0xFA,
        0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6,
        0x0E, 0x8E, 0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE,
        0x01, 0x81, 0x41, 0xC1, 0x21, 0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1,
        0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9, 0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9,
        0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55, 0xD5, 0x35, 0xB5, 0x75, 0xF5,
        0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD, 0x7D, 0xFD,
        0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3,
        0x0B, 0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB,
        0x07, 0x87, 0x47, 0xC7, 0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7,
        0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F, 0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF
    };

    return table[byte];
}


/** @brief Swap bytes
 *  @param b1
 *      Operand address
 *  @param b2
 *      Operand address
 */
static void byteswap(unsigned char *b1, unsigned char *b2)
{
    const unsigned char swap = *b1;

    *b1 = *b2;
    *b2 = swap;
}


/** @brief Swap standard C complex floating-point objects
 *  @param c1
 *      Operand address
 *  @param c2
 *      Operand address
 */
static void complexswap(complex double *c1, complex double *c2)
{
    const complex double swap = *c1;

    *c1 = *c2;
    *c2 = swap;
}


/** @brief Reverse the bit string composing @p x
 *  @param x
 *      Integer to reverse
 *  @returns The reverse of @p x
 */
static unsigned cfft_reverse(unsigned x)
{
    union {
        unsigned char byte[4];
        unsigned      word;
    } u = { .word = x };

    u.byte[0] = byte_reverse(u.byte[0]);
    u.byte[1] = byte_reverse(u.byte[1]);
    u.byte[2] = byte_reverse(u.byte[2]);
    u.byte[3] = byte_reverse(u.byte[3]);
    byteswap(&u.byte[0], &u.byte[3]);
    byteswap(&u.byte[1], &u.byte[2]);
    return u.word;
}


/** @brief Reorder @p data into its bit-reversal permutation
 *  @param len
 *      Length of @p data
 *  @param data
 *      Input data in "natural order"
 */
static void cfft_reorder(unsigned len, complex double data[])
{
    unsigned shift, i, j;

    shift = cfft_clz(len) + 1;
    for (i = 0; i < len; i++) {
        j = cfft_reverse(i) >> shift;
        j = (j < i) ? j : i;
        complexswap(&data[i], &data[j]);
    }
}


/** @brief Returns the principal root of unity for a polynomial of order @p size
 *  @param size
 *      Polynomial order (this can be negative!)
 *  @returns The principal root of z ^ @p size - 1
 *  @note "Principal" in this context means a unit multiple of the complex
 *      argument (i.e. exp(2πi / size))
 */
static complex double root_of_unity(int size)
{
    const double pi = 3.141592653589793, twopi = 2.0 * pi;
    double arg = twopi / (double)size;

    return cos(arg) + sin(arg) * I;
}


/** @brief Reduce an entire row of recursion size @p cur using the Danielson-
 *      Lanczos lemma
 *  @param len
 *      Input data buffer size
 *  @param data
 *      Data buffer
 *  @param princip
 *      Principal root of unity for this recursion level
 *  @param cur
 *      Current size of each Danielson-Lanczos block
 */
static void cfft_reduce(unsigned       len,
                        complex double data[],
                        complex double princip,
                        unsigned       cur)
{
    const unsigned half = cur / 2;
    complex double twiddle, *offs, even, odd;
    unsigned skip, i;

    for (skip = 0; skip < len; skip += cur) {
        twiddle = 1.0;
        offs = data + skip;
        for (i = 0; i < half; i++) {
            even = offs[i];
            odd = twiddle * offs[half + i];
            offs[i] = even + odd;
            offs[half + i] = even - odd;
            twiddle *= princip;
        }
    }
}


int cfft_compute(unsigned len, complex double data[])
{
    complex double princip;
    unsigned cur;

    if (cfft_popcount(len) != 1) {
        return -1;  /* you THOUGHT */
    }
    cfft_reorder(len, data);
    for (cur = 2; cur <= len; cur *= 2) {
        princip = root_of_unity(-(int)cur);
        cfft_reduce(len, data, princip, cur);
    }
    return 0;
}


/** @brief Normalize @p data to its length
 *  @param len
 *      Length of @p data
 *  @param data
 *      Data buffer
 */
static void cfft_normalize(unsigned len, complex double data[])
{
    const double norm = (double)len;
    unsigned i;

    for (i = 0; i < len; i++) {
        data[i] /= norm;
    }
}


int cfft_inverse(unsigned len, complex double data[])
{
    complex double princip;
    unsigned cur;

    if (cfft_popcount(len) != 1) {
        return -1;
    }
    cfft_reorder(len, data);
    for (cur = 2; cur <= len; cur *= 2) {
        princip = root_of_unity(cur);
        cfft_reduce(len, data, princip, cur);
    }
    cfft_normalize(len, data);
    return 0;
}
