
/** Reference implementation of BP-128 bitpacking compression.
 * 
 * These implementations are meant to be easily-readable implementations
 * of the pack and unpack kernels, along with optional data transformations.
 * 
 * Note that manual loop unrolling is recommended to improve performance for a
 * production implementation. 
 * 
 */

#include <cstdint>

// Here we use GCC vector for ease of readability. Each arithmetic
// operation is applied per-element to a length 4 SIMD vector.
typedef uint32_t uint32x4 __attribute__((vector_size(16), aligned(4)));
typedef int32_t int32x4 __attribute__((vector_size(16), aligned(4)));


/////////////////////////////////////////////////////////////////////
// Pack and Unpack kernels
////////////////////////////////////////////////////////////////////


// Pack 128 integers from `in` using `B` bits per integer into `out`
void pack128(const uint32x4 *in, uint32x4 *out, uint8_t B) {
    uint32x4 InReg, OutReg;
    uint32_t mask = (B == 32) ? 0xffffffff : ((1U << B) - 1);
    uint32x4 mask_vec = {mask, mask, mask, mask};
    for (int i = 0; i < 32; i++) {
        unsigned int shift = (i * B) & 31;
        InReg = *in++;

        // This would be the spot to put a transform encode:
        // InReg = transform_encode(InReg);

        InReg = InReg & mask;
        if (shift == 0) OutReg = InReg;
        else OutReg = OutReg | (InReg << shift);

        if (shift >= 32 - B) *out++ = OutReg;
        if (shift > 32 - B) OutReg = InReg >> (32 - shift);
    }
}

// Pack 128 integers from `in` using `B` bits per integer into `out`
void unpack128(const uint32x4 *in, uint32x4 *out, uint8_t B) {
    uint32x4 InReg, OutReg;
    uint32_t mask = (B == 32) ? 0xffffffff : ((1U << B) - 1);
    uint32x4 mask_vec = {mask, mask, mask, mask};
    for (int i = 0; i < 32; i++) {
        unsigned int shift = (i * B) & 31;
        if (shift == 0) OutReg = InReg = *in++;
        else OutReg = InReg >> shift;

        if (shift > 32 - B) {
            InReg = *in++;
            OutReg |= InReg << (32 - shift);
        }
        OutReg &= mask;
        
        // This would be the spot to put a transform decode:
        // OutReg = transform_decode(OutReg);

        *out++ = OutReg;
    }
}



/////////////////////////////////////////////////////////////////////
// Data transformation encode/decode helpers
////////////////////////////////////////////////////////////////////


// init holds the preceding set of 4 values,
// Or simply vec[0] for all 4 values at the start
// of a chunk of 128
uint32x4 delta_encode(uint32x4 vec, uint32x4 init) {
    return (uint32x4) {
        vec[0] - init[3], 
        vec[1] - vec[0], 
        vec[2] - vec[1], 
        vec[3] - vec[2]
    };
}

uint32x4 delta_decode(uint32x4 vec, uint32x4 init) {
    return (uint32x4) {
      init[3] + vec[0],
      init[3] + vec[0] + vec[1],
      init[3] + vec[0] + vec[1] + vec[2],
      init[3] + vec[0] + vec[1] + vec[2] + vec[3]
    };
}

uint32x4 zigzag_encode(int32x4 x) {
    // Need a signed integer to get arithmetic shift right
    return (uint32x4)(x >> 31) ^ (x << 1);
}

uint32x4 zigzag_decode(uint32x4 x) { 
    return (x >> 1) ^ -(x & 1); 
}