#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/*
 * Decode the CBF x-CBF_BYTE_OFFSET stream used by the frozen 9Z6F frames.
 *
 * Return values:
 *   0  success
 *  -1  compressed input ended early
 *  -2  decoded cumulative value does not fit signed 32 bits
 *  -3  trailing bytes remain after the declared element count
 */
int cbf_byte_offset_decode(const uint8_t *input,
                           size_t input_size,
                           int32_t *output,
                           size_t output_count,
                           size_t *bytes_used) {
    size_t in_pos = 0;
    int64_t previous = 0;

    for (size_t out_pos = 0; out_pos < output_count; ++out_pos) {
        if (in_pos + 1 > input_size) {
            return -1;
        }

        int8_t delta8;
        memcpy(&delta8, input + in_pos, sizeof(delta8));
        in_pos += sizeof(delta8);
        int64_t delta = delta8;

        if (delta8 == INT8_MIN) {
            if (in_pos + 2 > input_size) {
                return -1;
            }
            int16_t delta16;
            memcpy(&delta16, input + in_pos, sizeof(delta16));
            in_pos += sizeof(delta16);
            delta = delta16;

            if (delta16 == INT16_MIN) {
                if (in_pos + 4 > input_size) {
                    return -1;
                }
                int32_t delta32;
                memcpy(&delta32, input + in_pos, sizeof(delta32));
                in_pos += sizeof(delta32);
                delta = delta32;

                if (delta32 == INT32_MIN) {
                    if (in_pos + 8 > input_size) {
                        return -1;
                    }
                    int64_t delta64;
                    memcpy(&delta64, input + in_pos, sizeof(delta64));
                    in_pos += sizeof(delta64);
                    delta = delta64;
                }
            }
        }

        previous += delta;
        if (previous < INT32_MIN || previous > INT32_MAX) {
            return -2;
        }
        output[out_pos] = (int32_t)previous;
    }

    if (bytes_used != NULL) {
        *bytes_used = in_pos;
    }
    return in_pos == input_size ? 0 : -3;
}
