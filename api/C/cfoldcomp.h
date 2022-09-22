#ifndef FOLDCOMP_H
#define FOLD_COMP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdbool.h>

int compress(const char* buffer, size_t input_size, char* output, int* output_size);
int decompress(const char* input, size_t input_size, bool use_alt_order, char* output, int* output_size);

#ifdef __cplusplus
}
#endif

#endif