#ifndef __SIMPLE_H__
#define __SIMPLE_H__
#include "../lzma/easylzma/compress.h"
#include "../lzma/easylzma/decompress.h"
/* compress a chunk of memory and return a dynamically allocated buffer
 * if successful.  return value is an easylzma error code */
int simpleCompress(elzma_file_format format, const unsigned char * inData, size_t inLen, unsigned char ** outData, size_t * outLen);
/* decompress a chunk of memory and return a dynamically allocated buffer
 * if successful.  return value is an easylzma error code */
int simpleDecompress(elzma_file_format format, const unsigned char * inData, size_t inLen, unsigned char ** outData, size_t * outLen);
#endif