#include <bwc.h>

#include "index.h"
#include "archives_compress_bigwhoop.h"

// TODO: find out if BigWhoop also works with int data

namespace vistle {
namespace detail {

template<>
size_t compressBigWhoop<float>(float *src, const Index dim[3], float *compressed, const CompressionSettings &config)
{
    bwc_codec *coder = bwc_alloc_coder(dim[0], dim[1], dim[2], 1, config.p_bigWhoop_nPar, bwc_precision_single);
    bwc_stream *stream = bwc_init_stream(src, compressed, comp);

    bwc_create_compression(coder, stream, const_cast<char *>(config.p_bigWhoop_rate));
    size_t compressed_size = bwc_compress(coder, stream);

    bwc_free_codec(coder);

    return compressed_size;
}

template<>
size_t compressBigWhoop<double>(double *src, const Index dim[3], double *compressed, const CompressionSettings &config)
{
    bwc_codec *coder = bwc_alloc_coder(dim[0], dim[1], dim[2], 1, config.p_bigWhoop_nPar, bwc_precision_double);
    bwc_stream *stream = bwc_init_stream(src, compressed, comp);

    bwc_create_compression(coder, stream, const_cast<char *>(config.p_bigWhoop_rate));
    size_t compressed_size = bwc_compress(coder, stream);

    bwc_free_codec(coder);

    return compressed_size;
}

template<>
bool decompressBigWhoop<float>(float *dest, float *compressed, uint8_t layer)
{
    bwc_codec *decoder = bwc_alloc_decoder();
    bwc_stream *stream = bwc_init_stream(compressed, dest, decomp);

    bwc_create_decompression(decoder, stream, layer);
    bwc_decompress(decoder, stream);
    bwc_free_codec(decoder);

    return true;
}

template<>
bool decompressBigWhoop<double>(double *dest, double *compressed, uint8_t layer)
{
    bwc_codec *decoder = bwc_alloc_decoder();
    bwc_stream *stream = bwc_init_stream(compressed, dest, decomp);

    bwc_create_decompression(decoder, stream, layer);
    bwc_decompress(decoder, stream);
    bwc_free_codec(decoder);

    return true;
}

} // namespace detail
} // namespace vistle
