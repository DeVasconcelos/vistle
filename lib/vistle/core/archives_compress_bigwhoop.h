#ifndef VISTLE_CORE_ARCHIVES_COMPRESS_BIGWHOOP_H
#define VISTLE_CORE_ARCHIVES_COMPRESS_BIGWHOOP_H

#include <type_traits>

namespace vistle {
namespace detail {

struct BigWhoopParameters {
    uint8_t nPar;
    char *rate;
};

template<typename T>
size_t compressBigWhoop(T *toCompress, const Index dim[3], T *compressed, const BigWhoopParameters &parameters);

template<>
size_t compressBigWhoop(float *toCompress, const Index dim[3], float *compressed, const BigWhoopParameters &parameters);

template<>
size_t compressBigWhoop(double *toCompress, const Index dim[3], double *compressed,
                        const BigWhoopParameters &parameters);

template<typename T>
void decompressBigWhoop(T *toDecompress, const Index dim[3], T *decompressed, uint8_t layer);

template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
void decompressBigWhoop(T *toDecompress, const Index dim[3], T *decompressed, uint8_t layer);

} // namespace detail
} // namespace vistle

#endif
