#pragma once
#include "openfhe.h"
#include <array>
#include <vector>

namespace lbcrypto {

    struct SHIPMask4 {
        std::array<Ciphertext<DCRTPoly>, 4> ct;
    };

    struct SHIPBlindRotKeys {
        std::vector<std::map<uint32_t, SHIPMask4>> colKeys;
        std::vector<std::vector<Ciphertext<DCRTPoly>>> muxKeys;
    };

    struct SHIPPrecom {
        uint32_t slots = 0;
        size_t   h = 0;
        double   scale = 0.0;
        uint32_t theta = 0;

        std::vector<uint32_t> indices;
        std::vector<int>      signs;

        SHIPBlindRotKeys         brotKeys;
    };

    SHIPPrecom GenerateSHIPPrecom(CryptoContext<DCRTPoly> cc,
        PublicKey<DCRTPoly> pk,
        const lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
        uint32_t slots,
        const std::vector<uint32_t>& suppIndices,
        const std::vector<int>& signs,
        double scale,
        uint32_t theta
    );

}
