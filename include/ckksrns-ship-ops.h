#pragma once
#include "ckksrns-ship-precom.h"
#include <tuple>

namespace lbcrypto {

    struct OmegaPack {
        Plaintext w1_a;
        Plaintext w1_na;
        Plaintext w2_a;
        Plaintext w2_na;
    };

    std::tuple<OmegaPack, Plaintext> MakeOmegaPowersAndT0(CryptoContext<DCRTPoly> cc,
        ConstCiphertext<DCRTPoly>& ctxt,
        const SHIPPrecom& sp);

    Ciphertext<DCRTPoly> ApplyFusedMaskAndColumnRot(CryptoContext<DCRTPoly> cc,
        const OmegaPack& omega,
        const std::map<uint32_t, SHIPMask4>& colKeys,
        uint32_t numGiantSteps,
        uint32_t theta);

    Ciphertext<DCRTPoly> ApplyMuxRot(CryptoContext<DCRTPoly> cc,
        ConstCiphertext<DCRTPoly>& ct,
        const std::vector<Ciphertext<DCRTPoly>>& muxKeys);

    Ciphertext<DCRTPoly> ProductTree(CryptoContext<DCRTPoly> cc,
        std::vector<Ciphertext<DCRTPoly>> terms,
        Plaintext pt_t0);

}