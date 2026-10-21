#include "ckksrns-ship-precom.h"

namespace lbcrypto {

    static Plaintext MakeMaskPT(CryptoContext<DCRTPoly> cc,
        const std::vector<double>& v, double scale) {
        return cc->MakeCKKSPackedPlaintext(v, 1, 0, scale);
    }

    SHIPPrecom GenerateSHIPPrecom(CryptoContext<DCRTPoly> cc,
        PublicKey<DCRTPoly> pk,
        const lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
        uint32_t slots,
        const std::vector<uint32_t>& suppIndices,
        const std::vector<int>& signs,
        double scale,
        uint32_t theta) {

        if (slots == 0)
            OPENFHE_THROW("GenerateSHIPPrecom: slots must be > 0");
        if (suppIndices.size() != signs.size())
            OPENFHE_THROW("GenerateSHIPPrecom: suppIndices and signs size mismatch");
        
        const uint32_t N = slots * 2;
        SHIPPrecom sp; sp.slots = slots;
        sp.scale = scale;
        sp.indices = suppIndices;
        sp.signs = signs;
        sp.h = suppIndices.size();
        sp.theta = theta;
        sp.brotKeys.colKeys.resize(sp.h);

        const uint32_t numGiantSteps = N / theta;

        for (size_t k = 0; k < sp.h; ++k) {
            int s_j = (sp.signs[k] >= 0) ? +1 : -1;
            uint32_t j = sp.indices[k];

            uint32_t j_0 = j % theta;
            uint32_t j_1 = j / theta;

            std::vector<double> m1(slots), m2(slots), m3(slots), m4(slots);

            if (j < slots) {
                for (uint32_t i = 0; i < slots; ++i) {
                    if (i >= j) {
                        m1[i] = (1.0 + s_j) / 2.0;
                        m2[i] = (1.0 - s_j) / 2.0;
                        m3[i] = 0.0;
                        m4[i] = 0.0;
                    }
                    else {
                        m1[i] = 0.0;
                        m2[i] = 0.0;
                        m3[i] = (1.0 - s_j) / 2.0;
                        m4[i] = (1.0 + s_j) / 2.0;
                    }
                }
            }
            else {
                uint32_t j_prime = j - slots;
                for (uint32_t i = 0; i < slots; ++i) {
                    if (i >= j_prime) {
                        m1[i] = 0.0;
                        m2[i] = 0.0;
                        m3[i] = (1.0 - s_j) / 2.0;
                        m4[i] = (1.0 + s_j) / 2.0;
                    }
                    else {
                        m1[i] = (1.0 - s_j) / 2.0;
                        m2[i] = (1.0 + s_j) / 2.0;
                        m3[i] = 0.0;
                        m4[i] = 0.0;
                    }
                }
            }


            SHIPMask4 fusedMask;
            fusedMask.ct[0] = cc->Encrypt(pk, MakeMaskPT(cc, m1, scale));
            fusedMask.ct[1] = cc->Encrypt(pk, MakeMaskPT(cc, m2, scale));
            fusedMask.ct[2] = cc->Encrypt(pk, MakeMaskPT(cc, m3, scale));
            fusedMask.ct[3] = cc->Encrypt(pk, MakeMaskPT(cc, m4, scale));

            sp.brotKeys.colKeys[k][j_1] = std::move(fusedMask);

        }

        return sp;
    }

}
