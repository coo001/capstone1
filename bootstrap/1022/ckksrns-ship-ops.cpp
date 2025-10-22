#include "ckksrns-ship-ops.h"
#include <tuple>

namespace lbcrypto {

    static std::vector<double> ExtractCoeffsAndDecode(ConstCiphertext<DCRTPoly>& ctxt,
        bool getA) {
        // ¡¶¿€¡ﬂ

        OPENHE_THROW("ExtractCoeffsAndDecode: Not implemented.");
        return std::vector<double>();
    }

    std::tuple<OmegaPack, Plaintext> MakeOmegaPowersAndT0(CryptoContext<DCRTPoly> cc,
        ConstCiphertext<DCRTPoly>& ctxt,
        const SHIPPrecom& sp) {

        std::vector<double> a_coeffs = ExtractCoeffsAndDecode(ctxt, true);
        std::vector<double> b_coeffs = ExtractCoeffsAndDecode(ctxt, false);

        const double PI = 3.141592653589793;
        const double q0 = 0.0;

        std::vector<double> v_w1_a(sp.slots), v_w1_na(sp.slots);
        std::vector<double> v_w2_a(sp.slots), v_w2_na(sp.slots);
        std::vector<double> v_t0(sp.slots);

        for (uint32_t i = 0; i < sp.slots; ++i) {
            double a_i = a_coeffs[i];
            double a_i_upper = a_coeffs[i + sp.slots];

            v_w1_a[i] = std::cos(2 * PI * a_i / q0);
            v_w1_na[i] = std::cos(2 * PI * (-a_i) / q0);
            v_w2_a[i] = std::cos(2 * PI * a_i_upper / q0);
            v_w2_na[i] = std::cos(2 * PI * (-a_i_upper) / q0);

            v_t0[i] = (gamma / (4 * PI)) * std::cos(2 * PI * b_coeffs[i] / q0); //
        }

        OmegaPack omega = {
            cc->MakeCKKSPackedPlaintext(v_w1_a, 1, 0, sp.scale),
            cc->MakeCKKSPackedPlaintext(v_w1_na, 1, 0, sp.scale),
            cc->MakeCKKSPackedPlaintext(v_w2_a, 1, 0, sp.scale),
            cc->MakeCKKSPackedPlaintext(v_w2_na, 1, 0, sp.scale)
        };
        Plaintext pt_t0 = cc->MakeCKKSPackedPlaintext(v_t0, 1, 0, sp.scale);

        return std::make_tuple(std::move(omega), std::move(pt_t0));
    }

    Ciphertext<DCRTPoly> ApplyFusedMaskAndColumnRot(CryptoContext<DCRTPoly> cc,
        const OmegaPack& omega,
        const std::map<uint32_t, SHIPMask4>& colKeys,
        uint32_t numGiantSteps,
        uint32_t theta) {

        OPENHE_THROW("ApplyFusedMaskAndColumnRot: Not implemented.");
        return Ciphertext<DCRTPoly>();
    }

    Ciphertext<DCRTPoly> ApplyMuxRot(CryptoContext<DCRTPoly> cc,
        ConstCiphertext<DCRTPoly>& ct,
        const std::vector<Ciphertext<DCRTPoly>>& muxKeys) {

        OPENHE_THROW("ApplyMuxRot: Not implemented.");
        return ct->Clone();
    }


    Ciphertext<DCRTPoly> ProductTree(CryptoContext<DCRTPoly> cc,
        std::vector<Ciphertext<DCRTPoly>> terms,
        Plaintext pt_t0) {

        if (terms.empty())
            OPENHE_THROW("ProductTree: empty terms");

        auto first_term = cc->EvalMult(pt_t0, terms[0]);

        std::vector<Ciphertext<DCRTPoly>> cur;
        cur.reserve((terms.size() + 1) / 2);
        cur.push_back(std::move(first_term));

        cur.insert(cur.end(),
            std::make_move_iterator(terms.begin() + 1),
            std::make_move_iterator(terms.end()));

        while (cur.size() > 1) {
            std::vector<Ciphertext<DCRTPoly>> next;
            next.reserve((cur.size() + 1) / 2);
            for (size_t i = 0; i + 1 < cur.size(); i += 2) {
                auto m = cc->EvalMult(cur[i], cur[i + 1]);
                next.push_back(std::move(m));
            }
            if (cur.size() % 2 == 1)
                next.push_back(cur.back());
            cur.swap(next);
        }
        return cur.front();
    }

}