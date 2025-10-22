#include "ckksrns-ship-s2c.h"
#include "openfhe.h"
#include <tuple>

namespace lbcrypto {

    Ciphertext<DCRTPoly> EvalCoeffsToSlotsSHIP(
        ConstCiphertext<DCRTPoly>& ctxt,
        const SHIPPrecom& sp) {

        if (!ctxt) {
            OPENHE_THROW("EvalCoeffsToSlotsSHIP: input ciphertext is null");
        }
        if (sp.slots == 0 || sp.h == 0 || sp.indices.size() != sp.h || sp.signs.size() != sp.h) {
            OPENHE_THROW("EvalCoeffsToSlotsSHIP: invalid SHIPPrecom");
        }

        auto cc = ctxt->GetCryptoContext();

        auto [omega, pt_t0] = MakeOmegaPowersAndT0(cc, ctxt, sp);

        std::vector<Ciphertext<DCRTPoly>> pieces;
        pieces.reserve(sp.h);

        const uint32_t numGiantSteps = (sp.slots * 2) / sp.theta;

        for (size_t j = 0; j < sp.h; ++j) {

            auto masked_col_rot = ApplyFusedMaskAndColumnRot(cc,
                omega,
                sp.brotKeys.colKeys[j],
                numGiantSteps,
                sp.theta);

            auto rotated = ApplyMuxRot(cc,
                masked_col_rot,
                sp.brotKeys.muxKeys[j]);

            if (rotated->GetNoiseScaleDeg() >= 3) {
                cc->ModReduceInPlace(rotated);
            }

            pieces.push_back(std::move(rotated));
        }

        if (pieces.empty())
            OPENHE_THROW("EvalCoeffsToSlotsSHIP: no pieces produced for product tree");

        Ciphertext<DCRTPoly> acc = ProductTree(cc, std::move(pieces), std::move(pt_t0));

        auto conj = cc->EvalConj(acc);
        cc->EvalAddInPlace(acc, conj);

        return acc;
    }

}