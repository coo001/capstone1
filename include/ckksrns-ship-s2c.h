#pragma once
#include "ckksrns-ship-precom.h"
#include "ckksrns-ship-ops.h"

namespace lbcrypto {
    Ciphertext<DCRTPoly> EvalCoeffsToSlotsSHIP(
        ConstCiphertext<DCRTPoly>& ctxt,
        const SHIPPrecom& sp);
}