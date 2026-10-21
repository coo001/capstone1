#include "openfhe.h"

#include <algorithm>
#include <chrono>
#include <complex>
#include <vector>
#include <iostream>
#include "ckksrns-ship-precom.h"
#include "ckksrns-ship-s2c.h"
#include <type_traits>
#include <utility>
#include <stdexcept>

using namespace std;
using namespace lbcrypto;

static double mae_vec(const vector<complex<double>>& a, const vector<complex<double>>& b, double& maxe) {
    size_t L = min(a.size(), b.size());
    double mae = 0.0; maxe = 0.0;
    for (size_t i = 0; i < L; ++i) {
        double e = std::abs(a[i] - b[i]);
        mae += e; if (e > maxe) maxe = e;
    }
    return L ? (mae / L) : 0.0;
}


void run_minimal() {
    CCParams<CryptoContextCKKSRNS> p;

    const uint32_t ringDim = 1u << 14; 
    const uint32_t slots = ringDim / 2;

    p.SetRingDim(ringDim);
    p.SetBatchSize(slots);

    p.SetScalingTechnique(FLEXIBLEAUTO);
    p.SetScalingModSize(59);
    p.SetFirstModSize(60);

    p.SetKeySwitchTechnique(HYBRID);
    p.SetNumLargeDigits(3);

    p.SetSecurityLevel(HEStd_NotSet);

    vector<uint32_t> levelBudget = { 4, 4 };
    vector<uint32_t> dim = { 0, 0 };
    const SecretKeyDist skDist = UNIFORM_TERNARY;
    const uint32_t pipelineStages = 9;

    usint bootstrapDepth = FHECKKSRNS::GetBootstrapDepth(pipelineStages, levelBudget, skDist);
    usint levelsKeep = 22;
    usint totalDepth = bootstrapDepth + levelsKeep;
    p.SetMultiplicativeDepth(totalDepth);

    auto cc = GenCryptoContext(p);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    cout << "[*] ringDim=" << cc->GetRingDimension()
        << ", slots=" << slots
        << ", depth=" << totalDepth << " (keep=" << levelsKeep
        << ", need=" << bootstrapDepth << ")\n";

    cc->EvalBootstrapSetup(levelBudget, dim, slots, 0, true);

    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    cc->EvalBootstrapKeyGen(kp.secretKey, slots);

    vector<complex<double>> x(slots);
    for (size_t i = 0; i < slots; ++i)
        x[i] = { std::sin(i * 0.01), std::cos(i * 0.01) };

    auto pt = cc->MakeCKKSPackedPlaintext(x);
    auto ct0 = cc->Encrypt(kp.publicKey, pt);

    auto ct_work = ct0;
    int mults_before = 40;
    for (int i = 0; i < mults_before; ++i) {
        ct_work = cc->EvalMult(ct_work, ct_work);
    }

    cout << "[*] Before bootstrap: Level=" << ct_work->GetLevel()
        << ", NoiseScaleDeg=" << ct_work->GetNoiseScaleDeg()
        << ", ScalingFactor=" << ct_work->GetScalingFactor() << "\n";

    cout << "[*] Bootstrapping...\n";
    auto t0 = chrono::high_resolution_clock::now();
    //auto ct_bs = cc->EvalBootstrap(ct_work);
    uint32_t l = ct_work->GetSlots();
    auto& pre = GetBootPrecom(l);
    auto a = cc->EvalCoeffsToSlots(pre.m_U0hatTPreFFT, ct_work);
    auto ct_bs = cc->EvalSlotsToCoeffs(pre.m_U0PreFFT, a);
    auto t1 = chrono::high_resolution_clock::now();
    cout << "[*] Bootstrap time = " << chrono::duration<double>(t1 - t0).count() << " s\n";

    cout << "[*] After bootstrap : Level=" << ct_bs->GetLevel()
        << ", NoiseScaleDeg=" << ct_bs->GetNoiseScaleDeg()
        << ", ScalingFactor=" << ct_bs->GetScalingFactor() << "\n";

    auto towers_before = ct_work->GetElements()[0].GetNumOfElements();
    auto towers_after = ct_bs->GetElements()[0].GetNumOfElements();

    cout << "[*] Towers(before)=" << towers_before
        << ", Towers(after)=" << towers_after << "\n";

    cout << "[*] ScalingFactor(before)=" << ct_work->GetScalingFactor()
        << ", (after)=" << ct_bs->GetScalingFactor() << "\n";

    Plaintext out_before, out_after;
    cc->Decrypt(kp.secretKey, ct_work, &out_before);
    cc->Decrypt(kp.secretKey, ct_bs, &out_after);
    out_before->SetLength(16);
    out_after->SetLength(16);

    auto yb = out_before->GetCKKSPackedValue();
    auto ya = out_after->GetCKKSPackedValue();
    double maxe = 0.0;
    double mae = mae_vec(yb, ya, maxe);

    cout << "[OK] Consistency check: y_before[0]=" << yb[0]
        << ", y_after[0]=" << ya[0] << "\n";
    cout << "[*] MAE(after vs before) = " << mae
        << ", MAXE = " << maxe << "\n";

    auto ct_more = ct_bs;
    bool ok_more = true;
    int mults_after = 20;
    try {
        for (int i = 0; i < mults_after; ++i) {
            ct_more = cc->EvalMult(ct_more, ct_more);
        }
    }
    catch (const exception& e) {
        ok_more = false;
        cerr << "[!] Post-bootstrap extra ops failed: " << e.what() << "\n";
    }

    cout << (ok_more ? "[OK] Post-bootstrap extra ops succeeded"
        : "[X] Post-bootstrap extra ops failed")
        << " | Level=" << ct_more->GetLevel()
        << ", NoiseScaleDeg=" << ct_more->GetNoiseScaleDeg() << "\n";

    if (ok_more) {
        Plaintext out_more;
        cc->Decrypt(kp.secretKey, ct_more, &out_more);
        out_more->SetLength(8);
        auto ym = out_more->GetCKKSPackedValue();
        double maxe2 = 0.0, mae2 = mae_vec(ya, ym, maxe2);
        cout << "[*] MAE(after-extra vs after) = " << mae2
            << ", MAXE = " << maxe2 << "\n";
    }
}

void sophomr() {

    run_minimal();

}