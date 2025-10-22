#include "openfhe.h"

//#include "signal.h"
//#include "global.h"
//#include "param.h"
//#include "setup.h"
//#include "detect.h"
//#include "compress.h"
//#include "decode.h"
//#include "enco.h"
#include <algorithm>
#include <chrono>
#include <complex>
#include <vector>
#include <iostream>
#include "ckksrns-ship-precom.h"
#include "ckksrns-ship-s2c.h"

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

/*
void FHECKKSRNS::EvalBootstrapSetupSHIP(const BootstrapShipParams& p) {
    auto cc = GetCryptoContext().template GetContext<DCRTPoly>();
    auto kp = cc->KeyGenIfNeeded();    // 키 없는 경우 생성(또는 외부에서 보장)
    auto slots = cc->GetRingDimension() / 2;

    double sf = p.targetSf > 0 ? p.targetSf : GetDefaultBootstrapTargetScale();

    auto sp = std::make_shared<SHIPPrecom>(GenerateSHIPPrecom(
        cc, kp.publicKey, kp, slots,
        p.suppIdx, p.signs, sf, p.theta,
            p.blind   // 필요 시 생성 경로 분기
            ));
            impl_->shipPrecom = std::move(sp);
            impl_->shipParams = p;
}

*/

Ciphertext<DCRTPoly> EvalBootstrapStd(CryptoContext<DCRTPoly> cc,
    const ConstCiphertext<DCRTPoly>& ct) {
    auto ct_l = ct;
    return cc->EvalBootstrap(ct_l);
}

/*
Ciphertext<DCRTPoly> EvalBootstrapSHIP(CryptoContext<DCRTPoly> cc,
    ConstCiphertext<DCRTPoly>& ct,
    const SHIPPrecom& sp) {
    auto algo = cc->GetScheme();
    
    //auto compositeDegree =
    //    std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ct->GetCryptoParameters())->GetCompositeDegree();

    // (1) 필요하면 이곳에서 raise/modreduce (공개 API만 사용)
    auto raised = ct->Clone();
    //algo->ModReduceInternalInPlace(raised, compositeDegree);

    
    
    // (2) Coeffs -> Slots : SHIP
    auto ctSlots = EvalCoeffsToSlots(raised, sp);

    // [선택] 스케일 맞춤
    // if (std::abs(ctSlots->GetScalingFactor() - expected) > tol) { ... }

    // (3) Slots -> Coeffs : 기존 공개 API
    auto ctOut = EvalSlotsToCoeffs(ctSlots);

    // (4) [선택] ModRaise/Refresh 등
    return ctOut;
}

*/

/*
Ciphertext<DCRTPoly> FHECKKSRNS::EvalBootstrapSHIP(ConstCiphertext<DCRTPoly>& ciphertext, uint32_t numIterations,
    uint32_t precision) const {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ciphertext->GetCryptoParameters());

    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping only supported with Hybrid key switching.");
    auto st = cryptoParams->GetScalingTechnique();
#if NATIVEINT == 128
    if (st == FLEXIBLEAUTO || st == FLEXIBLEAUTOEXT)
        OPENFHE_THROW("128-bit CKKS Bootstrapping only supported for FIXEDMANUAL and FIXEDAUTO.");
#endif
    if (numIterations != 1 && numIterations != 2)
        OPENFHE_THROW("CKKS Bootstrapping only supported for 1 or 2 iterations.");

#ifdef BOOTSTRAPTIMING
    TimeVar t;
    double timeEncode(0.0);
    double timeModReduce(0.0);
    double timeDecode(0.0);
#endif

    auto cc = ciphertext->GetCryptoContext();
    uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();
    auto initSizeQ = ciphertext->GetElements()[0].GetNumOfElements();
    uint32_t compositeDegree = cryptoParams->GetCompositeDegree();

    if (numIterations > 1) {
        // Step 1: Get the input.
        uint32_t powerOfTwoModulus = 1 << precision;

        // Step 2: Scale up by powerOfTwoModulus, and extend the modulus to powerOfTwoModulus * q.
        // Note that we extend the modulus implicitly without any code calls because the value always stays 0.
        // We multiply by powerOfTwoModulus, and leave the last CRT value to be 0 (mod powerOfTwoModulus).
        auto ctScaledUp = cc->EvalMultNoCheck(ciphertext, powerOfTwoModulus);
        ctScaledUp->SetLevel(L0 - ctScaledUp->GetElements()[0].GetNumOfElements());

        // Step 3: Bootstrap the initial ciphertext.
        auto ctInitialBootstrap = cc->EvalBootstrap(ciphertext, numIterations - 1, precision);
        cc->GetScheme()->ModReduceInternalInPlace(ctInitialBootstrap, compositeDegree);

        // Step 4: Scale up by powerOfTwoModulus.
        cc->GetScheme()->MultByIntegerInPlace(ctInitialBootstrap, powerOfTwoModulus);

        // Step 5: Mod-down to powerOfTwoModulus * q
        // We mod down, and leave the last CRT value to be 0 because it's divisible by powerOfTwoModulus.
        auto ctBootstrappedScaledDown = ctInitialBootstrap->Clone();
        auto bootstrappingSizeQ = ctBootstrappedScaledDown->GetElements()[0].GetNumOfElements();

        // If we start with more towers, than we obtain from bootstrapping, return the original ciphertext.
        if (bootstrappingSizeQ <= initSizeQ) {
            return ciphertext->Clone();
        }

        // TODO: YSP Can be removed for FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
        if (st != COMPOSITESCALINGAUTO && st != COMPOSITESCALINGMANUAL) {
            for (auto& cv : ctBootstrappedScaledDown->GetElements())
                cv.DropLastElements(bootstrappingSizeQ - initSizeQ);
            ctBootstrappedScaledDown->SetLevel(L0 - ctBootstrappedScaledDown->GetElements()[0].GetNumOfElements());
        }

        // Step 6 and 7: Calculate the bootstrapping error by subtracting the original ciphertext from the bootstrapped ciphertext. Mod down to q is done implicitly.
        auto ctBootstrappingError = cc->EvalSub(ctBootstrappedScaledDown, ctScaledUp);

        // Step 8: Bootstrap the error.
        auto ctBootstrappedError = cc->EvalBootstrap(ctBootstrappingError, 1, 0);
        cc->GetScheme()->ModReduceInternalInPlace(ctBootstrappedError, compositeDegree);

        // Step 9: Subtract the bootstrapped error from the initial bootstrap to get even lower error.
        auto finalCiphertext = cc->EvalSub(ctInitialBootstrap, ctBootstrappedError);

        // Step 10: Scale back down by powerOfTwoModulus to get the original message.
        cc->EvalMultInPlace(finalCiphertext, 1.0 / powerOfTwoModulus);
        return finalCiphertext;
    }

    uint32_t slots = ciphertext->GetSlots();

    auto elementParamsRaised = *(cryptoParams->GetElementParams());
    // For FLEXIBLEAUTOEXT we raised ciphertext does not include extra modulus
    // as it is multiplied by auxiliary plaintext
    if (st == FLEXIBLEAUTOEXT)
        elementParamsRaised.PopLastParam();

    auto paramsQ = elementParamsRaised.GetParams();
    uint32_t sizeQ = paramsQ.size();
    std::vector<NativeInteger> moduli(sizeQ);
    std::vector<NativeInteger> roots(sizeQ);
    for (uint32_t i = 0; i < sizeQ; ++i) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i] = paramsQ[i]->GetRootOfUnity();
    }
    auto elementParamsRaisedPtr =
        std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(cc->GetCyclotomicOrder(), moduli, roots);

    double qDouble = GetBigModulus(cryptoParams);
    double powP = std::pow(2, cryptoParams->GetPlaintextModulus());
    int32_t deg = std::round(std::log2(qDouble / powP));
#if NATIVEINT != 128
    if (deg > static_cast<int32_t>(m_correctionFactor) && st != COMPOSITESCALINGAUTO && st != COMPOSITESCALINGMANUAL) {
        OPENFHE_THROW("Degree [" + std::to_string(deg) + "] must be less than or equal to the correction factor [" +
            std::to_string(m_correctionFactor) + "].");
    }
#endif
    uint32_t correction = m_correctionFactor - deg;
    double post = std::pow(2, static_cast<double>(deg));

    // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
    double pre = (compositeDegree > 1) ? cryptoParams->GetScalingFactorReal(0) / qDouble : 1. / post;
    uint64_t scalar = std::llround(post);

    //------------------------------------------------------------------------------
    // RAISING THE MODULUS
    //------------------------------------------------------------------------------

    // In FLEXIBLEAUTO, raising the ciphertext to a larger number
    // of towers is a bit more complex, because we need to adjust
    // it's scaling factor to the one that corresponds to the level
    // it's being raised to.
    // Increasing the modulus

    auto raised = ciphertext->Clone();
    auto algo = cc->GetScheme();
    algo->ModReduceInternalInPlace(raised, compositeDegree * (raised->GetNoiseScaleDeg() - 1));
    AdjustCiphertext(raised, correction);

    uint32_t N = cc->GetRingDimension();
    if (compositeDegree > 1) {
        // RNS basis extension from level 0 RNS limbs to the raised RNS basis
        auto& ctxtDCRT = raised->GetElements();
        ExtendCiphertext(ctxtDCRT, *cc, elementParamsRaisedPtr);
        raised->SetLevel(L0 - ctxtDCRT[0].GetNumOfElements());
    }
    else {
        if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
            auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(raised->GetKeyTag());

            // transform from a denser secret to a sparser one
            raised = KeySwitchSparse(raised, evalKeyMap.at(2 * N - 4));

            // Only level 0 ciphertext used here. Other towers ignored to make CKKS bootstrapping faster.
            auto& ctxtDCRT = raised->GetElements();
            for (auto& poly : ctxtDCRT) {
                poly.SetFormat(COEFFICIENT);
                DCRTPoly temp(elementParamsRaisedPtr, COEFFICIENT);
                temp = poly.GetElementAtIndex(0);
                temp.SetFormat(EVALUATION);
                poly = std::move(temp);
            }
            raised->SetLevel(L0 - ctxtDCRT[0].GetNumOfElements());

            // go back to a denser secret
            algo->KeySwitchInPlace(raised, evalKeyMap.at(2 * N - 2));
        }
        else {
            // Only level 0 ciphertext used here. Other towers ignored to make CKKS bootstrapping faster.
            auto& ctxtDCRT = raised->GetElements();
            for (auto& poly : ctxtDCRT) {
                poly.SetFormat(COEFFICIENT);
                DCRTPoly temp(elementParamsRaisedPtr, COEFFICIENT);
                temp = poly.GetElementAtIndex(0);
                temp.SetFormat(EVALUATION);
                poly = std::move(temp);
            }
            raised->SetLevel(L0 - ctxtDCRT[0].GetNumOfElements());
        }
    }

#ifdef BOOTSTRAPTIMING
    std::cerr << "\nNumber of levels at the beginning of bootstrapping: "
        << raised->GetElements()[0].GetNumOfElements() - 1 << std::endl;
#endif

    //------------------------------------------------------------------------------
    // SETTING PARAMETERS FOR APPROXIMATE MODULAR REDUCTION
    //------------------------------------------------------------------------------

    // Coefficients of the Chebyshev series interpolating 1/(2 Pi) Sin(2 Pi K x)
    std::vector<double> coefficients;
    double k = 0;

    if (cryptoParams->GetSecretKeyDist() == SPARSE_TERNARY) {
        coefficients = g_coefficientsSparse;
        // k = K_SPARSE;
        k = 1.0;  // do not divide by k as we already did it during precomputation
    }
    else if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
        coefficients = g_coefficientsSparseEncapsulated;
        k = 1.0;  // do not divide by k as we already did it during precomputation
    }
    else {
        // For larger composite degrees, larger K used to achieve a reasonable probability of failure
        if ((compositeDegree == 1) || ((compositeDegree == 2) && (N < (1 << 17)))) {
            coefficients = g_coefficientsUniform;
            k = K_UNIFORM;
        }
        else {
            coefficients = g_coefficientsUniformExt;
            k = K_UNIFORMEXT;
        }
    }

    cc->EvalMultInPlace(raised, pre * (1.0 / (k * N)));

    // no linear transformations are needed for Chebyshev series as the range has been normalized to [-1,1]
    double coeffLowerBound = -1;
    double coeffUpperBound = 1;

    auto& p = GetBootPrecom(slots);

    bool isLTBootstrap =
        (p.m_paramsEnc[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1) && (p.m_paramsDec[CKKS_BOOT_PARAMS::LEVEL_BUDGET] == 1);

    Ciphertext<DCRTPoly> ctxtDec;
    if (slots == cc->GetCyclotomicOrder() / 4) {
        //------------------------------------------------------------------------------
        // FULLY PACKED CASE
        //------------------------------------------------------------------------------

#ifdef BOOTSTRAPTIMING
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running CoeffToSlot
        //------------------------------------------------------------------------------

        // need to call internal modular reduction so it also works for FLEXIBLEAUTO
        algo->ModReduceInternalInPlace(raised, compositeDegree);

        // only one linear transform is needed as the other one can be derived
        auto ctxtEnc =
            (isLTBootstrap) ? EvalLinearTransform(p.m_U0hatTPre, raised) : EvalCoeffsToSlots(p.m_U0hatTPreFFT, raised);

        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ctxtEnc->GetKeyTag());
        auto conj = Conjugate(ctxtEnc, evalKeyMap);
        auto ctxtEncI = cc->EvalSub(ctxtEnc, conj);
        algo->MultByMonomialInPlace(ctxtEncI, 3 * slots);
        cc->EvalAddInPlaceNoCheck(ctxtEnc, conj);

        if (st == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                cc->ModReduceInPlace(ctxtEnc);
                cc->ModReduceInPlace(ctxtEncI);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
                algo->ModReduceInternalInPlace(ctxtEncI, compositeDegree);
            }
        }

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------

        // Evaluate Chebyshev series for the sine wave
        ctxtEnc = cc->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);
        ctxtEncI = cc->EvalChebyshevSeries(ctxtEncI, coefficients, coeffLowerBound, coeffUpperBound);

        // Double-angle iterations
        if (st != FIXEDMANUAL) {
            algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
            algo->ModReduceInternalInPlace(ctxtEncI, compositeDegree);
        }
        uint32_t numIter;
        if (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY)
            numIter = R_UNIFORM;
        else
            numIter = R_SPARSE;
        ApplyDoubleAngleIterations(ctxtEnc, numIter);
        ApplyDoubleAngleIterations(ctxtEncI, numIter);

        algo->MultByMonomialInPlace(ctxtEncI, slots);
        cc->EvalAddInPlaceNoCheck(ctxtEnc, ctxtEncI);

        if (st != COMPOSITESCALINGAUTO && st != COMPOSITESCALINGMANUAL) {
            // scale the message back up after Chebyshev interpolation
            algo->MultByIntegerInPlace(ctxtEnc, scalar);
        }

#ifdef BOOTSTRAPTIMING
        timeModReduce = TOC(t);
        std::cerr << "Approximate modular reduction time: " << timeModReduce / 1000.0 << " s" << std::endl;
        // Running SlotToCoeff
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running SlotToCoeff
        //------------------------------------------------------------------------------

        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (st != FIXEDMANUAL)
            algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);

        // Only one linear transform is needed
        ctxtDec = (isLTBootstrap) ? EvalLinearTransform(p.m_U0Pre, ctxtEnc) : EvalSlotsToCoeffs(p.m_U0PreFFT, ctxtEnc);
    }
    else {
        //------------------------------------------------------------------------------
        // SPARSELY PACKED CASE
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // Running PartialSum
        //------------------------------------------------------------------------------

        for (uint32_t j = 1; j < N / (2 * slots); j <<= 1)
            cc->EvalAddInPlaceNoCheck(raised, cc->EvalRotate(raised, j * slots));

#ifdef BOOTSTRAPTIMING
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running CoeffsToSlots
        //------------------------------------------------------------------------------

        algo->ModReduceInternalInPlace(raised, compositeDegree);

        auto ctxtEnc =
            (isLTBootstrap) ? EvalLinearTransform(p.m_U0hatTPre, raised) : EvalCoeffsToSlots(p.m_U0hatTPreFFT, raised);

        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(ctxtEnc->GetKeyTag());
        auto conj = Conjugate(ctxtEnc, evalKeyMap);
        cc->EvalAddInPlaceNoCheck(ctxtEnc, conj);

        if (st == FIXEDMANUAL) {
            while (ctxtEnc->GetNoiseScaleDeg() > 1) {
                cc->ModReduceInPlace(ctxtEnc);
            }
        }
        else {
            if (ctxtEnc->GetNoiseScaleDeg() == 2) {
                algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
            }
        }

#ifdef BOOTSTRAPTIMING
        timeEncode = TOC(t);
        std::cerr << "\nEncoding time: " << timeEncode / 1000.0 << " s" << std::endl;
        // Running Approximate Mod Reduction
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running Approximate Mod Reduction
        //------------------------------------------------------------------------------

        // Evaluate Chebyshev series for the sine wave
        ctxtEnc = cc->EvalChebyshevSeries(ctxtEnc, coefficients, coeffLowerBound, coeffUpperBound);

        // Double-angle iterations
        if (st != FIXEDMANUAL)
            algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);
        uint32_t numIter = (cryptoParams->GetSecretKeyDist() == UNIFORM_TERNARY) ? R_UNIFORM : R_SPARSE;
        ApplyDoubleAngleIterations(ctxtEnc, numIter);

        // TODO: YSP Can be extended to FLEXIBLE* scaling techniques as well as the closeness of 2^p to moduli is no longer needed
        if (st != COMPOSITESCALINGAUTO && st != COMPOSITESCALINGMANUAL) {
            // scale the message back up after Chebyshev interpolation
            algo->MultByIntegerInPlace(ctxtEnc, scalar);
        }

#ifdef BOOTSTRAPTIMING
        timeModReduce = TOC(t);
        std::cerr << "Approximate modular reduction time: " << timeModReduce / 1000.0 << " s" << std::endl;
        // Running SlotToCoeff
        TIC(t);
#endif

        //------------------------------------------------------------------------------
        // Running SlotsToCoeffs
        //------------------------------------------------------------------------------

        // In the case of FLEXIBLEAUTO, we need one extra tower
        // TODO: See if we can remove the extra level in FLEXIBLEAUTO
        if (st != FIXEDMANUAL)
            algo->ModReduceInternalInPlace(ctxtEnc, compositeDegree);

        // linear transform for decoding
        ctxtDec = (isLTBootstrap) ? EvalLinearTransform(p.m_U0Pre, ctxtEnc) : EvalSlotsToCoeffs(p.m_U0PreFFT, ctxtEnc);
        cc->EvalAddInPlaceNoCheck(ctxtDec, cc->EvalRotate(ctxtDec, slots));
    }

#if NATIVEINT != 128
    // 64-bit only: scale back the message to its original scale.
    uint64_t corFactor = static_cast<uint64_t>(1) << std::llround(correction);
    algo->MultByIntegerInPlace(ctxtDec, corFactor);
#endif

#ifdef BOOTSTRAPTIMING
    timeDecode = TOC(t);

    std::cout << "Decoding time: " << timeDecode / 1000.0 << " s" << std::endl;
#endif

    // If we start with more towers, than we obtain from bootstrapping, return the original ciphertext.
    if (ctxtDec->GetElements()[0].GetNumOfElements() <= initSizeQ)
        return ciphertext->Clone();
    return ctxtDec;
}

*/


void run_minimal() {
    // ===== 1) 파라미터 =====
    CCParams<CryptoContextCKKSRNS> p;

    const uint32_t ringDim = 1u << 14;   // 16384
    const uint32_t slots = ringDim / 2; // 8192

    p.SetRingDim(ringDim);
    p.SetBatchSize(slots);

    // 스케일/모듈러스 설정
    p.SetScalingTechnique(FLEXIBLEAUTO);
    p.SetScalingModSize(59);
    p.SetFirstModSize(60);

    // 키스위칭(부트스트랩 필수)
    p.SetKeySwitchTechnique(HYBRID);
    p.SetNumLargeDigits(3);

    // 보안레벨은 직접 파라미터 지정 모드
    p.SetSecurityLevel(HEStd_NotSet);

    // ---- 부트스트랩에 필요한 깊이 계산 ----
    vector<uint32_t> levelBudget = { 4, 4 }; // encode/decode FFT 예산
    vector<uint32_t> dim = { 0, 0 }; // auto baby/giant-step
    const SecretKeyDist skDist = UNIFORM_TERNARY;
    const uint32_t pipelineStages = 9;

    usint bootstrapDepth = FHECKKSRNS::GetBootstrapDepth(pipelineStages, levelBudget, skDist);
    usint levelsKeep = 22;                  // 부트스트랩 후 남길 레벨
    usint totalDepth = bootstrapDepth + levelsKeep;
    p.SetMultiplicativeDepth(totalDepth);

    // ===== 2) 컨텍스트 생성 & 기능 활성화 =====
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

    // ===== 3) 부트스트랩 사전계산 =====
    cc->EvalBootstrapSetup(levelBudget, dim, slots);
    //BootstrapShipParams sp;
    //EvalBootstrapSetupSHIP(sp);

    // ===== 4) 키 생성 =====
    auto kp = cc->KeyGen();
    cc->EvalMultKeyGen(kp.secretKey);
    cc->EvalBootstrapKeyGen(kp.secretKey, slots);

    // ===== 5) 입력 생성/암호화 =====
    vector<complex<double>> x(slots);
    for (size_t i = 0; i < slots; ++i)
        x[i] = { std::sin(i * 0.01), std::cos(i * 0.01) };

    auto pt = cc->MakeCKKSPackedPlaintext(x);
    auto ct0 = cc->Encrypt(kp.publicKey, pt);

    // ===== 6) (의도적) 레벨 소모 체인 =====
    auto ct_work = ct0;
    int mults_before = 40;   // 곱셈 20회로 레벨 소모
    for (int i = 0; i < mults_before; ++i) {
        ct_work = cc->EvalMult(ct_work, ct_work);      // 제곱
        // FLEXIBLEAUTO면 내부적으로 리스케일 진행
    }
    
    cout << "[*] Before bootstrap: Level=" << ct_work->GetLevel()
        << ", NoiseScaleDeg=" << ct_work->GetNoiseScaleDeg()
        << ", ScalingFactor=" << ct_work->GetScalingFactor() << "\n";

    // ===== 7) 부트스트랩 =====
    cout << "[*] Bootstrapping...\n";
    auto t0 = chrono::high_resolution_clock::now();
    //auto ct_bs = cc->EvalBootstrap(ct_work);
    auto ct_bs = EvalBootstrapStd(cc, ct_work);
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


    // ===== 8) 복호 & 오차(부트스트랩 전후 동일성) =====
    // 부트스트랩은 근사적 항등연산이므로 ct_work ≈ ct_bs
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

    // ===== 9) (핵심) 부트스트랩 이후 추가 연산이 다시 가능함을 확인 =====
    // 같은 강도의 연산(곱셈 6회 + 회전)을 부트스트랩 이후에도 한 번 더 시도
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

    // 오차도 한 번 더 체크(선택)
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


    /*
    chrono::high_resolution_clock::time_point clock_start, clock_end;
    updateGlobal();


    cout << "1. Setup Started!\n" << endl;

    auto PSsk = PSskGen(PSparam);
    auto PSpk = PSpkGen(PSparam, PSsk);

    CCParams<CryptoContextCKKSRNS> CKKSparam;
    initCKKSparam(CKKSparam);

    CryptoContext<DCRTPoly> context = GenCryptoContext(CKKSparam);
    enable(context);

    auto keyPair = context->KeyGen();
    context->EvalMultKeyGen(keyPair.secretKey);
    context->EvalRotateKeyGen(keyPair.secretKey, { 1, b_tilde1 });

    auto PSsk_enc = encryptPSsk_ckks(context, keyPair.publicKey, PSsk);

    //CKKSparam.SetMultiplicativeDepth(MultiplicativeDepth_comp);
    
    CryptoContext<DCRTPoly> context_comp = GenCryptoContext(CKKSparam);
    enable(context_comp);

    std::vector<uint32_t> dim1 = { 8, 8 }; 
    std::vector<uint32_t> levelBudget = { 6, 6 };
    context_comp->EvalBootstrapSetup(levelBudget, dim1, degree);

    auto keyPair_comp = context_comp->KeyGen();
    context_comp->EvalMultKeyGen(keyPair_comp.secretKey);

    std::vector<int32_t> rotIdx;
    for (int j = 1; j < degree_half / numrow_po2; j *= 2)
        rotIdx.push_back(static_cast<int32_t>(numrow_po2 * j));

    rotIdx.push_back(static_cast<int32_t>(degree_half));
    rotIdx.push_back(static_cast<int32_t>(degree_half / 2));
    rotIdx.push_back(1);
    rotIdx.push_back(2);

    std::sort(rotIdx.begin(), rotIdx.end());
    rotIdx.erase(std::unique(rotIdx.begin(), rotIdx.end()), rotIdx.end());
    context_comp->EvalRotateKeyGen(keyPair_comp.secretKey, rotIdx);
    context_comp->EvalBootstrapKeyGen(keyPair_comp.secretKey, degree);


    //CKKSparam.SetRingDim(degree_trace);
    CryptoContext<DCRTPoly> context_trace = GenCryptoContext(CKKSparam);
    enable(context_trace);
    auto keyPair_trace = context_trace->KeyGen();
    
    liftsk2(keyPair_comp, keyPair_trace);
    
    auto swk = context->GetScheme()->KeySwitchGen(keyPair.secretKey, keyPair_comp.secretKey);

    vector<int> step_comp = {1, b_tilde2, degree_half, degree_trace_half};
    for (int j = 1; j < degree_half / numrow_po2; j*=2) { 
        step_comp.push_back(numrow_po2 * j);
    }
    context_comp->EvalRotateKeyGen(keyPair_comp.secretKey, step_comp);

    updateTraceInfo_ckks(context_comp, context_trace, keyPair_comp, keyPair_trace);



    cout << "Setup Finished!\n" << endl;

    printParam(context, context_trace);
    saveKeys(context, context_comp, PSsk_enc, swk);

    
    cout << "2. Payload & Signal Simulation Started!\n" << endl;

    vector<vector<uint32_t>> payloads(num_transaction,vector<uint32_t>(payload_len));
    simulatePayloads(payloads);

    vector<int> pertinentIdx;
    sampleIdx(pertinentIdx);

    vector<vector<uint64_t>> signals_a(num_transaction,vector<uint64_t>(PSparam.n));
    vector<vector<uint64_t>> signals_b(num_transaction,vector<uint64_t>(PSparam.ell));
    simulateSignals(signals_a, signals_b, pertinentIdx, PSpk);
    
    std::vector<std::vector<std::complex<double>>> encoded_a(num_transaction);
    std::vector<std::vector<std::complex<double>>> encoded_a_minus(num_transaction);
    std::vector<std::vector<std::complex<double>>> encoded_b(num_transaction);

    for (size_t t = 0; t < static_cast<size_t>(num_transaction); ++t) {
        encoded_a[t].reserve(signals_a[t].size());
        encoded_a_minus[t].reserve(signals_a[t].size());
        encoded_b[t].reserve(signals_b[t].size());

        for (uint64_t val : signals_a[t]) {
            int m = static_cast<int>(symmetric_embed(val, ptxt_modulus));  // [-p/2, p/2)
            encoded_a[t].push_back(encodeInt(m, ptxt_modulus));
        }
        for (uint64_t val : signals_a[t]) {
            int m = static_cast<int>(symmetric_embed(-1*val, ptxt_modulus));  // [-p/2, p/2)
            encoded_a_minus[t].push_back(encodeInt(m, ptxt_modulus));
        }
        for (uint64_t val : signals_b[t]) {
            int m = static_cast<int>(symmetric_embed(-1*val, ptxt_modulus));  // [-p/2, p/2)
            encoded_b[t].push_back(encodeInt(m, ptxt_modulus));
        }
    }

    cout << "Payload & Signal Simulation Finished!\n" << endl;

    cout << "3. Server-side Digesting Started!\n" << endl;
    
    cout << "\t 3-1. Detection Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();

    vector<Ciphertext<DCRTPoly>> PV(numctxt);

    detect(PV, encoded_a, encoded_a_minus, encoded_b, context, PSsk_enc, swk);
    
    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Detection Finished!" << endl;
    cout << "\t Detection time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    
    cout << "\t 3-2. Compression Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> digest;

    compress(digest, PV, payloads, context_comp, context_trace, keyPair_comp.publicKey, keyPair_trace.publicKey, keyPair_comp.secretKey);

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Compression Finished!" << endl;
    cout << "\t Compression time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    Serial::SerializeToFile("data/digest.txt", digest, SerType::BINARY);
    

    cout << "4. Client-side Decoding\n" << endl;
    clock_start = chrono::high_resolution_clock::now();
    

    // 1. PV를 모두 복호화하여 탐지된 인덱스 목록을 찾습니다.
    vector<int> detected_indices;
    for (int i = 0; i < numctxt; ++i) {
        std::cerr << std::boolalpha << (PV[i]->GetCryptoContext().get() == context.get()) << "\n";
        lbcrypto::Plaintext pt;
        context->Decrypt(keyPair.secretKey, PV[i], &pt);
        cout << "?" << endl;
        auto v = pt->GetRealPackedValue();

        int base = i * degree;
        int last = std::min(degree, num_transaction - base);
        for (int m = 0; m < last; ++m) {
            // CKKS 근사값이므로 0/1 판정은 이런 식으로
            if (std::abs(v[m] - 1.0) < std::abs(v[m])) {
                detected_indices.push_back(base + m);
            }
            // 또는 임계값 if (v[m] > 0.5) ...
        }
    }

    cout << 1 << endl;

    std::sort(detected_indices.begin(), detected_indices.end());
    detected_indices.erase(std::unique(detected_indices.begin(), detected_indices.end()), detected_indices.end());

    vector<int64_t> digest_idx(num_pertinent, 0);
    vector<vector<int64_t>> digest_payload(payload_len, vector<int64_t>(num_pertinent, 0));

    auto modpow = [&](int64_t a, int64_t e) -> int64_t {
        const int64_t p = ptxt_modulus;
        long long base = (a % p + p) % p;
        long long res = 1 % p;
        while (e > 0) {
            if (e & 1) res = (res * base) % p;
            base = (base * base) % p;
            e >>= 1;
        }
        return (int64_t)res;
        };

    for (int r_idx = 0; r_idx < num_pertinent; ++r_idx) {
        for (int pertinent_idx : detected_indices) {
            int64_t pw = modpow(pertinent_idx + 1, r_idx + 1);
            digest_idx[r_idx] = (digest_idx[r_idx] + pw) % ptxt_modulus;

            // ★ 여기가 핵심: payload 슬롯에 대해 한 번 더 루프
            for (int p_len = 0; p_len < payload_len; ++p_len) {
                int64_t val = payloads[pertinent_idx][p_len] % ptxt_modulus;
                digest_payload[p_len][r_idx] =
                    (digest_payload[p_len][r_idx] + (val * pw) % ptxt_modulus) % ptxt_modulus;
            }
        }
    }
    cout << 2 << endl;

    // 3. 기존 decode 로직의 일부를 호출합니다.
    NTL::ZZ p(ptxt_modulus);
    NTL::ZZ_p::init(p);
    vector<NTL::ZZ_p> idx_ZZ_p;
    vector<int> decodedIdx;
    vector<vector<uint32_t>> decodedPayload;

    decodeIdx(idx_ZZ_p, digest_idx);
    decodePayload(decodedPayload, digest_payload, idx_ZZ_p);
    cout << 3 << endl;

    decodedIdx.resize(idx_ZZ_p.size());
    for (size_t i = 0; i < decodedIdx.size(); i++) {
        decodedIdx[i] = NTL::conv<int>(idx_ZZ_p[i]) - 1;
    }
    
    vector<int> decodedIdx;
    vector<vector<uint32_t>> decodedPayload;

    decode(decodedIdx, decodedPayload, digest, context_trace, keyPair_trace.secretKey);
    
    clock_end = chrono::high_resolution_clock::now();
    cout << "Client-side Decoding Finished!" << endl;
    cout << "Client-side Decoding Time: " << chrono::duration<double, milli>(clock_end - clock_start).count() << "ms\n" << endl;

    checkResult(decodedIdx, decodedPayload, pertinentIdx, payloads);
    */
}