#include "openfhe.h"

#include "signal.h"
#include "global.h"
#include "param.h"
#include "setup.h"
#include "detect.h"
#include "compress.h"
#include "decode.h"

void generate_coeff_rangeCheck(int r, uint32_t p) {
    std::vector<uint64_t> coeffs = { 1 }; // 초기 다항식: 상수항 1

    for (int i = 0; i <= r; ++i) {
        std::vector<uint64_t> new_factor = {
            static_cast<uint64_t>((p - (i * i % p)) % p), 0, 1
        };
        std::vector<uint64_t> new_coeffs(coeffs.size() + 2, 0);

        for (size_t deg1 = 0; deg1 < coeffs.size(); ++deg1) {
            for (size_t deg2 = 0; deg2 < new_factor.size(); ++deg2) {
                new_coeffs[deg1 + deg2] += (coeffs[deg1] * new_factor[deg2]) % p;
                new_coeffs[deg1 + deg2] %= p;
            }
        }

        coeffs = std::move(new_coeffs);
    }

    // 짝수 차수 계수 추출
    std::vector<uint32_t> filtered;
    for (size_t i = 0; i < coeffs.size(); i += 2) {
        filtered.push_back(static_cast<uint32_t>(coeffs[i]));
    }

    // 전역 변수에 8개 열로 나눠서 저장
    size_t k = 8;
    size_t rows = (filtered.size() + k - 1) / k;
    coeff_rangeCheck.assign(rows, std::vector<uint32_t>(k, 0));

    for (size_t i = 0; i < filtered.size(); ++i) {
        coeff_rangeCheck[i / k][i % k] = filtered[i];
    }
}

void sophomr() {
    using namespace std;
    using namespace lbcrypto;

    chrono::high_resolution_clock::time_point clock_start, clock_end;
    updateGlobal();
    generate_coeff_rangeCheck(range_r, ptxt_modulus);
    std::cout << coeff_rangeCheck.size() << std::endl;
    std::cout << coeff_rangeCheck[0].size() << std::endl;

    /*////////////////////////////////////////////////////////////////


            1. Setup


    */////////////////////////////////////////////////////////////////

    cout << "1. Setup Started!\n" << endl;

    /*
                1-1. PS Setup
    */

    auto PSsk = PSskGen(PSparam);
    auto PSpk = PSpkGen(PSparam, PSsk);

    /*
                1-2. HE Setup (Top Level)
    */

    CCParams<CryptoContextBFVRNS> BFVparam;
    initBFVparam(BFVparam);

    CryptoContext<DCRTPoly> context = GenCryptoContext(BFVparam);
    enable(context);

    auto keyPair = context->KeyGen();
    context->EvalMultKeyGen(keyPair.secretKey);
    context->EvalRotateKeyGen(keyPair.secretKey, {1, b_tilde1});

    auto PSsk_enc = encryptPSsk(context, keyPair.secretKey, PSsk); 

    /*
                1-3. HE Setup (Compression & Ring-Switching)
    */

    BFVparam.SetNumLargeDigits(NumLargeDigits_comp);
    BFVparam.SetMultiplicativeDepth(MultiplicativeDepth_comp);
    CryptoContext<DCRTPoly> context_comp = GenCryptoContext(BFVparam);
    enable(context_comp);
    auto keyPair_comp = context_comp->KeyGen();

	BFVparam.SetRingDim(degree_trace);
    CryptoContext<DCRTPoly> context_trace = GenCryptoContext(BFVparam);
    enable(context_trace);
    auto keyPair_trace = context_trace->KeyGen();

	liftsk(keyPair_comp, keyPair_trace);

	auto swk = context->GetScheme()->KeySwitchGen(keyPair.secretKey, keyPair_comp.secretKey);

    vector<int> step_comp = {1, b_tilde2, degree_half, degree_trace_half};
    for (int j = 1; j < degree_half / numrow_po2; j*=2) { 
        step_comp.push_back(numrow_po2 * j);
    }
    context_comp->EvalRotateKeyGen(keyPair_comp.secretKey, step_comp);

    updateTraceInfo(context_comp, context_trace, keyPair_comp, keyPair_trace);


    cout << "Setup Finished!\n" << endl;

    printParam(context, context_trace);
    saveKeys(context, context_comp, PSsk_enc, swk);



    /*////////////////////////////////////////////////////////////////


            2. Payload & Signal Simulation


    */////////////////////////////////////////////////////////////////

    cout << "2. Payload & Signal Simulation Started!\n" << endl;

    vector<vector<uint32_t>> payloads(num_transaction,vector<uint32_t>(payload_len));
    simulatePayloads(payloads);

    vector<int> pertinentIdx;
    vector<int> pertinentIdx2;
    vector<int> pertinentIdx3;
    sampleIdx(pertinentIdx);
    sampleIdx(pertinentIdx2);
    sampleIdx(pertinentIdx3);

    vector<vector<uint64_t>> signals1_a(num_transaction,vector<uint64_t>(PSparam.n));
    vector<vector<uint64_t>> signals1_b(num_transaction,vector<uint64_t>(PSparam.ell));
    vector<vector<uint64_t>> signals2_a(num_transaction,vector<uint64_t>(PSparam.n));
    vector<vector<uint64_t>> signals2_b(num_transaction,vector<uint64_t>(PSparam.ell));
    vector<vector<uint64_t>> signals3_a(num_transaction, vector<uint64_t>(PSparam.n));
    vector<vector<uint64_t>> signals3_b(num_transaction, vector<uint64_t>(PSparam.ell));
    simulateSignals(signals1_a, signals1_b, pertinentIdx, PSpk);
    simulateSignals(signals2_a, signals2_b, pertinentIdx2, PSpk);
    simulateSignals(signals3_a, signals3_b, pertinentIdx3, PSpk);

    cout << "Payload & Signal Simulation Finished!\n" << endl;



    /*////////////////////////////////////////////////////////////////


            3. Server-side Digesting


    */////////////////////////////////////////////////////////////////

    cout << "3. Server-side Digesting Started!\n" << endl;

    /*
                3-1. Detection
    */

    cout << "\t 3-1. Detection Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();

    vector<Ciphertext<DCRTPoly>> PV(numctxt);

    //std::vector<std::vector<uint64_t>> signals_a = AddVectors(signals1_a, signals2_a);
    //std::vector<std::vector<uint64_t>> signals_b = AddVectors(signals1_b, signals2_b);
    //signals_a = AddVectors(signals_a, signals3_a);
    //signals_b = AddVectors(signals_b, signals3_b);

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Addvectors time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    clock_start = chrono::high_resolution_clock::now();
    detect(PV, signals1_a, signals1_b, context, PSsk_enc, swk, keyPair_trace.secretKey);

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Detection Finished!" << endl;
    cout << "\t Detection time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    /*
                3-2. Compression
    */

    cout << "\t 3-2. Compression Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();

    Ciphertext<DCRTPoly> digest;

    compress(digest, PV, payloads, context_comp, context_trace, keyPair_comp.publicKey, keyPair_trace.publicKey, keyPair_comp.secretKey);

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Compression Finished!" << endl;
    cout << "\t Compression time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    Serial::SerializeToFile("data/digest.txt", digest, SerType::BINARY);



    /*////////////////////////////////////////////////////////////////


            4. Client-side Decoding


    */////////////////////////////////////////////////////////////////

    cout << "4. Client-side Decoding\n" << endl;
    clock_start = chrono::high_resolution_clock::now();

    vector<int> decodedIdx;
    vector<vector<uint32_t>> decodedPayload;

    decode(decodedIdx, decodedPayload, digest, context_trace, keyPair_trace.secretKey);

    clock_end = chrono::high_resolution_clock::now();
    cout << "Client-side Decoding Finished!" << endl;
    cout << "Client-side Decoding Time: " << chrono::duration<double, milli>(clock_end - clock_start).count() << "ms\n" << endl;

    checkResult(decodedIdx, decodedPayload, pertinentIdx, pertinentIdx2, pertinentIdx3, payloads);
}