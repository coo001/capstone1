#include "openfhe.h"

#include "signal.h"
#include "global.h"
#include "param.h"
#include "setup.h"
#include "detect.h"
#include "compress.h"
#include "decode.h"
#include "enco.h"

void sophomr() {
    using namespace std;
    using namespace lbcrypto;

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
    context->EvalRotateKeyGen(keyPair.secretKey, {1, b_tilde1});

    auto PSsk_enc = encryptPSsk_ckks(context, keyPair.publicKey, PSsk);

    //CKKSparam.SetMultiplicativeDepth(MultiplicativeDepth_comp);
    cout << 1 << endl;

    CryptoContext<DCRTPoly> context_comp = GenCryptoContext(CKKSparam);
    enable(context_comp);
    auto keyPair_comp = context_comp->KeyGen();
    context_comp->EvalMultKeyGen(keyPair_comp.secretKey);
    context_comp->EvalRotateKeyGen(keyPair_comp.secretKey, { 1,2 });

    //CKKSparam.SetRingDim(degree_trace);
    CryptoContext<DCRTPoly> context_trace = GenCryptoContext(CKKSparam);
    enable(context_trace);
    auto keyPair_trace = context_trace->KeyGen();
    cout << 2 << endl;
	liftsk2(keyPair_comp, keyPair_trace);
    cout << 3 << endl;
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

    cout << "Payload & Signal Simulation Finished!\n" << endl;

    cout << "3. Server-side Digesting Started!\n" << endl;
    
    cout << "\t 3-1. Detection Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();

    vector<Ciphertext<DCRTPoly>> PV(numctxt);

    detect(PV, signals_a, signals_b, context, PSsk_enc, swk); 
    uint32_t numTowers1 = PV[0]->GetElements()[0].GetNumOfElements();
    cout << "num towers : " << numTowers1 << endl;
    /*
    for (size_t i = 0; i < PV.size(); ++i) {
        std::vector<double> vals(degree, i); // i·Î Ã¤¿î º¤ÅÍ
        auto ptxt = context->MakeCKKSPackedPlaintext(vals);
        PV[i] = context->Encrypt(keyPair.publicKey, ptxt);
    }
    */
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

    vector<int> decodedIdx;
    vector<vector<uint32_t>> decodedPayload;

    decode(decodedIdx, decodedPayload, digest, context_trace, keyPair_trace.secretKey);

    clock_end = chrono::high_resolution_clock::now();
    cout << "Client-side Decoding Finished!" << endl;
    cout << "Client-side Decoding Time: " << chrono::duration<double, milli>(clock_end - clock_start).count() << "ms\n" << endl;

    checkResult(decodedIdx, decodedPayload, pertinentIdx, payloads);
}