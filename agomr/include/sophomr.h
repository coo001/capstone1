#include "openfhe.h"

#include "signal.h"
#include "global.h"
#include "param.h"
#include "setup.h"
#include "detect.h"
#include "compress.h"
#include "decode.h"
#include <random>
#include <array>

using Matrix = std::vector<std::vector<uint64_t>>;
using Vector = std::vector<uint64_t>;
using PRGSeed = std::array<uint64_t, 4>;

struct Clue {
    Matrix M;
    PRGSeed seed;
};

struct ClueMeta {
    Clue clue;
    std::vector<int> pertinentIdx;
    std::vector<Vector> idList;
};


std::vector<Vector> generateAllIDs(size_t numUsers, size_t I, uint64_t t) {
    Matrix allIDs(numUsers, Vector(I));
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis(0, t - 1);

    for (size_t i = 0; i < numUsers; ++i)
        for (size_t j = 0; j < I; ++j)
            allIDs[i][j] = dis(gen);

    return allIDs;
}

std::vector<ClueMeta> AGOMR3_GenClue_Batch(
    const std::vector<Vector>& allIDs,
    uint64_t t, size_t I, size_t G, size_t Gprime, size_t n, size_t ell, int totalData
) {
    std::vector<ClueMeta> allClues(totalData);

    for (int dataIdx = 0; dataIdx < totalData; ++dataIdx) {
        //std::cout << "dataIdx : " << dataIdx << std::endl;
        ClueMeta meta;
        meta.clue.M.resize(Gprime, std::vector<uint64_t>(n + ell, 0));

        // Step 1: Randomly select G indices
        std::vector<int> sampled;
        sampleIdx(sampled);
        meta.pertinentIdx = sampled;
        meta.idList.resize(G);

        for (size_t g = 0; g < G; ++g)
            meta.idList[g] = allIDs[sampled[g]];

        // Step 2: PRG Seed and Z
        PRGSeed seed;
        for (auto& s : seed) s = rand();
        meta.clue.seed = seed;

        std::mt19937_64 prng(seed[0] ^ seed[1] ^ seed[2] ^ seed[3]);
        std::uniform_int_distribution<uint64_t> dist(0, t - 1);

        Matrix Z(I, std::vector<uint64_t>(Gprime));
        for (size_t i = 0; i < I; i++)
            for (size_t j = 0; j < Gprime; j++)
                Z[i][j] = dist(prng);

        // Step 3: Compute id₁,...,id_G × Z = ˜id ∈ Z_t^{G × G'}
        Matrix idTilde(G, std::vector<uint64_t>(Gprime, 0));
        for (size_t g = 0; g < G; g++) {
            for (size_t j = 0; j < Gprime; j++) {
                for (size_t i = 0; i < I; i++) {
                    idTilde[g][j] += meta.idList[g][i] * Z[i][j];
                }
            }
        }

        // Step 4: Dummy clue vectors
        Matrix clueMatrix(G, std::vector<uint64_t>(n + ell, 1));
        for (size_t j = 0; j < n + ell; j++) {
            for (size_t k = 0; k < Gprime; k++) {
                meta.clue.M[k][j] = 0;
                for (size_t g = 0; g < G; g++) {
                    meta.clue.M[k][j] += idTilde[g][k] * clueMatrix[g][j];
                }
            }
        }

        allClues[dataIdx] = meta;
    }

    return allClues;
}

Matrix generatePRGMatrix(const PRGSeed& seed, size_t I, size_t Gprime, uint64_t t) {
    std::mt19937_64 prng(seed[0] ^ seed[1] ^ seed[2] ^ seed[3]);
    std::uniform_int_distribution<uint64_t> dist(0, t - 1);
    Matrix Z(I, std::vector<uint64_t>(Gprime));
    for (size_t i = 0; i < I; ++i)
        for (size_t j = 0; j < Gprime; ++j)
            Z[i][j] = dist(prng);
    return Z;
}

std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
RetrieveFromPlainID(const std::vector<uint64_t>& id, const Matrix& M, const PRGSeed& seed) {
    size_t Gprime = M.size();               // 행 개수
    size_t n_plus_ell = M[0].size();        // 열 개수
    size_t I = id.size();

    // Step 1: PRG(seed)로 Z ∈ Z_t^{I x G'} 생성
    Matrix Z = generatePRGMatrix(seed, I, Gprime, ptxt_modulus);
    
    // Step 2: id × Z → idZ ∈ Z_t^{G'}
    std::vector<uint64_t> idZ(Gprime, 0);
    for (size_t j = 0; j < Gprime; j++) {
        for (size_t i = 0; i < I; i++) {
            idZ[j] = (idZ[j] + id[i] * Z[i][j]) % ptxt_modulus;
        }
    }

    // Step 3: idZ × M → signal = idZ · M ∈ Z_t^{n + ell}
    std::vector<uint64_t> signal(n_plus_ell, 0);
    for (size_t j = 0; j < n_plus_ell; j++) {
        for (size_t i = 0; i < Gprime; i++) {
            signal[j] = (signal[j] + idZ[i] * M[i][j]) % ptxt_modulus;
        }
    }

    // Step 4: signal → sig_a, sig_b 분리
    size_t n = n_plus_ell - PSparam.ell;
    std::vector<uint64_t> sig_a(PSparam.n, 0);
    std::vector<uint64_t> sig_b(PSparam.ell, 0);

    for (size_t i = 0; i < static_cast<size_t>(PSparam.n); i++) {
        sig_a[i] = signal[i];
    }
    for (size_t i = 0; i < static_cast<size_t>(PSparam.ell); i++) {
        sig_b[i] = signal[n + i];
    }

    return { sig_a, sig_b };
}


std::pair<std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>> retrieve(
    const std::vector<uint64_t>& id,                   // id ∈ Z_t^I
    const std::vector<ClueMeta>& clue
) {
    std::vector<std::vector<uint64_t>> sig_a(num_transaction, std::vector<uint64_t>(PSparam.n, 0));
    std::vector<std::vector<uint64_t>> sig_b(num_transaction, std::vector<uint64_t>(PSparam.ell, 0));
    int t = 0;
    for (const auto& meta : clue) {
        const Clue& c = meta.clue;

        auto [signal_a, signal_b] = RetrieveFromPlainID(id, c.M, c.seed);

        for (size_t i = 0; i < static_cast<size_t>(PSparam.n); i++) {
            sig_a[t][i] += signal_a[i];
        }
        for (size_t i = 0; i < static_cast<size_t>(PSparam.ell); i++) {
            sig_b[t][i] += signal_b[i];
        }
        t++;
    }

    return { sig_a, sig_b };
}

void check(const std::vector<int>& decodedIdx,
    const std::vector<ClueMeta>& clueList,
    const Vector& my_id) {
    std::cout << "\n===== [Check Result] =====\n";

    for (int idx : decodedIdx) {
        const ClueMeta& meta = clueList[idx];
        bool match = false;

        for (const auto& id : meta.idList) {
            if (id == my_id) {
                match = true;
                break;
            }
        }

        if (match) {
            std::cout << "  [✓] idx " << idx << ": my_id is in pertinent id list\n";
        }
        else {
            std::cout << "  [✗] idx " << idx << ": my_id NOT in pertinent id list\n";
        }
    }

    std::cout << "==========================\n";
}


void sophomr() {
    using namespace std;
    using namespace lbcrypto;

    chrono::high_resolution_clock::time_point clock_start, clock_end;
    updateGlobal();

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

    cout << "2. Clue Generate Started!\n" << endl;

    vector<vector<uint32_t>> payloads(num_transaction,vector<uint32_t>(payload_len));
    simulatePayloads(payloads);

    std::vector<Vector> allIDs = generateAllIDs(num_transaction, 20, ptxt_modulus);
    std::vector<ClueMeta> clue = AGOMR3_GenClue_Batch(allIDs, ptxt_modulus, 20, num_pertinent, 16, PSparam.n, PSparam.ell, num_transaction);
   
    cout << "Clue Generate Finished!\n" << endl;



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

    // 평문 id 사용
    Vector my_id = allIDs[100];
    auto [sig_a, sig_b] = retrieve(my_id, clue);
    std::vector<std::vector<uint64_t>> signals_a = sig_a;
    std::vector<std::vector<uint64_t>> signals_b = sig_b;

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Clue Retrieve time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    clock_start = chrono::high_resolution_clock::now();
    detect(PV, signals_a, signals_b, context, PSsk_enc, swk, keyPair_trace.secretKey);

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

    check(decodedIdx, clue, my_id);
    //checkResult(decodedIdx, decodedPayload, pertinentIdx, payloads);
}