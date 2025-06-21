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

Matrix generatePRGMatrix(const PRGSeed& seed, size_t I, size_t Gprime, uint64_t t) {
    std::mt19937_64 prng(seed[0] ^ seed[1] ^ seed[2] ^ seed[3]);
    std::uniform_int_distribution<uint64_t> dist(0, t - 1);
    Matrix Z(I, std::vector<uint64_t>(Gprime));
    for (size_t i = 0; i < I; ++i)
        for (size_t j = 0; j < Gprime; ++j)
            Z[i][j] = dist(prng);
    return Z;
}

static uint64_t modInv(uint64_t a, uint64_t m) {
    int64_t t = 0, newt = 1;
    int64_t r = m, newr = a;
    while (newr != 0) {
        int64_t q = r / newr;
        std::tie(t, newt) = std::make_pair(newt, t - q * newt);
        std::tie(r, newr) = std::make_pair(newr, r - q * newr);
    }
    if (r > 1) throw std::runtime_error("modInv: not invertible");
    if (t < 0) t += m;
    return uint64_t(t);
}

static Matrix modInverseMatrix(const Matrix& A, uint64_t t) {
    size_t G = A.size();
    Matrix B = A; 
    Matrix I(G, std::vector<uint64_t>(G, 0));
    for (size_t i = 0; i < G; i++) I[i][i] = 1;

    for (size_t i = 0; i < G; i++) {
        size_t pivot = i;
        while (pivot < G && B[pivot][i] == 0) pivot++;
        if (pivot == G) throw std::runtime_error("modInverseMatrix: singular");
        if (pivot != i) {
            std::swap(B[pivot], B[i]);
            std::swap(I[pivot], I[i]);
        }
        uint64_t inv = modInv(B[i][i], t);
        for (size_t j = 0; j < G; j++) {
            B[i][j] = (B[i][j] * inv) % t;
            I[i][j] = (I[i][j] * inv) % t;
        }
        for (size_t r = 0; r < G; r++) if (r != i) {
            uint64_t factor = B[r][i];
            if (!factor) continue;
            for (size_t c = 0; c < G; c++) {
                B[r][c] = (B[r][c] + t - factor * B[i][c] % t) % t;
                I[r][c] = (I[r][c] + t - factor * I[i][c] % t) % t;
            }
        }
    }

    return I;
}

std::vector<ClueMeta> AGOMR3_GenClue_Batch(
    const std::vector<Vector>& allIDs,
    uint64_t t, size_t I, size_t G, size_t Gprime, size_t n, size_t ell, int totalData,
    const PSpk& PSpk
) {
    std::vector<ClueMeta> allClues(totalData);

    for (int dataIdx = 0; dataIdx < totalData; ++dataIdx) {
        if (dataIdx % 1000 == 0)
            std::cout << "dataIdx : " << dataIdx << std::endl;

        while (true) {
            try {
                ClueMeta meta;
                meta.clue.M.resize(Gprime, std::vector<uint64_t>(n + ell, 0));

                std::vector<int> sampled;
                sampleIdx(sampled, dataIdx < num_person);
                meta.pertinentIdx = sampled;
                meta.idList.resize(G);
                for (size_t g = 0; g < G; ++g)
                    meta.idList[g] = allIDs[sampled[g]];

                PRGSeed seed;
                for (auto& s : seed) s = rand();
                meta.clue.seed = seed;
                Matrix Z = generatePRGMatrix(seed, I, Gprime, t);

                Matrix idTilde(G, std::vector<uint64_t>(Gprime, 0));
                for (size_t g = 0; g < G; g++) {
                    for (size_t j = 0; j < Gprime; j++) {
                        for (size_t i = 0; i < I; i++) {
                            idTilde[g][j] += meta.idList[g][i] * Z[i][j];
                            idTilde[g][j] %= t;
                        }
                    }
                }

                Matrix clueMatrix(G, std::vector<uint64_t>(n + ell, 0));
                for (size_t g = 0; g < G; g++) {
                    Signal sig;
                    PSsignal(sig, PSpk, PSparam);

                    for (size_t j = 0; j < n; j++)
                        clueMatrix[g][j] = sig.a[j].ConvertToInt();

                    for (size_t j = 0; j < ell; j++)
                        clueMatrix[g][n + j] = sig.b[j].ConvertToInt();
                }

                Matrix A(G, std::vector<uint64_t>(G));
                for (size_t g = 0; g < G; g++)
                    for (size_t h = 0; h < G; h++)
                        A[g][h] = idTilde[g][h];

                // singular matrix 예외 가능성 있는 부분
                Matrix Ainv = modInverseMatrix(A, t);  // throws if singular

                for (size_t k = 0; k < G; k++) {
                    for (size_t j = 0; j < n + ell; j++) {
                        uint64_t acc = 0;
                        for (size_t h = 0; h < G; h++) {
                            acc = (acc + Ainv[k][h] * clueMatrix[h][j]) % t;
                        }
                        meta.clue.M[k][j] = acc;
                    }
                }

                allClues[dataIdx] = std::move(meta);
                break;  // 정상적으로 생성되었으면 반복 탈출

            }
            catch (const std::runtime_error& e) {
                if (std::string(e.what()).find("modInverseMatrix: singular") != std::string::npos) {
                    std::cerr << "[Warning] Singular matrix encountered at dataIdx " << dataIdx << ", retrying...\n";
                    continue;  // singular이면 다시 시도
                }
                else {
                    throw;  // 다른 오류는 재전파
                }
            }
        }
    }

    return allClues;
}

std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
RetrieveFromPlainID(const std::vector<uint64_t>& id, const Matrix& M, const PRGSeed& seed) {
    size_t Gprime = M.size();               
    size_t n_plus_ell = M[0].size();     
    size_t I = id.size();

    Matrix Z = generatePRGMatrix(seed, I, Gprime, ptxt_modulus);

    std::vector<uint64_t> idZ(Gprime, 0);
    for (size_t j = 0; j < Gprime; j++) {
        for (size_t i = 0; i < I; i++) {
            idZ[j] = (idZ[j] + id[i] * Z[i][j]) % ptxt_modulus;
        }
    }

    std::vector<uint64_t> signal(n_plus_ell, 0);
    for (size_t j = 0; j < n_plus_ell; j++) {
        for (size_t i = 0; i < Gprime; i++) {
            signal[j] = (signal[j] + idZ[i] * M[i][j]) % ptxt_modulus;
        }
    }

    std::vector<uint64_t> sig_a(PSparam.n, 0);
    std::vector<uint64_t> sig_b(PSparam.ell, 0);

    for (size_t i = 0; i < static_cast<size_t>(PSparam.n); i++) {
        sig_a[i] = signal[i];
    }
    for (size_t i = 0; i < static_cast<size_t>(PSparam.ell); i++) {
        sig_b[i] = signal[PSparam.n + i];
    }

    return { sig_a, sig_b };
}


void retrieve(
    std::vector<std::vector<uint64_t>>& sig_a,
    std::vector<std::vector<uint64_t>>& sig_b,
    const std::vector<uint64_t>& id,
    const std::vector<ClueMeta>& clue,
    const PSpk& PSpk
) {
    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(PSparam.q);

    int t = 0;
    for (const auto& meta : clue) {
        const Clue& c = meta.clue;

        auto [signal_a, signal_b] = RetrieveFromPlainID(id, c.M, c.seed);

        for (size_t i = 0; i < static_cast<size_t>(PSparam.n); i++) {
            sig_a[t][i] = signal_a[i];
        }
        for (size_t i = 0; i < static_cast<size_t>(PSparam.ell); i++) {
            sig_b[t][i] = signal_b[i];
        }
        t++;
    }
}

void check(const std::vector<int>& decodedIdx,
    const std::vector<ClueMeta>& clueList,
    const std::vector<Vector>& id_set) {
    std::cout << "\n===== [Check Result] =====\n";

    for (int idx : decodedIdx) {
        std::cout << "idx: " << idx << std::endl;
        /*const ClueMeta& meta = clueList[idx];
        bool match = true;

        for (const Vector& id : id_set) {
            bool found = false;
            for (const Vector& meta_id : meta.idList) {
                if (id == meta_id) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                match = false;
                break;
            }
        }

        if (match) {
            std::cout << "  [✓] idx " << idx << ": all my_ids are in pertinent id list\n";
        }
        else {
            std::cout << "  [✗] idx " << idx << ": some my_ids NOT in pertinent id list\n";
        }*/
    }

    std::cout << "==========================\n";
}

void printPertinentIndices(const std::vector<ClueMeta>& clueList) {
    std::cout << "\n===== [Pertinent Indices of ClueMeta Set] =====\n";
    for (size_t i = 0; i < clueList.size(); ++i) {
        std::cout << "ClueMeta[" << i << "] pertinentIdx: ";
        for (int idx : clueList[i].pertinentIdx) {
            std::cout << idx << " ";
        }
        std::cout << "\n";
    }
    std::cout << "===============================================\n";
}

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

lbcrypto::PublicKey<lbcrypto::DCRTPoly> MergePublicKeys(
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& cc,
    const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& pk1,
    const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& pk2) {
    using namespace lbcrypto;

    // new public key object to return
    auto pk_total = std::make_shared<PublicKeyImpl<DCRTPoly>>(cc);

    // 가져오기
    auto a1 = pk1->GetPublicElements()[0];
    auto b1 = pk1->GetPublicElements()[1];
    auto a2 = pk2->GetPublicElements()[0];
    auto b2 = pk2->GetPublicElements()[1];

    // 더하기
    auto a_total = a1 + a2;
    auto b_total = b1 + b2;

    std::vector<DCRTPoly> pubElements = { a_total, b_total };
    pk_total->SetPublicElements(pubElements);
    pk_total->SetKeyTag(pk1->GetKeyTag());

    return pk_total;
}

lbcrypto::KeyPair<lbcrypto::DCRTPoly> MergeKeyPairs(
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context,
    const lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair1,
    const lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair2) {
    using namespace std;
    using namespace lbcrypto;

    // 1. SecretKey 추출
    auto sk1 = keyPair1.secretKey->GetPrivateElement();
    auto sk2 = keyPair2.secretKey->GetPrivateElement();

    // 2. 둘 다 Format 변환 (NTT 도메인 등 연산 일관성 위해)
    sk1.SwitchFormat();
    sk2.SwitchFormat();

    // 3. 비밀키 합산
    auto sk_sum = sk1 + sk2;

    // 4. 다시 Format 돌리기
    sk_sum.SwitchFormat();

    // 5. 새로운 SecretKey 객체 생성
    auto sk_total = std::make_shared<PrivateKeyImpl<DCRTPoly>>(context);
    sk_total->SetPrivateElement(sk_sum);

    // 6. 새로운 PublicKey 생성
    auto pk_total = MergePublicKeys(context, keyPair1.publicKey, keyPair2.publicKey);

    // 7. 합산된 KeyPair 생성 및 반환
    KeyPair<DCRTPoly> keyPair_total;
    keyPair_total.secretKey = sk_total;
    keyPair_total.publicKey = pk_total;

    return keyPair_total;
}

void liftsk2(
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair_comp,
    const lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& fusedSecretKey   // 단일 trace 대신 fused
) {
    using namespace lbcrypto;

    // comp 쪽 비밀키 요소 가져오기
    auto sk_comp = keyPair_comp.secretKey->GetPrivateElement();
    sk_comp.SwitchFormat();

    // fused trace 키 하나로 리프팅
    auto sk_trace = fusedSecretKey->GetPrivateElement();
    sk_trace.SwitchFormat();

    for (size_t i = 0; i < sk_trace.GetNumOfElements(); i++) {
        auto temp_comp = sk_comp.GetElementAtIndex(i);
        auto temp_trace = sk_trace.GetElementAtIndex(i);
        for (size_t j = 0; j < sk_comp.GetElementAtIndex(i).GetLength(); j++) {
            temp_comp[j] = 0;
        }
        for (size_t j = 0; j < sk_trace.GetElementAtIndex(i).GetLength(); j++) {
            auto coeff_trace = temp_trace[j].ConvertToInt();
            if (coeff_trace == sk_trace.GetElementAtIndex(i).GetModulus() - 1) {
                temp_comp[j * dim_trace] = sk_comp.GetElementAtIndex(i).GetModulus() - 1;
            }
            else {
                temp_comp[j * dim_trace] = coeff_trace;
            }
        }
        sk_comp.SetElementAtIndex(i, temp_comp);
    }
    sk_comp.SwitchFormat();
    keyPair_comp.secretKey->SetPrivateElement(sk_comp);
}

void updateTraceInfo2(
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context_comp,
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context_trace,
    const lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair_comp,
    const lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair_trace,
    const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& pubkey
    ) {
    using namespace lbcrypto;

    // 1) comp 컨텍스트에서 인덱스용 plaintext 생성⋅암호화
    std::vector<int64_t> vec(degree, 0);
    for (int i = 0; i < degree_trace_half; i++) {
        vec[i] = i * dim_trace;
        vec[i + degree_half] = (i + degree_trace_half) * dim_trace;
    }
    auto ptxt_comp = context_comp->MakePackedPlaintext(vec);
    auto ctxt_comp = context_comp->Encrypt(keyPair_comp.secretKey, ptxt_comp);
    ctxt_comp = context_comp->Compress(ctxt_comp);

    // 2) trace 컨텍스트에서 빈 벡터 암호화(원래는 swap/shift 정보 추출용)
    std::vector<int64_t> vec_null(degree_trace);
    auto ptxt_null = context_trace->MakePackedPlaintext(vec_null);
    auto ctxt_trace = context_trace->Encrypt(pubkey, ptxt_null);
    ctxt_trace = context_trace->Compress(ctxt_trace);

    // 3) 내부 폴리 변환 후 comp→trace로 재배치
    auto poly_c = ctxt_comp->GetElements();
    auto poly_t = ctxt_trace->GetElements();
    for (int i = 0; i < 2; i++) {
        poly_c[i].SwitchFormat();
        poly_t[i].SwitchFormat();
    }
    for (int i = 0; i < 2; i++) {
        auto pc = poly_c[i].GetElementAtIndex(0);
        auto pt = poly_t[i].GetElementAtIndex(0);
        for (size_t k = 0; k < pt.GetLength(); k++) {
            pt[k] = pc[dim_trace * k].Mod(pt.GetModulus());
        }
        poly_t[i].SetElementAtIndex(0, pt);
    }
    for (int i = 0; i < 2; i++)
        poly_t[i].SwitchFormat();

    ctxt_trace->SetElements(poly_t);

    // 4) fusedSecretKey로 복호화 → swap/shift 정보 추출
    Plaintext ptxt_res;
    context_trace->Decrypt(keyPair_trace.secretKey, ctxt_trace, &ptxt_res);
    std::vector<int64_t> vec_res = ptxt_res->GetPackedValue();

    trace_swap = vec_res[0] / degree_trace_half;
    trace_shift = vec_res[0] % degree_trace_half;
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

    auto keyPair = context->KeyGen();;
    context->EvalMultKeyGen(keyPair.secretKey);
    context->EvalRotateKeyGen(keyPair.secretKey, { 1, b_tilde1 });

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

    auto keyPair_trace1 = context_trace->KeyGen();
    auto keyPair_trace2 = context_trace->MultipartyKeyGen(keyPair_trace1.publicKey, false, true);

    std::vector<lbcrypto::PrivateKey<lbcrypto::DCRTPoly>> secretKeys{
    keyPair_trace1.secretKey, keyPair_trace2.secretKey
    };

    auto keyPair_trace = context_trace->MultipartyKeyGen(secretKeys);
    auto pubkey = context_trace->MultiAddPubKeys(keyPair_trace1.publicKey, keyPair_trace2.publicKey, keyPair_trace2.publicKey->GetKeyTag());
    if (pubkey->GetCryptoContext() != context_trace) {
        std::cerr << "Warning: pubkey and context_trace do not match!" << std::endl;
    }

    liftsk(keyPair_comp, keyPair_trace);

    auto swk = context->GetScheme()->KeySwitchGen(keyPair.secretKey, keyPair_comp.secretKey);
    
    vector<int> step_comp = { 1, b_tilde2, degree_half, degree_trace_half };
    for (int j = 1; j < degree_half / numrow_po2; j *= 2) {
        step_comp.push_back(numrow_po2 * j);
    }
    context_comp->EvalRotateKeyGen(keyPair_comp.secretKey, step_comp);

    updateTraceInfo2(context_comp, context_trace, keyPair_comp, keyPair_trace, pubkey);

    cout << "Setup Finished!\n" << endl;

    printParam(context, context_trace);
    saveKeys(context, context_comp, PSsk_enc, swk);



    /*////////////////////////////////////////////////////////////////


            2. Payload & Signal Simulation


    */////////////////////////////////////////////////////////////////

    cout << "2. Clue Generate Started!\n" << endl;

    vector<vector<uint32_t>> payloads(num_transaction, vector<uint32_t>(payload_len));
    simulatePayloads(payloads);

    std::vector<Vector> allIDs = generateAllIDs(num_transaction, 10, ptxt_modulus);
    std::vector<ClueMeta> clue = AGOMR3_GenClue_Batch(allIDs, ptxt_modulus, 10, num_person, 10, PSparam.n, PSparam.ell, num_transaction, PSpk);

    cout << "Clue Generate Finished!\n" << endl;

    for (int g = 0; g < 1; g++) {
        num_set.push_back(100);
        num_set.push_back(101);

        id_num = num_set.size();
        range_r = id_num < 3 ? 39 : id_num * 14;
        generate_coeff_rangeCheck(range_r, ptxt_modulus);
        cout << "generat group : " << id_num << endl;
        
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

            // 평문 id set 사용
        std::vector<Vector> id_set(id_num);

        for (int i = 0; i < id_num; i++) {
            id_set[i] = allIDs[num_set[i]];
        }

        std::vector<std::vector<std::vector<uint64_t>>> signals_a(id_num,
            std::vector<std::vector<uint64_t>>(num_transaction, std::vector<uint64_t>(PSparam.n)));
        std::vector<std::vector<std::vector<uint64_t>>> signals_b(id_num,
            std::vector<std::vector<uint64_t>>(num_transaction, std::vector<uint64_t>(PSparam.ell)));

        for (int i = 0; i < id_num; i++) {
            retrieve(signals_a[i], signals_b[i], id_set[i], clue, PSpk);
        }

        clock_end = chrono::high_resolution_clock::now();
        cout << "\t Clue Retrieve time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

        clock_start = chrono::high_resolution_clock::now();
        detect(PV, signals_a, signals_b, context, PSsk_enc, swk);

        clock_end = chrono::high_resolution_clock::now();
        cout << "\t Detection Finished!" << endl;
        cout << "\t Detection time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

        /*
                    3-2. Compression
        */

        cout << "\t 3-2. Compression Started!\n" << endl;
        clock_start = chrono::high_resolution_clock::now();

        Ciphertext<DCRTPoly> digest;

        compress(digest, PV, payloads, context_comp, context_trace, keyPair_comp.publicKey, pubkey, keyPair_comp.secretKey);

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

        decode(decodedIdx, decodedPayload, digest, context_trace, secretKeys, keyPair_trace.secretKey);

        clock_end = chrono::high_resolution_clock::now();
        cout << "Client-side Decoding Finished!" << endl;
        cout << "Client-side Decoding Time: " << chrono::duration<double, milli>(clock_end - clock_start).count() << "ms\n" << endl;

        check(decodedIdx, clue, id_set);
    }

    return;
}