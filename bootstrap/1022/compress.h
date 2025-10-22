
void computeVandermonde(std::vector<std::vector<double>>& vandermonde) // uint64_t -> double
{
    for (size_t j = 0; j < vandermonde[0].size(); j++) {
        vandermonde[0][j] = static_cast<double>(j + 1);
    }
    for (size_t i = 1; i < vandermonde.size(); i++) {
        for (size_t j = 0; j < vandermonde[0].size(); j++) {
            // `% ptxt_modulus` ����
            vandermonde[i][j] = vandermonde[i - 1][j] * static_cast<double>(j + 1);
        }
    }
}

void ringSwitch(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& digest, 
                const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context_comp,
                const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context_trace,
                const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& publicKey_comp, 
                const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& publicKey_trace)
{
    using namespace std;
    using namespace lbcrypto;

    cout << "hi" << endl;

    vector<double> vec_mask1(degree,0);
    vector<double> vec_mask2(degree,0);
    for (int i = 0; i < degree_trace_half; i++) {
        vec_mask1[i] = dim_trace;
        vec_mask2[i + degree_half] = dim_trace;
    }
    auto ptxt_mask1 = context_comp->MakeCKKSPackedPlaintext(vec_mask1);
    auto ptxt_mask2 = context_comp->MakeCKKSPackedPlaintext(vec_mask2);
    //여기도 레벨조정해야해서함
    //size_t level = digest->GetElements()[0].GetNumOfElements();
    auto enc_ptxt1 = context_comp->Encrypt(publicKey_comp, ptxt_mask1);
    auto enc_ptxt2 = context_comp->Encrypt(publicKey_comp, ptxt_mask2);
    /*
    std::vector<DCRTPoly>& elements = enc_ptxt1->GetElements();
    for (auto& element : elements) {
        element.DropLastElements(MultiplicativeDepth-level+1);
    }
    std::vector<DCRTPoly>& elements2 = enc_ptxt2->GetElements();
    for (auto& element : elements2) {
        element.DropLastElements(MultiplicativeDepth-level+1);
    }
    */

    cout << "by" << endl;
    Ciphertext<DCRTPoly> temp;
    temp = context_comp->EvalRotate(digest, degree_trace_half);
    temp = context_comp->EvalMult(temp, enc_ptxt2);
    std::cout << 10 << std::endl;
    digest = context_comp->EvalMult(digest, enc_ptxt1);
    context_comp->EvalAddInPlace(digest, temp);

    //digest = context_comp->Compress(digest);

    vector<double> vec_null(degree_trace);
    auto ptxt_null = context_trace->MakeCKKSPackedPlaintext(vec_null);
    auto ctxt_trace = context_trace->Encrypt(publicKey_trace, ptxt_null);
    //ctxt_trace = context_trace->Compress(ctxt_trace);

    auto poly_comp = digest->GetElements();
    auto poly_trace = ctxt_trace->GetElements();
    for (int i = 0; i < 2; i++) {
        poly_comp[i].SwitchFormat();
        poly_trace[i].SwitchFormat();
    }
    for (int i = 0; i < 2; i++) {
        NativePoly poly_comp_ = poly_comp[i].GetElementAtIndex(0);
        NativePoly poly_trace_ = poly_trace[i].GetElementAtIndex(0);
        for (size_t k = 0; k < poly_trace_.GetLength(); k++) {
            poly_trace_[k] = poly_comp_[dim_trace * k].Mod(poly_trace_.GetModulus());
        }
        poly_trace[i].SetElementAtIndex(0, poly_trace_);
    }
    for (int i = 0; i < 2; i++) {
        poly_trace[i].SwitchFormat();
    }
    ctxt_trace->SetElements(poly_trace);

    digest = ctxt_trace;
}

void compress(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& digest, 
                std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& PV, 
                const std::vector<std::vector<uint32_t>>& payloads, 
                const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context_comp,
                const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context_trace, 
                const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& publicKey_comp, 
                const lbcrypto::PublicKey<lbcrypto::DCRTPoly>& publicKey_trace,
                const lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& privateKey_comp)
{
    using namespace std;
    using namespace lbcrypto;

    auto Dump = [&](auto ct, const char* tag) {
        std::cerr << tag << " | level=" << ct->GetLevel()
            << " | noise=" << ct->GetNoiseScaleDeg()
            << " | scale=" << ct->GetScalingFactor() << "\n";
        };

    for (int i = 0; i < numctxt; i++) {
        Dump(PV[i], "after detect+rangeCheck"); // ★ 여기서 level이 0이면 Compression 못 들어감
    }


    vector<double> vec_null(degree);
    auto ptxt_null = context_comp->MakeCKKSPackedPlaintext(vec_null);
    for (int i = 0; i < numctxt; i++) {
        auto ctxt_temp = context_comp->Encrypt(publicKey_comp, ptxt_null);
        ctxt_temp->SetElements(PV[i]->GetElements());
        PV[i] = ctxt_temp;
    }
    /*for (int i = 0; i < numctxt; i++) {
        // rangeCheck 출력은 0/1 근방이라 부트 입력 안정
        PV[i] = context_comp->EvalBootstrap(PV[i]);
        // 디버그(선택): 부트 직후 레벨 확인
        // DumpCt(PV[i], "after bootstrap@compress");
    }*/
    
    // Precompute Vandermonde
    vector<vector<double>> vandermonde(num_pertinent, vector<double>(payloads.size())); 
    computeVandermonde(vandermonde);
    
    // Baby-step
    vector<vector<Ciphertext<DCRTPoly>>> rotated_PV(numctxt, vector<Ciphertext<DCRTPoly>>(b_tilde2));
    for (int i = 0; i < numctxt; i++) {
        rotated_PV[i][0] = PV[i];
        for (int b = 0; b < b_tilde2-1; b++) {
            rotated_PV[i][b+1] = context_comp->EvalRotate(rotated_PV[i][b], 1);
        }
    }

    // Giant-step
    vector<Ciphertext<DCRTPoly>> giant(numctxt);
    vector<double> ptxt_vec(degree,0.0);
    std::cout << 0 << std::endl;

    for (int g_ = 0; g_ < g_tilde2 ; g_++) {
        int g = g_tilde2 - g_ - 1;

        for (int i = 0; i < numctxt; i++) {    
            int counter = i * degree;
            for (int b = 0; b < b_tilde2; b++) {                
                if ( g * b_tilde2 + b >= numrow_po2 ) {break;}

                int k2 = degree_half;
                int idxr = (numrow_po2 - g * b_tilde2) % numrow_po2; 
                int idxc1 = counter + b;
                int idxc2 = idxc1 + degree_half;

                for (int k1 = 0; k1 < degree_half; k1++) {
            
                    if (idxr >= numrow) {
                        ptxt_vec[k1] = 0;
                        ptxt_vec[k2] = 0;
                    } else if (idxr < num_pertinent) {
                        ptxt_vec[k1] = vandermonde[idxr][idxc1]; 
                        ptxt_vec[k2] = vandermonde[idxr][idxc2];
                    } else {
                        int idxr_q = idxr / num_pertinent - 1;
                        int idxr_r = idxr % num_pertinent;
                        ptxt_vec[k1] = (vandermonde[idxr_r][idxc1] * payloads[idxc1][idxr_q]);
                        ptxt_vec[k2] = (vandermonde[idxr_r][idxc2] * payloads[idxc2][idxr_q]);
                    }

                    k2++; idxr++; idxc1++; idxc2++;
                    if (idxr == numrow_po2) { idxr = 0; }
                    if (idxc2 == counter + degree) { 
                        idxc1 -= degree_half; 
                        idxc2 -= degree_half; 
                    }
                }

                //std::cout << "i : " << i << ", b : " << b << "\n";
                //auto& ctxt_to_mult = rotated_PV[i][b];
                //size_t level = ctxt_to_mult->GetElements()[0].GetNumOfElements();
                //std::cout << "level: " << level << std::endl;

                //double ctxt_scale = ctxt_to_mult->GetScalingFactor();

                /*20251009
                auto ptxt = context_comp->MakeCKKSPackedPlaintext(ptxt_vec);

                //ptxt->SetScalingFactor(ctxt_scale);
                //giant[i] = rotated_PV[i][b];
                auto enc_ptxt = context_comp->Encrypt(publicKey_comp, ptxt);
                std::vector<DCRTPoly>& elements = enc_ptxt->GetElements();
                for (auto& element : elements) {
                    element.DropLastElements(MultiplicativeDepth-level+2);
                }
                //enc_ptxt->SetLevel(enc_ptxt->GetLevel()+ level+1);
                //std::cout << "level: " << enc_ptxt->GetElements()[0].GetNumOfElements() << std::endl;
                //std::cout << "ptxt �� �� ��" << std::endl;
                if ( b == 0 ) {
                    giant[i] = context_comp->EvalMult(rotated_PV[i][0], enc_ptxt);
                } else {
                    auto temp = context_comp->EvalMult(rotated_PV[i][b], enc_ptxt);
                    context_comp->EvalAddInPlace(giant[i], temp);
                }
                */
                // 평문을 ct의 레벨에 맞춰 "깊이(depth)"만 지정해 생성
                // ct 레벨을 안전하게 가져와 depth 계산 (언더플로우 방지)
                // ct의 현재 레벨/스케일 확인
                auto& ct = rotated_PV[i][b];
                uint32_t lvl = ct->GetLevel();                    // 예: 2
                auto     sc = ct->GetScalingFactor();

                // ✅ 평문을 '같은 레벨'로 생성 (params=nullptr, slots=degree 주의)
                auto ptxt = context_comp->MakeCKKSPackedPlaintext(
                    ptxt_vec,
                    /*scaleDeg=*/1,                // 보통 1이면 충분. 필요시 조정
                    /*level   =*/ lvl,             // ★ ct와 동일 레벨
                    /*params  =*/ nullptr,
                    /*slots   =*/ static_cast<usint>(degree)
                );

                // (권장) 스케일을 ct에 맞춤 — 곱 뒤 정합 안정화
                ptxt->SetScalingFactor(sc);

                // 이제 ct × pt
                auto tmp = context_comp->EvalMult(ct, ptxt);
                if (b == 0) giant[i] = tmp;
                else        context_comp->EvalAddInPlace(giant[i], tmp);

             
            }
        }

        Ciphertext<DCRTPoly> sum = giant[0];
        for (int i = 1; i < numctxt; i++) {
            context_comp->EvalAddInPlace(sum, giant[i]);
        }

        if (g_ == 0) {
            digest = sum;
        } else {
            digest = context_comp->EvalRotate(digest, b_tilde2);
            context_comp->EvalAddInPlace(digest, sum);
        }
    }
    std::cout << 1 << std::endl;

    // === before EvalBootstrap ===
    Dump(digest, "pre-boot");

    if (digest->GetLevel() == 0) {
        throw std::runtime_error("EvalBootstrap input level=0. Move bootstrap earlier or adjust params.");
    }

    digest = context_comp->EvalBootstrap(digest);
    
    Ciphertext<DCRTPoly> temp;
    for (int j = 1; j < degree_half / numrow_po2; j*=2) {
        temp = context_comp->EvalRotate(digest, numrow_po2 * j);
        context_comp->EvalAddInPlace(digest, temp);
    }
    temp = context_comp->EvalRotate(digest, degree_half);
    context_comp->EvalAddInPlace(digest, temp);

    std::cout << 2 << std::endl;
    // Ring-Switching
    ringSwitch(digest, context_comp, context_trace, publicKey_comp, publicKey_trace);
}


auto ctxt_temp = cc->Encrypt(pk, ptxt0);
ctxt_temp->SetElements(PV[i]->GetElements());
PV[i] = ctxt_temp;

std::vector<DCRTPoly>& elements = enc_ptxt->GetElements();
for (auto& element : elements) {
    element.DropLastElements(MultiplicativeDepth - level + 2);
}