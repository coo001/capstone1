#include "openfhe.h"

void affine(std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>>& output, 
            const std::vector<std::vector<uint64_t>>& signals_a, 
            const std::vector<std::vector<uint64_t>>& signals_b, 
            const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context, 
            const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& PSsk_enc)
{   
    using namespace std;
    using namespace lbcrypto;
    chrono::high_resolution_clock::time_point clock_start, clock_end;
    // Baby-step
    clock_start = chrono::high_resolution_clock::now();
    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> PSsk_bs(b_tilde1);
    PSsk_bs[0] = PSsk_enc;

    for (int b = 0; b < b_tilde1 - 1; b++) {
        PSsk_bs[b+1] = context->EvalRotate(PSsk_bs[b], 1);
    }
    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Baby step time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    clock_start = chrono::high_resolution_clock::now();
    // Giant-step
    for (int i = 0; i < numctxt; i++) {
        int idx = i * degree;
        std::vector<int64_t> vec_ints(degree);

        for (int l = 0; l < PSparam.ell; l++) {
            for (int g_ = 0; g_ < g_tilde1 ; g_++) {
                int g = g_tilde1 - g_ - 1;
                lbcrypto::Ciphertext<lbcrypto::DCRTPoly> giant;

                for (int b = 0; b < b_tilde1; b++) {
                    int k = g * b_tilde1 + b;      
                    if ( k >= PSparam.n ) {break;}

                    for(int m = 0; m < degree; m++){
                        int idxc = (k + m) % PSparam.n;
                        if (idxc <= l) {
                            vec_ints[m] = signals_a[idx + m][l - idxc];
                        } else {
                            vec_ints[m] = PSparam.q - signals_a[idx + m][PSparam.n + l - idxc];
                            if (vec_ints[m] == PSparam.q) { vec_ints[m] = 0; }
                        }
                    }

                    int idx_rot = g * b_tilde1; 
                    rotate(vec_ints.begin(), vec_ints.begin() + degree_half - idx_rot, vec_ints.begin() + degree_half);
                    rotate(vec_ints.begin() + degree_half, vec_ints.begin() + degree - idx_rot, vec_ints.begin() + degree);

                    if ( b == 0 ) {
                        giant = context->EvalMult(PSsk_bs[0], context->MakePackedPlaintext(vec_ints));
                    } else {
                        auto temp = context->EvalMult(PSsk_bs[b], context->MakePackedPlaintext(vec_ints));
                        context->EvalAddInPlace(giant, temp);
                    }
                }

                if (g_ == 0) {
                    output[i][l] = giant;
                } else {
                    output[i][l] = context->EvalRotate(output[i][l], b_tilde1);
                    context->EvalAddInPlace(output[i][l], giant);
                }
            }
        }

        for (int l = 0; l < PSparam.ell; l++){
            for(int m = 0; m < degree; m++){
                vec_ints[m] = signals_b[idx + m][l];
            }
            output[i][l] = context->EvalSub(output[i][l], context->MakePackedPlaintext(vec_ints)); 
        }
    }
    clock_end = chrono::high_resolution_clock::now();
    cout << "\t Giant step time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;
}

void powers(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& inplace, const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context) 
{
    inplace[1] = inplace[0];
    for (size_t i = 2; i < inplace.size(); i++) {
        if(i % 2 == 0){
            inplace[i] = context->EvalSquare(inplace[i/2]);
        } else {
            inplace[i] = context->EvalMult(inplace[(i-1)/2], inplace[(i+1)/2]);
        }
    }
}

void patersonStockmeyer(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& inplace, 
                            int r,
                            //const std::vector<std::vector<uint32_t>>& coeff, 
                            const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context,
                            const lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& secretKey)
{
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> x;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> x2;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> x4;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> t;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> k;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> k2;
    lbcrypto::Plaintext ptxt;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> temp;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> temp2;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> tum;

    x = inplace;
    x2 = context->EvalSquare(inplace);
    x4 = context->EvalSquare(x2);

    std::vector<int64_t> vec(degree, (r + 1)*(r + 1));
    ptxt = context->MakePackedPlaintext(vec);
    temp = context->EvalMult(x2, ptxt);
    t = context->EvalSub(x4, temp);

    vec.assign(degree, 2 * r - 2);
    ptxt = context->MakePackedPlaintext(vec);
    k = context->EvalMult(x2, ptxt);
    k2 = context->EvalMult(x2, ptxt);
    vec.assign(degree, 2 * r + 2);
    ptxt = context->MakePackedPlaintext(vec);
    temp = context->EvalMult(x, ptxt);
    context->EvalAddInPlace(k, temp);
    context->EvalSubInPlace(k2, temp);
    vec.assign(degree, 2 * r - 2);
    ptxt = context->MakePackedPlaintext(vec);
    context->EvalAddInPlace(k, ptxt);
    context->EvalAddInPlace(k2, ptxt);

    vec.assign(degree, 4);
    ptxt = context->MakePackedPlaintext(vec);
    tum = context->EvalMult(x2, ptxt);
    vec.assign(degree, 2);
    ptxt = context->MakePackedPlaintext(vec);
    context->EvalAddInPlace(tum, ptxt);

    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> a;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> b;
    a = t; b = t;
    for (int i = 1; i <= r / 4; ++i) {
        if (i == 1) { //a1, b1 설정
            vec.assign(degree, 3 * r - 2);
            ptxt = context->MakePackedPlaintext(vec);
            temp = context->EvalMult(x2, ptxt);
            vec.assign(degree, r*r - r - 2);
            ptxt = context->MakePackedPlaintext(vec);
            temp2 = context->EvalMult(x, ptxt);
            context->EvalSubInPlace(temp, temp2);
            context->EvalAddInPlace(a, temp);
            
            temp = context->EvalAdd(temp2, temp2);
            b = context->EvalAdd(a, temp);

            vec.assign(degree, 2 * r * (r - 1));
            ptxt = context->MakePackedPlaintext(vec);
            context->EvalAddInPlace(a, ptxt);
            context->EvalAddInPlace(b, ptxt);
        }
        else { // a_i, b_i 설정
            context->EvalSubInPlace(k, tum);
            context->EvalAddInPlace(a, k);
            context->EvalSubInPlace(k2, tum);
            context->EvalAddInPlace(b, k2);
        }

        inplace = context->EvalMult(inplace, a);
        inplace = context->EvalMult(inplace, b);
    }
}

void FLT(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& inplace, const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context) 
{
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> curr = inplace;

    auto pow = ptxt_modulus - 1;
    bool flag_first = true;

    while(pow > 0){
        if(pow % 2 == 0){
            pow /= 2;
            context->EvalSquareInPlace(curr);
        } else {
            pow -= 1;
            if (flag_first) {
                inplace = curr;
                flag_first = false;
                continue;
            }
            inplace = context->EvalMult(inplace, curr);
        }
    }
    //std::cout << inplace << std::endl;
}

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> product(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& input, 
                                                    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context)
{
    while(input.size() != 1){
        for(size_t i = 0; i < input.size()/2; i++){
            input[i] = context->EvalMult(input[i], input[input.size()/2+i]);
        }
        if(input.size() % 2 == 0)
            input.resize(input.size()/2);
        else{
            input[input.size()/2] = input[input.size()-1];
            input.resize(input.size()/2+1);
        }
    }

    return input[0];
}

void rangeCheck(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& output, 
                    std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>>& input, 
                    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context, 
                    const lbcrypto::EvalKey<lbcrypto::DCRTPoly>& swk,
                    const lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& sk)
{ 
    using namespace std;
    using namespace lbcrypto;
    chrono::duration<double> clockg = chrono::duration<double>::zero();
    chrono::duration<double> clockh = chrono::duration<double>::zero();
    chrono::high_resolution_clock::time_point clock_start, clock_end;

    //std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>> ig = input;

    // One-vector for Logical Negation
    std::vector<int64_t> vec_one(degree, 1);
    auto ptxt_one = context->MakePackedPlaintext(vec_one);

    for (int i = 0; i < numctxt; i++){
        for(int j = 0; j < PSparam.ell; j++){
            
            context->EvalSquareInPlace(input[i][j]);

            // Zeroize
            clock_start = chrono::high_resolution_clock::now();
            patersonStockmeyer(input[i][j], range_r, context, sk);
            //patersonStockmeyer(input[i][j], coeff_rangeCheck, context, sk);

            clock_end = chrono::high_resolution_clock::now();
            clockg += chrono::duration<double>(clock_end - clock_start);

            // Fermat's Little Theorem
            clock_start = chrono::high_resolution_clock::now();
            FLT(input[i][j], context);
            clock_end = chrono::high_resolution_clock::now();
            clockh += chrono::duration<double>(clock_end - clock_start);
        
            // Logical Negation
            input[i][j] = context->EvalSub(ptxt_one, input[i][j]);
        }
    }

    cout << "\t Zeroize Time: " << clockg.count() << "sec\n" << endl;
    cout << "\t FLT Time: " << clockh.count() << "sec\n" << endl;
    
    for (int i = 0; i < numctxt; i++){
        output[i] = product(input[i], context);
        output[i] = context->Compress(output[i], 2);
        context->GetScheme()->KeySwitchInPlace(output[i], swk);
    }
} 

std::vector<std::vector<uint64_t>> AddVectors(
    const std::vector<std::vector<uint64_t>>& v1,
    const std::vector<std::vector<uint64_t>>& v2)
{
    std::vector<std::vector<uint64_t>> result(v1.size());

    for (size_t i = 0; i < v1.size(); ++i) {
        result[i].resize(v1[i].size());
        for (size_t j = 0; j < v1[i].size(); ++j) {
            result[i][j] = (v1[i][j] + v2[i][j]) % ptxt_modulus;
        }
    }

    return result;
}

void detect(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& output, 
                const std::vector<std::vector<uint64_t>>& signals_a, 
                const std::vector<std::vector<uint64_t>>& signals_b, 
                const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context, 
                const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& PSsk_enc, 
                const lbcrypto::EvalKey<lbcrypto::DCRTPoly>& swk,
                const lbcrypto::PrivateKey<lbcrypto::DCRTPoly>& sk)
{     
    using namespace std;
    using namespace lbcrypto;
    
    chrono::high_resolution_clock::time_point clock_start, clock_end;

    vector<vector<Ciphertext<DCRTPoly>>> PV_ell(numctxt, vector<Ciphertext<DCRTPoly>>(PSparam.ell));
    
    cout << "\t\t Affine Transform Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();
    
    affine(PV_ell, signals_a, signals_b, context, PSsk_enc);
    
    clock_end = chrono::high_resolution_clock::now();
    cout << "\t\t Affine Transform Finished!" << endl;
    cout << "\t\t Affine Transform time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;


    cout << "\t\t RangeCheck Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();
   
    rangeCheck(output, PV_ell, context, swk, sk);

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t\t RangeCheck Finished!" << endl;
    cout << "\t\t RangeCheck time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;
}