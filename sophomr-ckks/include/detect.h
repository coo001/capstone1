#include "enco.h"

void affine(std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>>& output, 
            const  std::vector<std::vector<std::complex<double>>>& signals_a,
            const  std::vector<std::vector<std::complex<double>>>& signals_a_minus,
            const  std::vector<std::vector<std::complex<double>>>& signals_b,
            const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context, 
            const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& PSsk_enc)
{
    // Baby-step
    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> PSsk_bs(b_tilde1);
    PSsk_bs[0] = PSsk_enc;
    for (int b = 0; b < b_tilde1 - 1; b++) {
        PSsk_bs[b + 1] = context->EvalRotate(PSsk_bs[b], 1);
    }

    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sk_j;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sk_sq;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sk_j_plus_1;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sk_j_minus_1;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> term_pos;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> term_neg;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> term_j;
    lbcrypto::Plaintext exp_aj_plain;
    lbcrypto::Plaintext exp_minus_aj_plain;
    lbcrypto::Plaintext half_plain;
    lbcrypto::Plaintext one_plain;
    one_plain = context->MakeCKKSPackedPlaintext(std::vector<double>{1.0});
    half_plain = context->MakeCKKSPackedPlaintext(std::vector<double>{0.5});

    // Giant-step
    for (int i = 0; i < numctxt; i++) {
        int idx = i * degree;
        std::vector<std::complex<double>> vec_ints(degree);
        std::vector<std::complex<double>> vec_ints2(degree);

        for (int l = 0; l < PSparam.ell; l++) {
            for (int g_ = 0; g_ < g_tilde1; g_++) {
                int g = g_tilde1 - g_ - 1;
                lbcrypto::Ciphertext<lbcrypto::DCRTPoly> giant;

                for (int b = 0; b < b_tilde1; b++) {
                    int k = g * b_tilde1 + b;
                    if (k >= PSparam.n) { break; }

                    sk_j = PSsk_bs[b];  // 암호문: packed form에서 j번째 위치로 회전됨
                    sk_j_plus_1 = context->EvalAdd(sk_j, 1.0);
                    sk_j_minus_1 = context->EvalSub(sk_j, 1.0);

                    for (int m = 0; m < degree; m++) {
                        int idxc = (k + m) % PSparam.n;
                        if (idxc <= l) {
                            vec_ints[m] = signals_a[idx + m][l - idxc];
                            vec_ints2[m] = signals_a_minus[idx + m][l - idxc];
                        }
                        else {
                            vec_ints[m] = signals_a_minus[idx + m][PSparam.n + l - idxc];
                            vec_ints2[m] = signals_a[idx + m][PSparam.n + l - idxc];
                        }
                    }

                    int idx_rot = g * b_tilde1;
                    rotate(vec_ints.begin(), vec_ints.begin() + degree_half - idx_rot, vec_ints.begin() + degree_half);
                    rotate(vec_ints.begin() + degree_half, vec_ints.begin() + degree - idx_rot, vec_ints.begin() + degree);
                    rotate(vec_ints2.begin(), vec_ints2.begin() + degree_half - idx_rot, vec_ints2.begin() + degree_half);
                    rotate(vec_ints2.begin() + degree_half, vec_ints2.begin() + degree - idx_rot, vec_ints2.begin() + degree);

                    exp_aj_plain = context->MakeCKKSPackedPlaintext(vec_ints);
                    exp_minus_aj_plain = context->MakeCKKSPackedPlaintext(vec_ints2);

                    term_pos = context->EvalMult(sk_j, sk_j_plus_1);          // sk*(sk+1)
                    term_pos = context->EvalMult(term_pos, exp_aj_plain);
                    term_pos = context->EvalMult(term_pos, half_plain);

                    term_neg = context->EvalMult(sk_j, sk_j_minus_1);         // sk*(sk−1)
                    term_neg = context->EvalMult(term_neg, exp_minus_aj_plain);
                    term_neg = context->EvalMult(term_neg, half_plain);

                    term_j = context->EvalAdd(term_pos, term_neg);
                    sk_sq = context->EvalMult(sk_j, sk_j); // sk^2
                    sk_sq = context->EvalSub(one_plain, sk_sq); // 1 - sk^2

                    term_j = context->EvalAdd(term_j, sk_sq);  // 최종 term_j

                    if (b == 0) {
                        giant = term_j;
                    }
                    else {
                        giant = context->EvalMult(giant, term_j);
                        context->RescaleInPlace(giant);
                    }
                }

                if (g_ == 0) {
                    output[i][l] = giant;
                }
                else {
                    output[i][l] = context->EvalRotate(output[i][l], b_tilde1);
                    output[i][l] = context->EvalMult(output[i][l], giant);
                    context->RescaleInPlace(output[i][l]);
                }
            }
        }

        for (int l = 0; l < PSparam.ell; l++) {
            for (int m = 0; m < degree; m++) {
                vec_ints[m] = signals_b[idx + m][l];
            }
            output[i][l] = context->EvalMult(output[i][l], context->MakeCKKSPackedPlaintext(vec_ints));
        }
    }
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

void powers(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& inplace, const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context)
{
    inplace[1] = inplace[0];
    for (size_t i = 2; i < inplace.size(); i++) {
        if (i % 2 == 0) {
            inplace[i] = context->EvalSquare(inplace[i / 2]);
        }
        else {
            inplace[i] = context->EvalMult(inplace[(i - 1) / 2], inplace[(i + 1) / 2]);
        }
    }
}

std::vector<std::vector<double>> generate_coeff_rangecheck_cos_based(double theta, int degree) {
    std::vector<double> coeffs(degree + 1, 0.0);
    std::vector<std::vector<double>> coeff_mat;

    coeffs[0] = 0.5;
    int k = 10;
    if (degree >= 1)
        coeffs[1] = k / 4.0;
    if (degree >= 3)
        coeffs[3] = -std::pow(k, 3) / 48.0;
    if (degree >= 5)
        coeffs[5] = std::pow(k, 5) / 480.0;
    if (degree >= 7)
        coeffs[7] = -std::pow(k, 7) / 4032.0;
    if (degree >= 9)
        coeffs[9] = std::pow(k, 9) / 46080.0;

    int numbs = static_cast<int>(std::ceil(std::sqrt(degree + 1)));
    int numgs = static_cast<int>(std::ceil((degree + 1.0) / numbs));
    coeff_mat.resize(numgs, std::vector<double>(numbs, 0));

    for (int g = 0; g < numgs; ++g) {
        for (int b = 0; b < numbs; ++b) {
            int idx = g * numbs + b;
            if (idx <= degree) {
                coeff_mat[g][b] = coeffs[idx];
            }
        }
    }

    return coeff_mat;
}

void patersonStockmeyer(lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& inplace,
    const std::vector<std::vector<double>>& coeff,
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context)
{
    int numbs = coeff[0].size();
    int numgs = coeff.size();

    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> bs(numbs + 1);
    bs[0] = inplace;
    powers(bs, context);

    std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> gs(numgs);
    gs[0] = bs[numbs];
    powers(gs, context);

    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> hoist;

    for (int g = 0; g < numgs; g++) {
        lbcrypto::Ciphertext<lbcrypto::DCRTPoly> bsSum;
        for (int b = 1; b < numbs; b++) {
            std::vector<double> vec(degree, coeff[g][b]);
            auto ptxt = context->MakeCKKSPackedPlaintext(vec);

            auto temp = context->EvalMult(bs[b], ptxt);
            if (b == 1) {
                bsSum = temp;
            }
            else {
                context->EvalAddInPlace(bsSum, temp);
            }
        }

        std::vector<double> vec(degree, coeff[g][0]);
        auto ptxt = context->MakeCKKSPackedPlaintext(vec);
        context->EvalAddInPlace(bsSum, ptxt);

        if (g == 0) {
            inplace = bsSum;
        }
        else if (g == 1) {
            hoist = context->EvalMultNoRelin(bsSum, gs[g]);
        }
        else {
            bsSum = context->EvalMultNoRelin(bsSum, gs[g]);
            context->EvalAddInPlace(hoist, bsSum);
        }
    }

    context->RelinearizeInPlace(hoist);
    context->EvalAddInPlace(inplace, hoist);
}

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> InverseCKKS(
    const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& c,
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context,
    size_t r) 
{
    auto one_plain = context->MakeCKKSPackedPlaintext(std::vector<std::complex<double>>(1, 1.0));
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> one_minus_c = context->EvalSub(one_plain, c);

    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sum;
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> term;

    for (size_t k = 1; k <= r; ++k) {
        if (k == 1) {
            term = one_minus_c;
            sum = context->EvalAdd(term, one_plain);
            continue;
        }
        term = context->EvalMult(term, one_minus_c);
        term = context->Rescale(term);

        sum = context->EvalAdd(sum, term);
    }

    return sum;
}

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> EvalRealPart(
    const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& c,
    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context)
{
    using namespace lbcrypto;

    auto c_conj = InverseCKKS(c, context, 5);
    auto c_real = context->EvalAdd(c, c_conj);
    c_real = context->EvalMult(c_real, 0.5);

    return c_real;
}

void rangeCheck(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& output, 
                    std::vector<std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>>& input, 
                    const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context, 
                    const lbcrypto::EvalKey<lbcrypto::DCRTPoly>& swk) 
{ 
    // One-vector for Logical Negation
    std::vector<double> vec_one(degree, 1);
    auto ptxt_one = context->MakeCKKSPackedPlaintext(vec_one);


    double theta = std::cos(2.0 * M_PI * static_cast<double>(range) / static_cast<double>(ptxt_modulus));
    auto coeffs = generate_coeff_rangecheck_cos_based(theta, 10);

    for (int i = 0; i < numctxt; i++){
        for(int j = 0; j < PSparam.ell; j++){
            input[i][j] = EvalRealPart(input[i][j], context);
            patersonStockmeyer(input[i][j], coeffs, context);
        }
    }
    
    for (int i = 0; i < numctxt; i++){
        output[i] = product(input[i], context);
        context->GetScheme()->KeySwitchInPlace(output[i], swk);
    }
} 

void detect(std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>>& output, 
                const  std::vector<std::vector<std::complex<double>>>& signals_a,
                const  std::vector<std::vector<std::complex<double>>>& signals_a_minus,
                const  std::vector<std::vector<std::complex<double>>>& signals_b,
                const lbcrypto::CryptoContext<lbcrypto::DCRTPoly>& context, 
                const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& PSsk_enc, 
                const lbcrypto::EvalKey<lbcrypto::DCRTPoly>& swk)
{     
    using namespace std;
    using namespace lbcrypto;
    
    chrono::high_resolution_clock::time_point clock_start, clock_end;

    vector<vector<Ciphertext<DCRTPoly>>> PV_ell(numctxt, vector<Ciphertext<DCRTPoly>>(PSparam.ell));

    cout << "\t\t Affine Transform Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();
    
    affine(PV_ell, signals_a, signals_a_minus, signals_b, context, PSsk_enc);
    //uint32_t numTowers1 = PV_ell[0][0]->GetElements()[0].GetNumOfElements();
    //cout << "num towers : " << numTowers1 << endl;
    clock_end = chrono::high_resolution_clock::now();
    cout << "\t\t Affine Transform Finished!" << endl;
    cout << "\t\t Affine Transform time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;

    cout << "\t\t RangeCheck Started!\n" << endl;
    clock_start = chrono::high_resolution_clock::now();
   
    rangeCheck(output, PV_ell, context, swk);

    clock_end = chrono::high_resolution_clock::now();
    cout << "\t\t RangeCheck Finished!" << endl;
    cout << "\t\t RangeCheck time: " << chrono::duration<double>(clock_end - clock_start).count() << "sec\n" << endl;
}