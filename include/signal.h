
struct param{
    int n; int q; int ell; int h; double sigma;
    param(int n, int q, int ell, int h, double sigma)
    : n(n), q(q), ell(ell), h(h), sigma(sigma)
    {}
};

typedef NativeVector PSsk;

struct Signal{
    NativeVector a;
    NativeVector b;
};

typedef Signal PSpk;

PSsk PSskGen(const param& param) 
{
    lbcrypto::TernaryUniformGeneratorImpl<NativeVector> tug;
    PSsk sk = tug.GenerateVector(param.n, param.q, param.h);
    return sk;
}

std::vector<std::vector<int64_t>> expand_nomod(const NativeVector& a, int n) {
    std::vector<std::vector<int64_t>> res(n);

    for (int cnt = 0; cnt < n; cnt++) {
        res[cnt].resize(n);
        int ind = 0;

        // Forward part: a[cnt], a[cnt-1], ..., a[0]
        for (int i = cnt; i >= 0 && ind < n; i--) {
            res[cnt][ind++] = static_cast<int64_t>(a[i].ConvertToInt());
        }

        // Backward part: -a[n-1], -a[n-2], ..., -a[cnt+1]
        for (int i = n - 1; i > cnt && ind < n; i--) {
            res[cnt][ind++] = -static_cast<int64_t>(a[i].ConvertToInt());
        }
    }

    return res;
}

NativeVector ringMult_nomod(NativeVector a, NativeVector b, int n, int q, int ell) 
{
    NativeVector res = NativeVector(ell);

    std::vector<std::vector<int64_t>> a_expanded = expand_nomod(a, n);
    for (int i = 0; i < ell; i++) {
        int64_t temp = 0;
        for (int j = 0; j < n; j++) {
            int64_t aj = a_expanded[i][j];
            int64_t bj = static_cast<int64_t>(b[j].ConvertToInt());
            temp += aj * bj;
        }
        res[i] = lbcrypto::NativeInteger(temp);   // 실제론 오버플로 주의해야 함
    }

    return res;
}



PSpk PSpkGen(const param& param, const PSsk& sk) 
{   
    PSpk pk;

    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> dug;
    dug.SetModulus(param.q);

    pk.a = dug.GenerateVector(param.n);
    pk.b = ringMult_nomod(pk.a, sk, param.n, param.q, param.n);

    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> dgg(param.sigma);
    for (int i = 0; i < param.n; i++) {
        pk.b[i].ModAddFastEq(dgg.GenerateInteger(param.q), param.q);
    }

    return pk;
}

void PSsignal(Signal& sig, const PSpk& pk, const param& param) 
{
    lbcrypto::TernaryUniformGeneratorImpl<NativeVector> tug;
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> dgg(param.sigma);

    NativeVector r = tug.GenerateVector(param.n, param.q, param.h);
    sig.a = ringMult_nomod(pk.a, r, param.n, param.q, param.n);
    sig.b = ringMult_nomod(pk.b, r, param.n, param.q, param.ell);

    for(int i = 0; i < param.n; i++){
        sig.a[i] += dgg.GenerateInteger(param.q);
        //sig.a[i].ModAddFastEq(dgg.GenerateInteger(param.q), param.q);
    }
    for(int i = 0; i < param.ell; i++){
        sig.b[i] += dgg.GenerateInteger(param.q);
        //sig.b[i].ModAddFastEq(dgg.GenerateInteger(param.q), param.q);
    }
}