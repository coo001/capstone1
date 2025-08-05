#ifndef ENCO_H
#define ENCO_H

#include <vector>
#include <complex>
#include <cmath>
#include <cstdint>

constexpr double PI = 3.14159265358979323846;

inline double symmetric_embed(uint64_t x, uint64_t p) {
    return (x < p / 2) ? (double)x : (double)x - (double)p;
}

// 1. 정수 m을 받으면 exp(-2*pi*i*m/n)을 반환
std::complex<double> encodeInt(int m, int n) {
    double angle = 2.0 * PI * static_cast<double>(m) / static_cast<double>(n);
    return std::complex<double>(cos(angle), sin(angle));
}

// 2. 복소수 하나를 받아 가장 가까운 m을 찾는 decode
int decodeComplex(const std::complex<double>& c, int n) {
    double angle = std::arg(c); // [-pi, pi]
    if (angle < 0) angle += 2 * PI; // [0, 2pi)
    int m = static_cast<int>(std::round(angle * n / (2.0 * PI))) % n;
    return m;
}

// 3. 정수 벡터에 대해 encode
std::vector<std::complex<double>> encodeVector(const std::vector<double>& vec, int n) {
    std::vector<std::complex<double>> result;
    result.reserve(vec.size());
    for (int m : vec) {
        result.push_back(encodeInt(m, n));
    }
    return result;
}

// 4. 복소수 벡터를 decode
std::vector<int> decodeVector(const std::vector<std::complex<double>>& vec, int n) {
    std::vector<int> result;
    result.reserve(vec.size());
    for (const auto& c : vec) {
        result.push_back(decodeComplex(c, n));
    }
    return result;
}

#endif