import numpy as np

def compute_coeffs(r, p):
    # 초기 다항식: 1 (상수항만)
    coeffs = [1]  # degree 0

    for i in range(r+1):
        new_factor = [(- (i**2)) % p, 0, 1]  # X^2 - i^2  → 계수 순서: 상수항, X^1항(0), X^2항
        # 다항식 곱셈: coeffs * new_factor
        new_coeffs = [0] * (len(coeffs) + 2)
        for deg1, c1 in enumerate(coeffs):
            for deg2, c2 in enumerate(new_factor):
                new_coeffs[deg1 + deg2] = (new_coeffs[deg1 + deg2] + c1 * c2) % p
        coeffs = new_coeffs

    # 이제 coeffs는 [a_0, a_1, ..., a_{2(r+1)}]
    # 필요에 따라 coeff 배열 포맷에 맞게 변환 가능
    return coeffs

r = 39
p = 786433
k = 8
arr = compute_coeffs(r,p)[::2]
reshaped = [arr[i*k:(i+1)*k] for i in range(k)]

for row in reshaped:
    print('{', end='')
    print(','.join(str(x) for x in row), end='')
    print('},')