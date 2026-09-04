#ifndef INF_HPP_
#define INF_HPP_

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace inf {

class integer;
struct primality;
integer prime(int k, primality check);

template<typename int_u>
struct qr
{
    int_u q, r;
};

/**
 * Marsaglia XorShift64* — used by integer::rand and prime generation.
 */
class XorShift64
{
public:
    explicit XorShift64(uint64_t seed)
        : state(seed ? seed : 0x9E3779B97F4A7C15ULL)
    {}
    uint64_t next()
    {
        uint64_t x = state;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        state = x;
        return x * 0x2545F4914F6CDD1DULL;
    }
private:
    uint64_t state;
};

inline int clz32(uint32_t x)
{
    if (x == 0)
        return 32;
    int n = 0;
    if (x <= 0x0000FFFFu) { n += 16; x <<= 16; }
    if (x <= 0x00FFFFFFu) { n += 8;  x <<= 8;  }
    if (x <= 0x0FFFFFFFu) { n += 4;  x <<= 4;  }
    if (x <= 0x3FFFFFFFu) { n += 2;  x <<= 2;  }
    if (x <= 0x7FFFFFFFu) { n += 1; }
    return n;
}

inline int clz64(uint64_t x)
{
    if (x == 0)
        return 64;
    if (x > 0xFFFFFFFFull)
        return clz32((uint32_t)(x >> 32));
    return 32 + clz32((uint32_t)x);
}

inline void mul_u64(uint64_t a, uint64_t b, uint64_t& lo, uint64_t& hi)
{
    uint64_t a_lo = (uint32_t)a, a_hi = a >> 32;
    uint64_t b_lo = (uint32_t)b, b_hi = b >> 32;
    uint64_t p0 = a_lo * b_lo;
    uint64_t p1 = a_lo * b_hi;
    uint64_t p2 = a_hi * b_lo;
    uint64_t p3 = a_hi * b_hi;
    uint64_t mid = (p0 >> 32) + (uint32_t)p1 + (uint32_t)p2;
    lo = (p0 & 0xFFFFFFFFull) | (mid << 32);
    hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
}

inline uint32_t ntt_mod_pow(uint32_t a, uint32_t e, uint32_t m)
{
    uint64_t r = 1, b = a;
    while (e)
    {
        if (e & 1)
            r = r * b % m;
        b = b * b % m;
        e >>= 1;
    }
    return (uint32_t)r;
}

inline uint32_t ntt_mod_inv(uint32_t a, uint32_t m)
{
    return ntt_mod_pow(a, m - 2, m);
}

inline void ntt_bit_reverse(std::vector<uint32_t>& a)
{
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i)
    {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }
}

inline void ntt_transform(std::vector<uint32_t>& a, uint32_t mod, uint32_t root, bool invert)
{
    const size_t n = a.size();
    ntt_bit_reverse(a);
    for (size_t len = 2; len <= n; len <<= 1)
    {
        uint32_t wlen = ntt_mod_pow(root, (mod - 1) / (uint32_t)len, mod);
        if (invert)
            wlen = ntt_mod_inv(wlen, mod);
        for (size_t i = 0; i < n; i += len)
        {
            uint32_t w = 1;
            for (size_t j = 0; j < len / 2; ++j)
            {
                uint32_t u = a[i + j];
                uint32_t v = (uint32_t)((uint64_t)a[i + j + len / 2] * w % mod);
                uint32_t x = u + v;
                if (x >= mod)
                    x -= mod;
                uint32_t y = u + mod - v;
                if (y >= mod)
                    y -= mod;
                a[i + j] = x;
                a[i + j + len / 2] = y;
                w = (uint32_t)((uint64_t)w * wlen % mod);
            }
        }
    }
    if (invert)
    {
        uint32_t ninv = ntt_mod_inv((uint32_t)n, mod);
        for (size_t i = 0; i < n; ++i)
            a[i] = (uint32_t)((uint64_t)a[i] * ninv % mod);
    }
}

/**
 * big integer — little-endian base-2^32 limbs, sign +1 / -1
 */
class integer
{
private:
    std::vector<uint32_t> limbs;
    int_fast32_t sign;

    void trim()
    {
        while (!limbs.empty() && limbs.back() == 0)
            limbs.pop_back();
        if (limbs.empty())
            sign = 1;
    }

    static int cmp_abs(const integer& a, const integer& b)
    {
        if (a.limbs.size() != b.limbs.size())
            return a.limbs.size() < b.limbs.size() ? -1 : 1;
        for (int i = (int)a.limbs.size() - 1; i >= 0; --i)
        {
            if (a.limbs[i] != b.limbs[i])
                return a.limbs[i] < b.limbs[i] ? -1 : 1;
        }
        return 0;
    }

    bool even() const
    {
        return limbs.empty() || (limbs[0] & 1u) == 0;
    }

    static integer from_limbs(std::vector<uint32_t> v, int_fast32_t s = 1)
    {
        integer r;
        r.limbs = std::move(v);
        r.sign = s;
        r.trim();
        return r;
    }

    friend integer schoolbook(const integer& a, const integer& b)
    {
        if (a.limbs.empty() || b.limbs.empty())
            return integer();
        integer r;
        r.sign = a.sign * b.sign;
        r.limbs.assign(a.limbs.size() + b.limbs.size(), 0);
        for (size_t i = 0; i < a.limbs.size(); ++i)
        {
            uint64_t carry = 0;
            for (size_t j = 0; j < b.limbs.size(); ++j)
            {
                uint64_t cur = (uint64_t)r.limbs[i + j]
                    + (uint64_t)a.limbs[i] * b.limbs[j]
                    + carry;
                r.limbs[i + j] = (uint32_t)cur;
                carry = cur >> 32;
            }
            r.limbs[i + b.limbs.size()] = (uint32_t)carry;
        }
        r.trim();
        return r;
    }

    integer low_limbs(size_t m) const
    {
        integer r;
        r.sign = 1;
        if (m == 0 || limbs.empty())
            return r;
        size_t n = std::min(m, limbs.size());
        r.limbs.assign(limbs.begin(), limbs.begin() + (std::ptrdiff_t)n);
        r.trim();
        return r;
    }

    integer high_limbs(size_t m) const
    {
        integer r;
        r.sign = 1;
        if (limbs.size() <= m)
            return r;
        r.limbs.assign(limbs.begin() + (std::ptrdiff_t)m, limbs.end());
        r.trim();
        return r;
    }

    static integer shl_limbs(const integer& a, size_t m)
    {
        if (a.limbs.empty() || m == 0)
            return a;
        integer r;
        r.sign = a.sign;
        r.limbs.assign(m, 0);
        r.limbs.insert(r.limbs.end(), a.limbs.begin(), a.limbs.end());
        return r;
    }

    void read(const std::string& str)
    {
        sign = 1;
        limbs.clear();
        size_t pos = 0;
        while (pos < str.size() && (str[pos] == '-' || str[pos] == '+'))
        {
            if (str[pos] == '-')
                sign = -sign;
            ++pos;
        }
        while (pos < str.size() && str[pos] == '0')
            ++pos;
        if (pos >= str.size())
        {
            sign = 1;
            return;
        }
        const uint32_t DEC = 1000000000u;
        size_t first = (str.size() - pos) % 9;
        if (first == 0)
            first = 9;
        auto take = [&](size_t from, size_t len) -> uint32_t {
            uint32_t x = 0;
            for (size_t i = 0; i < len; ++i)
                x = x * 10u + (uint32_t)(str[from + i] - '0');
            return x;
        };
        *this = integer((int64_t)take(pos, first));
        pos += first;
        while (pos < str.size())
        {
            *this = mul(*this, (int_fast64_t)DEC);
            *this += integer((int64_t)take(pos, 9));
            pos += 9;
        }
        if (sign < 0 && !limbs.empty())
            this->sign = -1;
        trim();
    }

    static uint32_t mont_n0inv(uint32_t n0)
    {
        uint32_t x = n0;
        x *= 2u - n0 * x;
        x *= 2u - n0 * x;
        x *= 2u - n0 * x;
        x *= 2u - n0 * x;
        x *= 2u - n0 * x;
        return (uint32_t)(0u - x);
    }

    static integer mont_redc(integer T, const integer& n, uint32_t n0inv, size_t L)
    {
        T.limbs.resize(2 * L + 2, 0);
        for (size_t i = 0; i < L; ++i)
        {
            uint32_t m = T.limbs[i] * n0inv;
            uint64_t carry = 0;
            for (size_t j = 0; j < n.limbs.size(); ++j)
            {
                uint64_t cur = (uint64_t)T.limbs[i + j]
                    + (uint64_t)m * n.limbs[j] + carry;
                T.limbs[i + j] = (uint32_t)cur;
                carry = cur >> 32;
            }
            size_t k = i + n.limbs.size();
            while (carry)
            {
                if (k >= T.limbs.size())
                    T.limbs.push_back(0);
                uint64_t cur = (uint64_t)T.limbs[k] + carry;
                T.limbs[k] = (uint32_t)cur;
                carry = cur >> 32;
                ++k;
            }
        }
        integer r;
        r.sign = 1;
        if (T.limbs.size() > L)
            r.limbs.assign(T.limbs.begin() + (std::ptrdiff_t)L, T.limbs.end());
        r.trim();
        if (cmp_abs(r, n) >= 0)
            r = sub(r, n);
        return r;
    }

    static bool to_u64(const integer& n, uint64_t& out)
    {
        if (n.sign < 0)
            return false;
        if (n.limbs.empty())
        {
            out = 0;
            return true;
        }
        if (n.limbs.size() > 2)
            return false;
        if (n.limbs.size() == 1)
        {
            out = n.limbs[0];
            return true;
        }
        out = (uint64_t)n.limbs[0] | ((uint64_t)n.limbs[1] << 32);
        return true;
    }

    integer shift_left_bits(uint64_t bits) const
    {
        if (limbs.empty() || bits == 0)
            return *this;
        const size_t limb_shift = (size_t)(bits / 32);
        const int bit_shift = (int)(bits % 32);
        integer r;
        r.sign = sign;
        r.limbs.assign(limb_shift, 0);
        if (bit_shift == 0)
        {
            r.limbs.insert(r.limbs.end(), limbs.begin(), limbs.end());
            return r;
        }
        uint32_t carry = 0;
        for (size_t i = 0; i < limbs.size(); ++i)
        {
            uint64_t cur = ((uint64_t)limbs[i] << bit_shift) | carry;
            r.limbs.push_back((uint32_t)cur);
            carry = (uint32_t)(cur >> 32);
        }
        if (carry)
            r.limbs.push_back(carry);
        r.trim();
        return r;
    }

    integer shift_right_bits(uint64_t bits) const
    {
        if (limbs.empty() || bits == 0)
            return *this;
        const size_t limb_shift = (size_t)(bits / 32);
        const int bit_shift = (int)(bits % 32);
        if (limb_shift >= limbs.size())
            return integer();
        integer r;
        r.sign = sign;
        if (bit_shift == 0)
        {
            r.limbs.assign(limbs.begin() + (std::ptrdiff_t)limb_shift, limbs.end());
            r.trim();
            return r;
        }
        uint32_t carry = 0;
        for (int i = (int)limbs.size() - 1; i >= (int)limb_shift; --i)
        {
            uint32_t cur = (limbs[(size_t)i] >> bit_shift) | carry;
            carry = (uint32_t)((uint64_t)limbs[(size_t)i] << (32 - bit_shift));
            r.limbs.push_back(cur);
        }
        std::reverse(r.limbs.begin(), r.limbs.end());
        r.trim();
        return r;
    }

    static integer newton_reciprocal(const integer& v, int p)
    {
        const int vb = (int)bit_length_abs(v);
        uint32_t top = v.limbs.back();
        if (v.limbs.size() >= 2 && vb > 32)
        {
            int extra = vb - 32;
            integer t = v.shift_right_bits((uint64_t)extra);
            top = t.limbs.empty() ? 1u : t.limbs[0];
        }
        if (top == 0)
            top = 1;
        integer x((int64_t)((uint64_t(1) << 32) / top));
        int k = vb;
        int guard = 0;
        while (k < p && guard < 64)
        {
            int nk = k * 2;
            if (nk > p + 2)
                nk = p + 2;
            integer inner = (integer(1) << (k + 1)) - v * x;
            if (inner.sign < 0)
                inner = integer();
            integer x2k = x * inner;
            if (2 * k > nk)
                x = x2k.shift_right_bits((uint64_t)(2 * k - nk));
            else if (nk > 2 * k)
                x = x2k.shift_left_bits((uint64_t)(nk - 2 * k));
            else
                x = x2k;
            k = nk;
            ++guard;
        }
        if (k > p)
            x = x.shift_right_bits((uint64_t)(k - p));
        else if (k < p)
            x = x.shift_left_bits((uint64_t)(p - k));
        return x;
    }

    static qr<integer> newton_div(const integer& lhs, const integer& rhs)
    {
        integer u = abs(lhs);
        integer v = abs(rhs);
        const int ub = (int)bit_length_abs(u);
        const int p = ub + 2;
        integer rec = newton_reciprocal(v, p);
        integer q = (u * rec).shift_right_bits((uint64_t)p);
        integer prod = q * v;
        int fix = 0;
        while (prod > u && fix < 8)
        {
            q -= 1;
            prod = q * v;
            ++fix;
        }
        integer diff = u - prod;
        while (diff >= v && fix < 16)
        {
            q += 1;
            prod = q * v;
            diff = u - prod;
            ++fix;
        }
        if (prod > u || diff >= v)
            return AlgoD(lhs, rhs);
        qr<integer> res;
        q.sign = lhs.sign * rhs.sign;
        if (q.limbs.empty())
            q.sign = 1;
        integer r = u - prod;
        r.sign = lhs.sign;
        r.trim();
        if (r.limbs.empty())
            r.sign = 1;
        q.trim();
        res.q = q;
        res.r = r;
        return res;
    }

    static int64_t bit_length_abs(const integer& num)
    {
        if (num.limbs.empty())
            return 0;
        int64_t bits = (int64_t)(num.limbs.size() - 1) * 32;
        uint32_t m = num.limbs.back();
        bits += 32 - clz32(m);
        return bits;
    }

    static int jacobi(integer a, integer n)
    {
        if (n.sign < 0 || n.limbs.empty() || n.even())
            return 0;
        int t = 1;
        a %= n;
        if (a.sign < 0)
            a += n;
        while (!a.limbs.empty())
        {
            while (a.even())
            {
                a = a.shift_right_bits(1);
                uint32_t r = n.limbs[0] & 7u;
                if (r == 3 || r == 5)
                    t = -t;
            }
            std::swap(a, n);
            if ((a.limbs[0] & 3u) == 3 && (n.limbs[0] & 3u) == 3)
                t = -t;
            a %= n;
        }
        return n == 1 ? t : 0;
    }

    static bool miller_rabin_base(const integer& n, const integer& a)
    {
        integer nm1 = n - 1;
        integer d = nm1;
        int s = 0;
        while (d.even())
        {
            d = d.shift_right_bits(1);
            ++s;
        }
        integer x = modexp(a, d, n);
        if (x == 1 || x == nm1)
            return true;
        for (int i = 1; i < s; ++i)
        {
            x = (x * x) % n;
            if (x == nm1)
                return true;
            if (x == 1)
                return false;
        }
        return false;
    }

    static integer norm_mod(integer x, const integer& n)
    {
        x %= n;
        if (x.sign < 0)
            x += n;
        return x;
    }

    friend struct primality;
    friend integer prime(int k, primality check);

    static bool lucas_strong(const integer& n)
    {
        int D = 5;
        int sgn = 1;
        integer Dd;
        for (int it = 0; it < 64; ++it)
        {
            Dd = integer(sgn * D);
            int j = jacobi(Dd, n);
            if (j == -1)
                break;
            if (j == 0 && abs(Dd) != n)
                return false;
            D += 2;
            sgn = -sgn;
            if (it == 20)
            {
                integer s = sqrt(n);
                if (s * s == n)
                    return false;
            }
            if (it == 63)
                return false;
        }
        integer P = 1;
        integer Q = (integer(1) - Dd) / 4;
        integer inv2 = (n + 1) / 2;
        integer k = n + 1;
        integer d = k;
        int s = 0;
        while (d.even())
        {
            d = d.shift_right_bits(1);
            ++s;
        }
        auto lucas_at = [&](const integer& idx, integer& U, integer& V) {
            U = 1;
            V = P % n;
            integer Qk = Q % n;
            if (idx == 1)
                return;
            integer bits = idx;
            int bl = (int)bit_length_abs(bits);
            for (int i = bl - 2; i >= 0; --i)
            {
                U = norm_mod(U * V, n);
                V = norm_mod(V * V - Qk * 2, n);
                Qk = norm_mod(Qk * Qk, n);
                bool bit = false;
                if ((int)bits.limbs.size() > i / 32)
                    bit = (bits.limbs[(size_t)(i / 32)] >> (i % 32)) & 1u;
                if (bit)
                {
                    integer U2 = norm_mod((P * U + V) * inv2, n);
                    integer V2 = norm_mod((Dd * U + P * V) * inv2, n);
                    U = U2;
                    V = V2;
                    Qk = norm_mod(Qk * Q, n);
                }
            }
        };
        integer U, V;
        lucas_at(d, U, V);
        if (U == 0 || V == 0)
            return true;
        for (int r = 1; r < s; ++r)
        {
            V = norm_mod(V * V - Q * 2, n);
            if (V == 0)
                return true;
        }
        return false;
    }

    static integer extract_bits(const integer& x, int64_t start, int64_t len)
    {
        if (len <= 0 || x.limbs.empty())
            return integer();
        integer t = start > 0 ? x.shift_right_bits((uint64_t)start) : x;
        integer hi = t.shift_right_bits((uint64_t)len);
        return t - hi.shift_left_bits((uint64_t)len);
    }

    static integer fermat_reduce(integer x, int N)
    {
        if (N <= 0)
            return integer();
        integer twoN = integer(1).shift_left_bits((uint64_t)N);
        integer fermat = twoN + 1;
        for (;;)
        {
            if (x.sign < 0)
                x += fermat;
            if (cmp_abs(x, twoN) < 0)
            {
                x.sign = 1;
                return x;
            }
            if (x == twoN)
                return twoN;
            integer hi = x.shift_right_bits((uint64_t)N);
            integer lo = x - hi.shift_left_bits((uint64_t)N);
            x = lo - hi;
        }
    }

    static integer fermat_shift(integer x, int64_t s, int N)
    {
        x = fermat_reduce(std::move(x), N);
        int64_t period = 2 * (int64_t)N;
        s %= period;
        if (s < 0)
            s += period;
        if (s == 0)
            return x;
        if (s >= N)
        {
            x = fermat_reduce(-x, N);
            s -= N;
            if (s == 0)
                return x;
        }
        return fermat_reduce(x.shift_left_bits((uint64_t)s), N);
    }

    static void fermat_fft(std::vector<integer>& a, int N, bool invert)
    {
        const size_t n = a.size();
        for (size_t i = 1, j = 0; i < n; ++i)
        {
            size_t bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a[i], a[j]);
        }
        for (size_t len = 2; len <= n; len <<= 1)
        {
            int64_t wlen = (int64_t)(2 * N / (int)len);
            if (invert)
                wlen = 2 * (int64_t)N - wlen;
            for (size_t i = 0; i < n; i += len)
            {
                int64_t w = 0;
                for (size_t j = 0; j < len / 2; ++j)
                {
                    integer u = a[i + j];
                    integer v = fermat_shift(a[i + j + len / 2], w, N);
                    a[i + j] = fermat_reduce(u + v, N);
                    a[i + j + len / 2] = fermat_reduce(u - v, N);
                    w += wlen;
                    if (w >= 2 * (int64_t)N)
                        w -= 2 * (int64_t)N;
                }
            }
        }
        if (invert)
        {
            int m = 0;
            for (size_t t = n; t > 1; t >>= 1)
                ++m;
            int64_t invshift = 2 * (int64_t)N - m;
            for (size_t i = 0; i < n; ++i)
                a[i] = fermat_shift(a[i], invshift, N);
        }
    }

public:
    friend integer mul(const integer lhs, const int_fast64_t rhs)
    {
        if (rhs == 0 || lhs.limbs.empty())
            return integer();
        integer res = lhs;
        uint64_t multiplier;
        if (rhs < 0)
        {
            res.sign = -res.sign;
            multiplier = (uint64_t)(-(rhs + 1)) + 1;
        }
        else
            multiplier = (uint64_t)rhs;
        if (multiplier > 0xFFFFFFFFull)
            return lhs * integer(rhs);
        uint64_t carry = 0;
        for (size_t i = 0; i < res.limbs.size() || carry; ++i)
        {
            if (i == res.limbs.size())
                res.limbs.push_back(0);
            uint64_t cur = (uint64_t)res.limbs[i] * (uint32_t)multiplier + carry;
            res.limbs[i] = (uint32_t)cur;
            carry = cur >> 32;
        }
        res.trim();
        return res;
    }

    /**
     * Karatsuba multiplication (schoolbook at 32 limbs).
     */
    friend integer karatsuba(const integer x, const integer y)
    {
        integer a = abs(x);
        integer b = abs(y);
        size_t n = std::max(a.limbs.size(), b.limbs.size());
        size_t mmin = std::min(a.limbs.size(), b.limbs.size());
        if (n <= 32 || mmin <= 16)
        {
            integer r = schoolbook(a, b);
            r.sign = x.sign * y.sign;
            if (r.limbs.empty())
                r.sign = 1;
            return r;
        }
        size_t m = (n + 1) / 2;
        integer a0 = a.low_limbs(m);
        integer a1 = a.high_limbs(m);
        integer b0 = b.low_limbs(m);
        integer b1 = b.high_limbs(m);
        integer z0 = karatsuba(a0, b0);
        integer z2 = karatsuba(a1, b1);
        integer z1 = karatsuba(a0 + a1, b0 + b1) - z0 - z2;
        integer res = z0 + shl_limbs(z1, m) + shl_limbs(z2, m + m);
        res.sign = x.sign * y.sign;
        if (res.limbs.empty())
            res.sign = 1;
        return res;
    }

    friend integer schonhage_strassen(const integer x, const integer y);

    static bool ntt_length_ok(const integer& x, const integer& y)
    {
        size_t need = x.limbs.size() * 2 + y.limbs.size() * 2 + 2;
        size_t n = 1;
        while (n < need)
        {
            if (n > (size_t(1) << 22))
                return false;
            n <<= 1;
        }
        return n <= (size_t(1) << 23);
    }

    /**
     * Cooley–Tukey NTT multiply: 16-bit digits, three 31-bit primes, CRT.
     */
    friend integer ntt(const integer x, const integer y)
    {
        if (x.limbs.empty() || y.limbs.empty())
            return integer();
        std::vector<uint32_t> da, db;
        da.reserve(x.limbs.size() * 2);
        db.reserve(y.limbs.size() * 2);
        for (uint32_t w : x.limbs)
        {
            da.push_back(w & 0xFFFFu);
            da.push_back(w >> 16);
        }
        for (uint32_t w : y.limbs)
        {
            db.push_back(w & 0xFFFFu);
            db.push_back(w >> 16);
        }
        while (!da.empty() && da.back() == 0)
            da.pop_back();
        while (!db.empty() && db.back() == 0)
            db.pop_back();
        size_t need = da.size() + db.size();
        size_t n = 1;
        while (n < need)
            n <<= 1;
        /* 998244353 = 119·2^23+1, g=3; max length 2^23 */
        if (n > (size_t(1) << 23))
            return schonhage_strassen(x, y);

        constexpr uint32_t P1 = 998244353u;
        constexpr uint32_t P2 = 897581057u;
        constexpr uint32_t P3 = 754974721u;
        constexpr uint32_t G1 = 3u;
        constexpr uint32_t G2 = 3u;
        constexpr uint32_t G3 = 11u;

        auto conv = [&](uint32_t mod, uint32_t gen) {
            std::vector<uint32_t> a(n, 0), b(n, 0);
            std::copy(da.begin(), da.end(), a.begin());
            std::copy(db.begin(), db.end(), b.begin());
            ntt_transform(a, mod, gen, false);
            ntt_transform(b, mod, gen, false);
            for (size_t i = 0; i < n; ++i)
                a[i] = (uint32_t)((uint64_t)a[i] * b[i] % mod);
            ntt_transform(a, mod, gen, true);
            return a;
        };

        std::vector<uint32_t> c1 = conv(P1, G1);
        std::vector<uint32_t> c2 = conv(P2, G2);
        std::vector<uint32_t> c3 = conv(P3, G3);

        const uint32_t inv12 = ntt_mod_inv(P1 % P2, P2);
        const uint64_t p1p2 = (uint64_t)P1 * P2;
        const uint32_t inv123 = ntt_mod_inv((uint32_t)(p1p2 % P3), P3);

        std::vector<uint64_t> digits(n + 8, 0);
        for (size_t i = 0; i < n; ++i)
        {
            uint64_t x1 = c1[i];
            uint64_t t2 = (uint64_t)((c2[i] + P2 - (uint32_t)(x1 % P2)) % P2) * inv12 % P2;
            uint64_t a = x1 + t2 * P1;
            uint64_t a_mod = a % P3;
            uint64_t t3 = (uint64_t)((c3[i] + P3 - (uint32_t)a_mod) % P3) * inv123 % P3;
            uint64_t lo, hi;
            mul_u64(t3, p1p2, lo, hi);
            uint64_t s = lo + a;
            if (s < lo)
                ++hi;
            /* hi:lo is the CRT digit; store as 96-bit then fold below */
            digits[i] += s;
            if (digits[i] < s)
                ++hi;
            if (hi)
            {
                /* hi < 2^32 in practice for our size; spill into higher 16-bit slots */
                uint64_t spill = hi;
                size_t k = i + 4; /* 4 × 16-bit = 64-bit */
                while (spill)
                {
                    if (k >= digits.size())
                        digits.resize(k + 4, 0);
                    uint64_t nxt = digits[k] + (spill & 0xFFFFu);
                    digits[k] = nxt;
                    spill = (spill >> 16) + (nxt >> 16);
                    digits[k] &= 0xFFFFu;
                    /* actually handle carry properly below */
                    ++k;
                    break;
                }
                /* simpler: add hi << 64 as +hi at digit i+4 in base 2^16 */
                size_t pos = i + 4;
                if (pos >= digits.size())
                    digits.resize(pos + 2, 0);
                digits[pos] += hi;
            }
        }
        uint64_t carry = 0;
        for (size_t i = 0; i < digits.size(); ++i)
        {
            uint64_t cur = digits[i] + carry;
            digits[i] = cur & 0xFFFFu;
            carry = cur >> 16;
        }
        while (carry)
        {
            digits.push_back(carry & 0xFFFFu);
            carry >>= 16;
        }
        integer res;
        res.sign = x.sign * y.sign;
        res.limbs.reserve((digits.size() + 1) / 2);
        for (size_t i = 0; i < digits.size(); i += 2)
        {
            uint32_t lo = (uint32_t)digits[i];
            uint32_t hi = (i + 1 < digits.size()) ? (uint32_t)digits[i + 1] : 0;
            res.limbs.push_back(lo | (hi << 16));
        }
        res.trim();
        return res;
    }

    /**
     * Schönhage–Strassen: FFT over Z/(2^N+1)Z (Fermat ring). Twiddles are shifts.
     *
     * A. Schönhage, V. Strassen. Schnelle Multiplikation großer Zahlen.
     * Computing 7 (1971) 281–292.
     */
    friend integer schonhage_strassen(const integer x, const integer y)
    {
        if (x.limbs.empty() || y.limbs.empty())
            return integer();
        integer a = abs(x);
        integer b = abs(y);
        const int64_t ba = bit_length_abs(a);
        const int64_t bb = bit_length_abs(b);
        const int64_t outbits = ba + bb;
        if (outbits <= 64)
        {
            integer r = schoolbook(a, b);
            r.sign = x.sign * y.sign;
            if (r.limbs.empty())
                r.sign = 1;
            return r;
        }

        int64_t M = 4;
        while (M * M < outbits)
            M <<= 1;

        int64_t Lmax = M / 2;
        int64_t inbits = (std::max(ba, bb) + Lmax - 1) / Lmax;
        if (inbits < 1)
            inbits = 1;
        int64_t La = (ba + inbits - 1) / inbits;
        int64_t Lb = (bb + inbits - 1) / inbits;
        int64_t L = std::max(La, Lb);
        if (L < 1)
            L = 1;
        while (M < 2 * L)
            M <<= 1;

        int logM = 0;
        for (int64_t t = M; t > 1; t >>= 1)
            ++logM;

        int64_t Nmin = 2 * inbits + logM + 2;
        int64_t N = ((Nmin + M - 1) / M) * M;
        if (N < M)
            N = M;
        while (N >= outbits && M < (int64_t(1) << 20))
        {
            M <<= 1;
            Lmax = M / 2;
            inbits = (std::max(ba, bb) + Lmax - 1) / Lmax;
            if (inbits < 1)
                inbits = 1;
            La = (ba + inbits - 1) / inbits;
            Lb = (bb + inbits - 1) / inbits;
            L = std::max(La, Lb);
            while (M < 2 * L)
                M <<= 1;
            logM = 0;
            for (int64_t t = M; t > 1; t >>= 1)
                ++logM;
            Nmin = 2 * inbits + logM + 2;
            N = ((Nmin + M - 1) / M) * M;
        }
        if (N > 0x7FFFFFFF)
            return karatsuba(x, y);

        std::vector<integer> fa((size_t)M), fb((size_t)M);
        for (int64_t i = 0; i < M; ++i)
        {
            fa[(size_t)i] = fermat_reduce(extract_bits(a, i * inbits, inbits), (int)N);
            fb[(size_t)i] = fermat_reduce(extract_bits(b, i * inbits, inbits), (int)N);
        }
        fermat_fft(fa, (int)N, false);
        fermat_fft(fb, (int)N, false);
        for (int64_t i = 0; i < M; ++i)
            fa[(size_t)i] = fermat_reduce(fa[(size_t)i] * fb[(size_t)i], (int)N);
        fermat_fft(fa, (int)N, true);

        integer res;
        for (int64_t i = 0; i < M; ++i)
        {
            integer c = fa[(size_t)i];
            if (c == integer(1).shift_left_bits((uint64_t)N))
                c = -integer(1);
            if (!c.limbs.empty())
                res += c.shift_left_bits((uint64_t)i * (uint64_t)inbits);
        }
        res.sign = x.sign * y.sign;
        if (res.limbs.empty())
            res.sign = 1;
        return res;
    }

    friend int_fast64_t mod(const integer& lhs, const int_fast64_t& rhs)
    {
        if (rhs == 0)
            throw std::runtime_error("modulo by zero");
        int_fast64_t den = rhs < 0 ? -rhs : rhs;
        if (den > 0x7FFFFFFF)
            return (int_fast64_t)(lhs % integer(rhs));
        uint64_t rem = 0;
        for (int i = (int)lhs.limbs.size() - 1; i >= 0; --i)
            rem = (((rem << 32) | lhs.limbs[(size_t)i]) % (uint64_t)den);
        int_fast64_t r = (int_fast64_t)rem;
        return r * lhs.sign;
    }

    friend integer div(const integer& lhs, const int_fast64_t& rhs)
    {
        if (rhs == 0)
            throw std::runtime_error("division by zero");
        integer num = lhs;
        int_fast64_t den = rhs;
        if (den < 0)
        {
            den = -den;
            num.sign = -num.sign;
        }
        if ((uint64_t)den > 0xFFFFFFFFull)
            return AlgoD(lhs, integer(rhs)).q;
        uint64_t rem = 0;
        for (int i = (int)num.limbs.size() - 1; i >= 0; --i)
        {
            uint64_t cur = (rem << 32) | num.limbs[(size_t)i];
            num.limbs[(size_t)i] = (uint32_t)(cur / (uint64_t)den);
            rem = cur % (uint64_t)den;
        }
        num.trim();
        return num;
    }

    /**
     * Donald Knuth's Algorithm D (TAOCP Vol. 2, Seminumerical Algorithms, §4.3.1).
     */
    friend qr<integer> AlgoD(const integer& lhs, const integer& rhs)
    {
        qr<integer> res;
        if (rhs.limbs.empty())
            throw std::runtime_error("division by zero");
        integer u = abs(lhs);
        integer v = abs(rhs);
        if (cmp_abs(u, v) < 0)
        {
            res.q = integer();
            res.r = lhs;
            return res;
        }
        if (v.limbs.size() == 1)
        {
            uint32_t d = v.limbs[0];
            uint64_t rem = 0;
            integer q;
            q.limbs.resize(u.limbs.size());
            for (int i = (int)u.limbs.size() - 1; i >= 0; --i)
            {
                uint64_t cur = (rem << 32) | u.limbs[(size_t)i];
                q.limbs[(size_t)i] = (uint32_t)(cur / d);
                rem = cur % d;
            }
            q.sign = lhs.sign * rhs.sign;
            q.trim();
            res.q = q;
            res.r = integer((int64_t)rem);
            res.r.sign = lhs.sign;
            if (res.r.limbs.empty())
                res.r.sign = 1;
            return res;
        }
        int lz = clz32(v.limbs.back());
        u = u.shift_left_bits((uint64_t)lz);
        v = v.shift_left_bits((uint64_t)lz);
        const size_t n = v.limbs.size();
        u.limbs.push_back(0);
        const size_t m = u.limbs.size() - n - 1;
        integer q;
        q.limbs.assign(m + 1, 0);
        const uint32_t vn1 = v.limbs[n - 1];
        const uint32_t vn2 = n >= 2 ? v.limbs[n - 2] : 0;
        for (int j = (int)m; j >= 0; --j)
        {
            uint64_t u_hi = ((uint64_t)u.limbs[(size_t)j + n] << 32) | u.limbs[(size_t)j + n - 1];
            uint64_t qhat;
            uint64_t rhat;
            if (u.limbs[(size_t)j + n] == vn1)
            {
                qhat = 0xFFFFFFFFull;
                rhat = (uint64_t)u.limbs[(size_t)j + n - 1] + vn1;
            }
            else
            {
                qhat = u_hi / vn1;
                rhat = u_hi % vn1;
            }
            uint32_t u_n2 = (j + (int)n >= 2) ? u.limbs[(size_t)j + n - 2] : 0;
            while (true)
            {
                uint64_t left_lo, left_hi;
                mul_u64(qhat, vn2, left_lo, left_hi);
                uint64_t right_hi = rhat >> 32;
                uint64_t right_lo = (rhat << 32) | u_n2;
                if (left_hi < right_hi || (left_hi == right_hi && left_lo <= right_lo))
                    break;
                --qhat;
                rhat += vn1;
                if (rhat >= (1ull << 32))
                    break;
            }
            uint64_t carry = 0;
            uint32_t borrow = 0;
            for (size_t i = 0; i < n; ++i)
            {
                uint64_t prod = (uint64_t)v.limbs[i] * qhat + carry;
                carry = prod >> 32;
                uint32_t plo = (uint32_t)prod;
                uint64_t left = u.limbs[(size_t)j + i];
                uint64_t right = (uint64_t)plo + borrow;
                if (left < right)
                {
                    u.limbs[(size_t)j + i] = (uint32_t)(left + (1ull << 32) - right);
                    borrow = 1;
                }
                else
                {
                    u.limbs[(size_t)j + i] = (uint32_t)(left - right);
                    borrow = 0;
                }
            }
            uint64_t left = u.limbs[(size_t)j + n];
            uint64_t right = carry + borrow;
            bool over = left < right;
            u.limbs[(size_t)j + n] = (uint32_t)(left - right);
            if (over)
            {
                uint32_t c = 0;
                for (size_t i = 0; i < n; ++i)
                {
                    uint64_t s = (uint64_t)u.limbs[(size_t)j + i] + v.limbs[i] + c;
                    u.limbs[(size_t)j + i] = (uint32_t)s;
                    c = (uint32_t)(s >> 32);
                }
                u.limbs[(size_t)j + n] = (uint32_t)((uint64_t)u.limbs[(size_t)j + n] + c);
                --qhat;
            }
            q.limbs[(size_t)j] = (uint32_t)qhat;
        }
        q.sign = lhs.sign * rhs.sign;
        q.trim();
        u.limbs.resize(n + 1);
        u.trim();
        integer r = u.shift_right_bits((uint64_t)lz);
        r.sign = lhs.sign;
        r.trim();
        if (r.limbs.empty())
            r.sign = 1;
        res.q = q;
        res.r = r;
        return res;
    }

    friend integer sub(const integer& lhs, const integer& rhs)
    {
        if (lhs.sign == rhs.sign)
        {
            if (cmp_abs(lhs, rhs) >= 0)
            {
                integer res = lhs;
                uint32_t borrow = 0;
                size_t n = std::max(lhs.limbs.size(), rhs.limbs.size());
                res.limbs.resize(n, 0);
                for (size_t i = 0; i < n; ++i)
                {
                    uint64_t left = res.limbs[i];
                    uint64_t right = (uint64_t)(i < rhs.limbs.size() ? rhs.limbs[i] : 0) + borrow;
                    if (left < right)
                    {
                        res.limbs[i] = (uint32_t)(left + (1ull << 32) - right);
                        borrow = 1;
                    }
                    else
                    {
                        res.limbs[i] = (uint32_t)(left - right);
                        borrow = 0;
                    }
                }
                res.trim();
                return res;
            }
            integer t = sub(rhs, lhs);
            t.sign = -t.sign;
            if (t.limbs.empty())
                t.sign = 1;
            return t;
        }
        return add(lhs, -rhs);
    }

    friend integer add(const integer& lhs, const integer& rhs)
    {
        if (lhs.sign == rhs.sign)
        {
            integer res = rhs;
            size_t n = std::max(lhs.limbs.size(), rhs.limbs.size());
            res.limbs.resize(n, 0);
            uint64_t carry = 0;
            for (size_t i = 0; i < n || carry; ++i)
            {
                if (i == res.limbs.size())
                    res.limbs.push_back(0);
                uint64_t cur = (uint64_t)res.limbs[i] + carry
                    + (i < lhs.limbs.size() ? lhs.limbs[i] : 0);
                res.limbs[i] = (uint32_t)cur;
                carry = cur >> 32;
            }
            return res;
        }
        return sub(lhs, -rhs);
    }

    integer(): sign(1) {}
    integer(const char* str)
    {
        read(str);
    }
    integer(const std::string& str)
    {
        read(str);
    }
    integer(const integer& num):
        limbs(num.limbs), sign(num.sign)
    {}
    integer(integer&& num) noexcept:
        limbs(std::move(num.limbs)), sign(num.sign)
    {
        num.sign = 1;
    }
    template<typename int_t>
    integer(const int_t& num)
    {
        static_assert(std::is_integral<int_t>::value, "num is non-integral.");
        *this = num;
    }

    integer& operator=(const char* str)
    {
        read(str);
        return *this;
    }
    integer& operator=(const std::string& str)
    {
        read(str);
        return *this;
    }
    integer& operator=(const integer& num)
    {
        if (this != &num)
        {
            sign = num.sign;
            limbs = num.limbs;
        }
        return *this;
    }
    integer& operator=(integer&& num) noexcept
    {
        if (this != &num)
        {
            sign = num.sign;
            limbs = std::move(num.limbs);
            num.sign = 1;
        }
        return *this;
    }
    template<typename int_t>
    integer& operator=(const int_t rhs)
    {
        static_assert(std::is_integral<int_t>::value, "digits is non-integral.");
        limbs.clear();
        sign = 1;
        if (rhs == 0)
            return *this;
        using u_t = typename std::make_unsigned<typename std::conditional<
            std::is_same<int_t, bool>::value, unsigned, int_t>::type>::type;
        u_t mag;
        if (std::is_signed<int_t>::value && rhs < 0)
        {
            sign = -1;
            int_t tmp = rhs;
            mag = (u_t)(-(tmp + 1)) + 1;
        }
        else
            mag = (u_t)rhs;
        uint64_t v = (uint64_t)mag;
        while (v)
        {
            limbs.push_back((uint32_t)v);
            v >>= 32;
        }
        return *this;
    }

    integer operator+(const integer& rhs) const
    {
        return add(*this, rhs);
    }
    template<typename int_t,
             typename = std::enable_if_t<std::is_arithmetic<int_t>::value>>
    integer operator+(const int_t& rhs) const
    {
        return *this + integer(rhs);
    }
    template<typename int_t>
    friend integer operator+(const int_t lhs, const integer& rhs)
    {
        return integer(lhs) + rhs;
    }
    integer& operator+=(const integer& rhs)
    {
        *this = *this + rhs;
        return *this;
    }
    integer operator++()
    {
        *this += 1;
        return *this;
    }
    integer operator++(int)
    {
        integer tmp = *this;
        *this += 1;
        return tmp;
    }
    integer operator-(const integer& rhs) const
    {
        return sub(*this, rhs);
    }
    template<typename int_t,
             typename = std::enable_if_t<std::is_arithmetic<int_t>::value>>
    integer operator-(const int_t& rhs) const
    {
        return *this - integer(rhs);
    }
    template<typename int_t>
    friend integer operator-(const int_t lhs, const integer& rhs)
    {
        return integer(lhs) - rhs;
    }
    integer& operator-=(const integer& rhs)
    {
        *this = sub(*this, rhs);
        return *this;
    }
    integer operator--()
    {
        *this -= 1;
        return *this;
    }
    integer operator--(int)
    {
        integer tmp = *this;
        *this -= 1;
        return tmp;
    }
    integer operator*(const integer& rhs) const
    {
        const size_t n = std::max(limbs.size(), rhs.limbs.size());
        if (n <= 32)
            return schoolbook(*this, rhs);
        if (n <= 256)
            return karatsuba(*this, rhs);
        if (n <= 4096 && ntt_length_ok(*this, rhs))
            return ntt(*this, rhs);
        return schonhage_strassen(*this, rhs);
    }
    integer& operator*=(const integer& rhs)
    {
        *this = *this * rhs;
        return *this;
    }
    integer operator/(const integer& rhs) const
    {
        if (rhs.limbs.size() >= 64 && limbs.size() >= 64)
            return newton_div(*this, rhs).q;
        return AlgoD(*this, rhs).q;
    }
    integer& operator/=(const integer& rhs)
    {
        *this = *this / rhs;
        return *this;
    }
    integer operator%(const integer& rhs) const
    {
        if (rhs.limbs.size() >= 64 && limbs.size() >= 64)
            return newton_div(*this, rhs).r;
        return AlgoD(*this, rhs).r;
    }
    integer& operator%=(const integer rhs)
    {
        *this = *this % rhs;
        return *this;
    }
    template<typename int_t>
    integer operator*(const int_t& rhs) const
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        return mul(*this, (int_fast64_t)rhs);
    }
    template<typename int_t>
    friend integer operator*(const int_t lhs, const integer& rhs)
    {
        return rhs * lhs;
    }
    template<typename int_t>
    integer& operator*=(int_t rhs)
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        *this = mul(*this, (int_fast64_t)rhs);
        return *this;
    }
    template<typename int_t>
    integer operator/(const int_t& rhs) const
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        return div(*this, (int_fast64_t)rhs);
    }
    template<typename int_t>
    integer& operator/=(const int_t& rhs)
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        *this = div(*this, (int_fast64_t)rhs);
        return *this;
    }
    template<typename int_t>
    integer operator%(const int_t& rhs) const
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        return integer(mod(*this, (int_fast64_t)rhs));
    }
    template<typename int_t>
    integer& operator%=(const int_t& rhs)
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        *this = integer(mod(*this, (int_fast64_t)rhs));
        return *this;
    }
    bool operator<(const integer& rhs) const
    {
        if (sign != rhs.sign)
            return sign < rhs.sign;
        int c = cmp_abs(*this, rhs);
        if (c == 0)
            return false;
        return (c < 0) == (sign > 0);
    }
    template<typename int_t>
    bool operator<(const int_t& rhs) const
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        return *this < integer(rhs);
    }
    template<typename int_t>
    friend bool operator<(const int_t& lhs, const integer& rhs)
    {
        static_assert(std::is_integral<int_t>::value, "rhs is non-integral.");
        return integer(lhs) < rhs;
    }
    template<typename int_t>
    bool operator>(const int_t& rhs) const
    {
        return rhs < *this;
    }
    bool operator>(const integer& rhs) const
    {
        return rhs < *this;
    }
    template<typename int_t>
    bool operator<=(const int_t& rhs) const
    {
        return !(rhs < *this);
    }
    bool operator<=(const integer& rhs) const
    {
        return !(rhs < *this);
    }
    template<typename int_t>
    bool operator>=(const int_t& rhs) const
    {
        return !(*this < rhs);
    }
    bool operator>=(const integer& rhs) const
    {
        return !(*this < rhs);
    }
    template<typename int_t>
    bool operator==(const int_t& rhs) const
    {
        return !(*this < rhs) && !(rhs < *this);
    }
    bool operator==(const integer& rhs) const
    {
        return !(*this < rhs) && !(rhs < *this);
    }
    template<typename int_t>
    bool operator!=(const int_t& rhs) const
    {
        return *this < rhs || rhs < *this;
    }
    bool operator!=(const integer& rhs) const
    {
        return *this < rhs || rhs < *this;
    }
    integer operator-() const
    {
        integer res = *this;
        if (!res.limbs.empty())
            res.sign = -sign;
        return res;
    }
    template<typename int_t>
    integer operator<<(const int_t& bits) const
    {
        integer b(bits);
        if (b.sign < 0 || b.limbs.empty())
            return *this;
        uint64_t n;
        if (!to_u64(b, n))
            throw std::runtime_error("shift amount too large");
        return shift_left_bits(n);
    }
    template<typename int_t>
    integer operator>>(const int_t& bits) const
    {
        integer b(bits);
        if (b.sign < 0 || b.limbs.empty())
            return *this;
        uint64_t n;
        if (!to_u64(b, n))
            return integer();
        return shift_right_bits(n);
    }
    template<typename int_t>
    integer& operator<<=(const int_t& bits)
    {
        *this = *this << bits;
        return *this;
    }
    template<typename int_t>
    integer& operator>>=(const int_t& bits)
    {
        *this = *this >> bits;
        return *this;
    }
    template<typename int_t,
             typename = std::enable_if_t<std::is_arithmetic<int_t>::value>>
    operator int_t() const
    {
        uint64_t res = 0;
        for (int i = (int)limbs.size() - 1; i >= 0; --i)
            res = (res << 32) | limbs[(size_t)i];
        if (sign < 0)
            return (int_t)(-(int64_t)res);
        return (int_t)res;
    }

    friend int_fast16_t size(const integer& num)
    {
        return (int_fast16_t)log10(num);
    }
    const std::string c_str() const
    {
        std::stringstream ss;
        ss << *this;
        std::string str;
        ss >> str;
        return str;
    }
    friend integer abs(const integer& num)
    {
        integer res(num);
        res.sign = 1;
        return res;
    }
    static integer rand(int_fast64_t size_f)
    {
        static XorShift64 eng(0xC0FFEEULL ^ 0x9E3779B97F4A7C15ULL);
        if (size_f <= 0)
            return integer();
        integer res;
        size_t nlimbs = (size_t)((size_f + 31) / 32);
        res.limbs.resize(nlimbs);
        for (size_t i = 0; i < nlimbs; ++i)
            res.limbs[i] = (uint32_t)eng.next();
        int top = (int)((size_f - 1) % 32);
        uint32_t mask = (top == 31) ? 0xFFFFFFFFu : ((1u << (top + 1)) - 1u);
        res.limbs.back() &= mask;
        res.limbs.back() |= (1u << top);
        res.trim();
        return res;
    }
    friend std::istream& operator>>(std::istream& stream, integer& rhs)
    {
        std::string s;
        stream >> s;
        rhs.read(s);
        return stream;
    }
    friend std::ostream& operator<<(std::ostream& stream, const integer& rhs)
    {
        if (rhs.limbs.empty())
        {
            stream << 0;
            return stream;
        }
        if (rhs.sign == -1)
            stream << '-';
        const uint32_t DEC = 1000000000u;
        integer tmp = abs(rhs);
        std::vector<uint32_t> chunks;
        while (!tmp.limbs.empty())
        {
            uint64_t rem = 0;
            for (int i = (int)tmp.limbs.size() - 1; i >= 0; --i)
            {
                uint64_t cur = (rem << 32) | tmp.limbs[(size_t)i];
                tmp.limbs[(size_t)i] = (uint32_t)(cur / DEC);
                rem = cur % DEC;
            }
            tmp.trim();
            chunks.push_back((uint32_t)rem);
        }
        stream << chunks.back();
        for (int i = (int)chunks.size() - 2; i >= 0; --i)
            stream << std::setw(9) << std::setfill('0') << chunks[(size_t)i];
        return stream;
    }
    friend int_fast64_t log10(const integer& num)
    {
        if (num.limbs.empty())
            return 0;
        integer t = abs(num);
        int_fast64_t digits = 0;
        const uint32_t DEC = 1000000000u;
        while (t.limbs.size() > 1)
        {
            t = div(t, (int_fast64_t)DEC);
            digits += 9;
        }
        uint32_t m = t.limbs.empty() ? 0 : t.limbs[0];
        while (m)
        {
            ++digits;
            m /= 10;
        }
        return digits;
    }
    friend int_fast64_t log2(const integer& num)
    {
        return bit_length_abs(num);
    }
    friend void debug(const integer& num)
    {
        std::cout << num.sign << std::endl;
        for (auto v : num.limbs)
            std::cout << v << std::endl;
    }
    friend integer min(const integer& lhs, const integer& rhs)
    {
        return lhs < rhs ? lhs : rhs;
    }
    friend integer max(const integer& lhs, const integer& rhs)
    {
        return lhs < rhs ? rhs : lhs;
    }

    /**
     * Lehmer's GCD (Euclidean fallback on a single limb).
     */
    friend integer gcd(integer lhs, integer rhs)
    {
        lhs = abs(lhs);
        rhs = abs(rhs);
        if (cmp_abs(lhs, rhs) < 0)
            std::swap(lhs, rhs);
        while (!rhs.limbs.empty() && rhs.limbs.size() > 1)
        {
            uint64_t xa = 0, xb = 0;
            size_t n = lhs.limbs.size();
            xa = lhs.limbs[n - 1];
            if (n >= 2)
                xa = (xa << 32) | lhs.limbs[n - 2];
            if (rhs.limbs.size() == n)
            {
                xb = rhs.limbs[n - 1];
                if (n >= 2)
                    xb = (xb << 32) | rhs.limbs[n - 2];
            }
            else if (rhs.limbs.size() + 1 == n && n >= 2)
                xb = rhs.limbs[n - 2];
            int sh = clz64(xa);
            xa <<= sh;
            xb <<= sh;
            int64_t A = 1, B = 0, C = 0, D = 1;
            int64_t xx = (int64_t)(xa >> 1);
            int64_t yy = (int64_t)(xb >> 1);
            bool progress = false;
            while (yy + C != 0 && yy + D != 0)
            {
                int64_t den1 = yy + C;
                int64_t den2 = yy + D;
                if (den1 == 0 || den2 == 0)
                    break;
                int64_t q1 = (xx + A) / den1;
                int64_t q2 = (xx + B) / den2;
                if (q1 != q2)
                    break;
                int64_t nA = C, nB = D, nC = A - q1 * C, nD = B - q1 * D;
                int64_t nx = yy, ny = xx - q1 * yy;
                A = nA;
                B = nB;
                C = nC;
                D = nD;
                xx = nx;
                yy = ny;
                progress = true;
            }
            if (!progress || B == 0)
            {
                integer r = lhs % rhs;
                lhs = rhs;
                rhs = r;
            }
            else
            {
                integer na = abs(lhs * A + rhs * B);
                integer nb = abs(lhs * C + rhs * D);
                lhs = na;
                rhs = nb;
            }
        }
        if (rhs.limbs.empty())
            return lhs;
        uint64_t a, b;
        if (!to_u64(lhs, a) || !to_u64(rhs, b))
        {
            while (!rhs.limbs.empty())
            {
                integer r = lhs % rhs;
                lhs = rhs;
                rhs = r;
            }
            return lhs;
        }
        while (b)
        {
            uint64_t t = a % b;
            a = b;
            b = t;
        }
        return integer((int64_t)a);
    }
    friend integer lcm(const integer& lhs, const integer& rhs)
    {
        if (lhs.limbs.empty() || rhs.limbs.empty())
            return integer();
        return abs(lhs) / gcd(lhs, rhs) * abs(rhs);
    }
    friend integer sqrt(const integer& num)
    {
        if (num.sign < 0)
            throw std::runtime_error("sqrt of negative");
        if (num.limbs.empty() || num == 1)
            return num;
        integer x = integer(1).shift_left_bits((uint64_t)((bit_length_abs(num) + 1) / 2));
        integer y = (x + num / x) / 2;
        while (y < x)
        {
            x = y;
            y = (x + num / x) / 2;
        }
        return x;
    }
    template<typename int_t>
    friend integer pow(const integer& lhs, const int_t& rhs)
    {
        integer e(rhs);
        if (e.sign < 0)
            throw std::runtime_error("negative exponent");
        integer res = 1;
        integer b = lhs;
        while (!e.limbs.empty())
        {
            if (!e.even())
                res *= b;
            e = e.shift_right_bits(1);
            if (!e.limbs.empty())
                b *= b;
        }
        return res;
    }
    friend integer modinv(const integer& lhs, const integer& rhs)
    {
        integer a = lhs % rhs;
        if (a.sign < 0)
            a += rhs;
        integer m = abs(rhs);
        integer t = 0, newt = 1;
        integer r = m, newr = a;
        while (!newr.limbs.empty())
        {
            integer q = r / newr;
            integer tmp = newt;
            newt = t - q * newt;
            t = tmp;
            tmp = newr;
            newr = r - q * newr;
            r = tmp;
        }
        if (r > 1)
            throw std::runtime_error("inverse does not exist");
        if (t.sign < 0)
            t += m;
        return t;
    }

    /**
     * Modular exponentiation — Montgomery reduction when the modulus is odd.
     */
    friend integer modexp(const integer& lhs, const integer& rhs, const integer& modv)
    {
        if (modv == 1)
            return 0;
        if (modv.limbs.empty())
            throw std::runtime_error("modulo by zero");
        integer e = rhs;
        if (e.sign < 0)
            throw std::runtime_error("negative exponent");
        integer m = abs(modv);
        integer b = lhs % m;
        if (b.sign < 0)
            b += m;
        if (m.even())
        {
            integer r = 1;
            while (!e.limbs.empty())
            {
                if (!e.even())
                    r = (r * b) % m;
                e = e.shift_right_bits(1);
                if (!e.limbs.empty())
                    b = (b * b) % m;
            }
            return r;
        }
        size_t L = m.limbs.size();
        uint32_t n0inv = mont_n0inv(m.limbs[0]);
        integer R2 = (integer(1).shift_left_bits((uint64_t)(64 * L))) % m;
        integer a_bar = mont_redc(b * R2, m, n0inv, L);
        integer r_bar = mont_redc(integer(1) * R2, m, n0inv, L);
        while (!e.limbs.empty())
        {
            if (!e.even())
                r_bar = mont_redc(r_bar * a_bar, m, n0inv, L);
            e = e.shift_right_bits(1);
            if (!e.limbs.empty())
                a_bar = mont_redc(a_bar * a_bar, m, n0inv, L);
        }
        return mont_redc(r_bar, m, n0inv, L);
    }
    friend std::vector<integer> sieve(const integer& lhs, const integer& rhs)
    {
        std::vector<integer> res;
        integer lower = lhs;
        integer upper = rhs;
        if (upper < 2 || upper < lower)
            return res;
        if (lower < 2)
            lower = 2;
        uint64_t lo = 0, hi = 0;
        if (to_u64(lower, lo) && to_u64(upper, hi) && hi - lo < 20000000ull)
        {
            std::vector<char> mark((size_t)(hi - lo + 1), 0);
            for (uint64_t i = 2; i * i <= hi; ++i)
            {
                uint64_t start = ((lo + i - 1) / i) * i;
                if (start < i * i)
                    start = i * i;
                for (uint64_t j = start; j <= hi; j += i)
                    mark[(size_t)(j - lo)] = 1;
            }
            for (uint64_t i = lo; i <= hi; ++i)
                if (!mark[(size_t)(i - lo)])
                    res.push_back(integer((int64_t)i));
            return res;
        }
        for (integer i = lower; i <= upper; i += 1)
        {
            bool ok = i >= 2;
            for (integer d = 2; d * d <= i; d += 1)
            {
                if (i % d == 0)
                {
                    ok = false;
                    break;
                }
            }
            if (ok)
                res.push_back(i);
        }
        return res;
    }
};

/**
 * Primality — small primes, then Baillie–PSW (Miller–Rabin base 2 + strong Lucas).
 */
struct primality
{
    integer k;
    primality(integer f): k(std::move(f)) {}
    bool operator()(integer n) const
    {
        if (n < 2)
            return false;
        if (n == 2 || n == 3)
            return true;
        if (n.even())
            return false;
        static const uint32_t small[] = {
            5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97
        };
        for (uint32_t p : small)
        {
            if (n == integer((int)p))
                return true;
            if (n % (int)p == 0)
                return false;
        }
        integer s = sqrt(n);
        if (s * s == n)
            return false;
        if (!integer::miller_rabin_base(n, 2))
            return false;
        return integer::lucas_strong(n);
    }
};

/* namespace-visible names for the in-class friends (ADL still works too) */
integer add(const integer& lhs, const integer& rhs);
integer sub(const integer& lhs, const integer& rhs);
integer mul(const integer lhs, const int_fast64_t rhs);
integer schoolbook(const integer& a, const integer& b);
integer karatsuba(const integer x, const integer y);
integer ntt(const integer x, const integer y);
integer schonhage_strassen(const integer x, const integer y);
integer div(const integer& lhs, const int_fast64_t& rhs);
int_fast64_t mod(const integer& lhs, const int_fast64_t& rhs);
qr<integer> AlgoD(const integer& lhs, const integer& rhs);
integer abs(const integer& num);
integer min(const integer& lhs, const integer& rhs);
integer max(const integer& lhs, const integer& rhs);
integer gcd(integer lhs, integer rhs);
integer lcm(const integer& lhs, const integer& rhs);
integer sqrt(const integer& num);
integer modinv(const integer& lhs, const integer& rhs);
integer modexp(const integer& lhs, const integer& rhs, const integer& modv);
std::vector<integer> sieve(const integer& lhs, const integer& rhs);
int_fast16_t size(const integer& num);
int_fast64_t log10(const integer& num);
int_fast64_t log2(const integer& num);
void debug(const integer& num);

inline integer prime(int k, primality check)
{
    if (k <= 1)
        return 2;
    if (k == 2)
        return 3;
    integer res = integer::rand(k);
    if (res.even())
        res += 1;
    while (!check(res))
        res += 2;
    return res;
}

} /* namespace inf */

namespace std
{
template<>
struct hash<inf::integer>
{
    size_t operator()(const inf::integer& obj) const
    {
        return hash<std::string>()(obj.c_str());
    }
};
}

#endif
