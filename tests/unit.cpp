#include "../inf.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <vector>

using inf::integer;

static void test_native_grid()
{
    integer x = "3425";
    x = "3425";
    int r1;
    int t1;
    integer t1f;
    int r2;
    int t2;
    integer t2f;
    integer a;
    integer b;
    integer c;
    int k;
    std::vector<integer> v;
    std::unordered_map<integer, integer> umap;
    for (int i = 0; i < 100; i++)
    {
        a = i;
        for (int j = 1; j < 100; j++)
        {
            r1 = rand();
            t1 = -j - r1;
            t1f = -j - r1;
            r2 = rand();
            t2 = i - r2;
            t2f = i - r2;
            assert((t1f - t2f) == integer(t1 - t2));
            t1 -= t2;
            t1f -= t2f;
            assert((t1f) == integer(t1));
            assert(inf::max(t1f, t2f) == integer(std::max(t1, t2)));
            assert(inf::min(t1f, t2f) == integer(std::min(t1, t2)));
            assert(inf::min(t1f, t2f) <= integer(std::max(t1, t2)));
            b = j;
            v.push_back(b);
            assert(v[v.size() - 1] == b);
            umap.insert({b, b});
            assert(umap[b] == b);
            assert((a / b) == integer(i / j));
            assert((a * b) == integer(i * j));
            assert((a - b) == integer(i - j));
            assert((a + b) == integer(i + j));
            assert((a % b) == integer(i % j));
            assert(pow(a, 3) / pow(b, 3) == integer((int)std::pow((double)i, 3) / (int)std::pow((double)j, 3)));
            assert(pow(a, 2) * pow(b, 2) == integer((int)std::pow((double)i, 2) * (int)std::pow((double)j, 2)));
            assert(pow(a, 3) - pow(b, 3) == integer((int)std::pow((double)i, 3) - (int)std::pow((double)j, 3)));
            assert(pow(a, 3) + pow(b, 3) == integer((int)std::pow((double)i, 3) + (int)std::pow((double)j, 3)));
            assert(pow(a, 3) % pow(b, 3) == integer((int)std::pow((double)i, 3) % (int)std::pow((double)j, 3)));
            x = integer::rand(i * j);
            assert(log2(x) <= (int_fast64_t)i * j);
            assert(log10(x) <= (int_fast64_t)((i * j) * 0.30103) + 2);
        }
        assert((--a) == integer(--i));
        assert((a--) == integer(i--));
        assert((++a) == integer(++i));
        assert((a++) == integer(i++));
        assert((a << 3) == integer(i << 3));
        assert((a >> 4) == integer(i >> 4));
        assert((a * a) == integer(i * i));
        assert((a - a) == integer(i - i));
        assert((a + a) == integer(i + i));
        c = sqrt(a + 1);
        assert(pow(c, 2) <= a + 1 && pow((c + 1), 2) >= a + 1);
        k = i * (i + 2);
        a = i * (i + 2);
        a >>= 1;
        k >>= 1;
        assert(a == integer(k));
        a >>= integer(2);
        k >>= 2;
        assert(a == integer(k));
        a <<= 3;
        k <<= 3;
        assert(a == integer(k));
        a <<= integer(4);
        k <<= 4;
        assert(a == integer(k));
        a += i;
        k += i;
        assert(a == integer(k));
        a += i;
        k += i;
        assert(a == integer(k));
        a -= i;
        k -= i;
        assert(a == integer(k));
        a -= i;
        k -= i;
        assert(a == integer(k));
        a *= i;
        k *= i;
        assert(a == integer(k));
        a *= i;
        k *= i;
        assert(a == integer(k));
        a /= i + 1;
        k /= i + 1;
        assert(a == integer(k));
        a /= i + 1;
        k /= i + 1;
        assert(a == integer(k));
        a %= i + 1;
        k %= i + 1;
        assert(a == integer(k));
        a %= i + 1;
        k %= i + 1;
        assert(a == integer(k));
    }
}

static void test_gcd_modexp()
{
    assert(gcd(integer(54), integer(24)) == 6);
    assert(lcm(integer(21), integer(6)) == 42);
    assert(modinv(integer(3), integer(11)) == 4);
    assert(modexp(integer(2), integer(10), integer(17)) == 4);
    assert(modexp(integer(3), integer(5), integer(13)) == 9);
    integer g = gcd(integer("123456789123456789"), integer("987654321"));
    assert((integer("123456789123456789") % g) == 0);
}

static void test_big_shifts()
{
    integer a = 1;
    a <<= 80;
    assert(log2(a) == 81);
    a >>= 16;
    assert(a == (integer(1) << 64));
    integer b = "0";
    for (int i = 0; i < 40; ++i)
        b = (b << 1) + 1;
    assert(b == ((integer(1) << 40) - 1));
}

static std::string repeat_digit(char d, int n)
{
    return std::string((size_t)n, d);
}

static void test_ntt_multiply()
{
    /* ~300+ 32-bit limbs forces ntt */
    integer a(repeat_digit('9', 3000));
    integer b(repeat_digit('8', 3000));
    integer p = ntt(a, b);
    integer q = karatsuba(a, b);
    assert(p == q);
    assert(p / a == b);
    assert(p / b == a);
    integer c = integer::rand(9000);
    integer d = integer::rand(9000);
    integer r = c * d;
    assert(r / c == d);
    assert(r / d == c);
}

static void test_ss_multiply()
{
    integer a(repeat_digit('9', 400));
    integer b(repeat_digit('8', 400));
    integer s = schonhage_strassen(a, b);
    integer n = ntt(a, b);
    assert(s == n);
    assert(s / a == b);
    integer c = integer::rand(2000);
    integer d = integer::rand(1800);
    integer ss = schonhage_strassen(c, d);
    assert(ss / c == d);
    assert(schonhage_strassen(-c, d) == -ss);
}

int main()
{
    test_native_grid();
    test_gcd_modexp();
    test_big_shifts();
    test_ntt_multiply();
    test_ss_multiply();
    std::cout << "Success!\n";
    return 0;
}
