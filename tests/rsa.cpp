#include "../inf.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using inf::integer;

namespace {

constexpr size_t kHashLen = 32;

uint32_t rotr32(uint32_t x, int n)
{
    return (x >> n) | (x << (32 - n));
}

// FIPS 180-4, SHA-256. Constants and compression follow that spec.
std::array<uint8_t, 32> sha256(const uint8_t* data, size_t len)
{
    static const uint32_t K[64] = {
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
    };
    uint32_t h[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
    };
    std::vector<uint8_t> msg(data, data + len);
    uint64_t bitlen = (uint64_t)len * 8;
    msg.push_back(0x80);
    while ((msg.size() % 64) != 56)
        msg.push_back(0);
    for (int i = 7; i >= 0; --i)
        msg.push_back((uint8_t)((bitlen >> (i * 8)) & 0xff));
    for (size_t off = 0; off < msg.size(); off += 64)
    {
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
        {
            w[i] = ((uint32_t)msg[off + 4 * i] << 24) | ((uint32_t)msg[off + 4 * i + 1] << 16)
                | ((uint32_t)msg[off + 4 * i + 2] << 8) | (uint32_t)msg[off + 4 * i + 3];
        }
        for (int i = 16; i < 64; ++i)
        {
            uint32_t s0 = rotr32(w[i - 15], 7) ^ rotr32(w[i - 15], 18) ^ (w[i - 15] >> 3);
            uint32_t s1 = rotr32(w[i - 2], 17) ^ rotr32(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i)
        {
            uint32_t S1 = rotr32(e, 6) ^ rotr32(e, 11) ^ rotr32(e, 25);
            uint32_t ch = (e & f) ^ ((~e) & g);
            uint32_t t1 = hh + S1 + ch + K[i] + w[i];
            uint32_t S0 = rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22);
            uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t t2 = S0 + maj;
            hh = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        h[0] += a;
        h[1] += b;
        h[2] += c;
        h[3] += d;
        h[4] += e;
        h[5] += f;
        h[6] += g;
        h[7] += hh;
    }
    std::array<uint8_t, 32> out{};
    for (int i = 0; i < 8; ++i)
    {
        out[4 * i] = (uint8_t)(h[i] >> 24);
        out[4 * i + 1] = (uint8_t)(h[i] >> 16);
        out[4 * i + 2] = (uint8_t)(h[i] >> 8);
        out[4 * i + 3] = (uint8_t)h[i];
    }
    return out;
}

std::array<uint8_t, 32> sha256(const std::vector<uint8_t>& v)
{
    return sha256(v.data(), v.size());
}

integer os2ip(const std::vector<uint8_t>& bytes)
{
    integer x = 0;
    for (uint8_t b : bytes)
        x = x * 256 + (int)b;
    return x;
}

std::vector<uint8_t> i2osp(integer x, size_t k)
{
    std::vector<uint8_t> out(k, 0);
    for (int i = (int)k - 1; i >= 0; --i)
    {
        out[(size_t)i] = (uint8_t)(int)(x % 256);
        x /= 256;
    }
    return out;
}

// RFC 8017 PKCS #1 v2.2, B.2.1 — MGF1 with SHA-256.
std::vector<uint8_t> mgf1(const uint8_t* seed, size_t seed_len, size_t mask_len)
{
    std::vector<uint8_t> mask;
    mask.reserve(mask_len);
    uint32_t counter = 0;
    while (mask.size() < mask_len)
    {
        std::vector<uint8_t> buf(seed, seed + seed_len);
        buf.push_back((uint8_t)(counter >> 24));
        buf.push_back((uint8_t)(counter >> 16));
        buf.push_back((uint8_t)(counter >> 8));
        buf.push_back((uint8_t)counter);
        auto h = sha256(buf);
        size_t take = std::min(h.size(), mask_len - mask.size());
        mask.insert(mask.end(), h.begin(), h.begin() + (std::ptrdiff_t)take);
        ++counter;
    }
    return mask;
}

std::vector<uint8_t> random_bytes(size_t n)
{
    std::vector<uint8_t> out;
    out.reserve(n);
    while (out.size() < n)
    {
        integer r = integer::rand(256);
        auto chunk = i2osp(r, 32);
        size_t take = std::min(chunk.size(), n - out.size());
        out.insert(out.end(), chunk.begin(), chunk.begin() + (std::ptrdiff_t)take);
    }
    return out;
}

// RFC 8017 PKCS #1 v2.2, §7.1 — RSAES-OAEP encode / decode (SHA-256, empty label).
std::vector<uint8_t> oaep_encode(const std::string& msg, size_t k)
{
    const size_t mLen = msg.size();
    if (k < 2 * kHashLen + 2 || mLen > k - 2 * kHashLen - 2)
        throw std::runtime_error("OAEP: message too long");
    auto lHash = sha256(nullptr, 0);
    std::vector<uint8_t> db;
    db.insert(db.end(), lHash.begin(), lHash.end());
    db.insert(db.end(), k - mLen - 2 * kHashLen - 2, 0);
    db.push_back(0x01);
    db.insert(db.end(), msg.begin(), msg.end());
    auto seed = random_bytes(kHashLen);
    auto dbMask = mgf1(seed.data(), seed.size(), k - kHashLen - 1);
    std::vector<uint8_t> maskedDB(db.size());
    for (size_t i = 0; i < db.size(); ++i)
        maskedDB[i] = db[i] ^ dbMask[i];
    auto seedMask = mgf1(maskedDB.data(), maskedDB.size(), kHashLen);
    std::vector<uint8_t> maskedSeed(kHashLen);
    for (size_t i = 0; i < kHashLen; ++i)
        maskedSeed[i] = seed[i] ^ seedMask[i];
    std::vector<uint8_t> em;
    em.push_back(0x00);
    em.insert(em.end(), maskedSeed.begin(), maskedSeed.end());
    em.insert(em.end(), maskedDB.begin(), maskedDB.end());
    return em;
}

std::string oaep_decode(const std::vector<uint8_t>& em, size_t k)
{
    if (em.size() != k || k < 2 * kHashLen + 2 || em[0] != 0)
        throw std::runtime_error("OAEP: decrypt error");
    std::vector<uint8_t> maskedSeed(em.begin() + 1, em.begin() + 1 + (std::ptrdiff_t)kHashLen);
    std::vector<uint8_t> maskedDB(em.begin() + 1 + (std::ptrdiff_t)kHashLen, em.end());
    auto seedMask = mgf1(maskedDB.data(), maskedDB.size(), kHashLen);
    std::vector<uint8_t> seed(kHashLen);
    for (size_t i = 0; i < kHashLen; ++i)
        seed[i] = maskedSeed[i] ^ seedMask[i];
    auto dbMask = mgf1(seed.data(), seed.size(), k - kHashLen - 1);
    std::vector<uint8_t> db(maskedDB.size());
    for (size_t i = 0; i < maskedDB.size(); ++i)
        db[i] = maskedDB[i] ^ dbMask[i];
    auto lHash = sha256(nullptr, 0);
    if (!std::equal(lHash.begin(), lHash.end(), db.begin()))
        throw std::runtime_error("OAEP: decrypt error");
    size_t i = kHashLen;
    while (i < db.size() && db[i] == 0)
        ++i;
    if (i >= db.size() || db[i] != 0x01)
        throw std::runtime_error("OAEP: decrypt error");
    ++i;
    return std::string(db.begin() + (std::ptrdiff_t)i, db.end());
}

integer gen_prime(int bits, inf::primality& check)
{
    integer p = inf::prime(bits, check);
    while (log2(p) != bits)
        p = inf::prime(bits, check);
    return p;
}

size_t modulus_bytes(const integer& n)
{
    return (size_t)((log2(n) + 7) / 8);
}

} // namespace

struct PublicKey
{
    integer n, e;

    integer encrypt(const std::string& plaintext) const
    {
        size_t k = modulus_bytes(n);
        auto em = oaep_encode(plaintext, k);
        return modexp(os2ip(em), e, n);
    }
};

class RSA
{
    integer p, q, d, dp, dq, qinv;
public:
    PublicKey pub;

    explicit RSA(int bits = 2048)
    {
        if (bits < 16)
            bits = 16;
        inf::primality check(8);
        const integer ee = 65537;
        const int half = bits / 2;
        for (;;)
        {
            p = gen_prime(half, check);
            q = gen_prime(half, check);
            while (q == p)
                q = gen_prime(half, check);
            integer diff = p > q ? p - q : q - p;
            if (log2(diff) < 900 && bits >= 2048)
                continue;
            integer phi = (p - 1) * (q - 1);
            if (gcd(ee, phi) != 1)
                continue;
            pub.n = p * q;
            pub.e = ee;
            d = modinv(ee, phi);
            dp = d % (p - 1);
            dq = d % (q - 1);
            qinv = modinv(q, p);
            int64_t nbits = log2(pub.n);
            if (bits < 2048 || (nbits >= bits - 1 && nbits <= bits))
                break;
        }
    }

    std::string decrypt(const integer& ciphertext) const
    {
        integer c1 = ciphertext % p;
        integer c2 = ciphertext % q;
        integer m1 = modexp(c1, dp, p);
        integer m2 = modexp(c2, dq, q);
        integer h = (qinv * (m1 - m2)) % p;
        if (h < 0)
            h += p;
        integer m = m2 + h * q;
        size_t k = modulus_bytes(pub.n);
        return oaep_decode(i2osp(m, k), k);
    }
};

struct Message
{
    enum Kind { Announce, Cipher };
    Kind kind;
    std::string from;
    std::string to;
    PublicKey pub;
    integer ciphertext;
};

class Channel
{
    std::map<std::string, std::vector<Message>> boxes_;
public:
    void join(const std::string& name)
    {
        boxes_[name];
    }

    void send(Message m)
    {
        if (m.to.empty())
        {
            for (auto& kv : boxes_)
                kv.second.push_back(m);
            return;
        }
        boxes_[m.to].push_back(m);
    }

    std::vector<Message> recv(const std::string& name)
    {
        std::vector<Message> out;
        std::swap(out, boxes_[name]);
        return out;
    }
};

class Agent
{
    std::string name_;
    RSA keys_;
    std::map<std::string, PublicKey> peers_;
    Channel& ch_;

    void ingest(const Message& m)
    {
        if (m.kind == Message::Announce && m.from != name_)
            peers_[m.from] = m.pub;
    }

public:
    Agent(std::string name, Channel& ch, int bits = 2048)
        : name_(std::move(name)), keys_(bits), ch_(ch)
    {
        ch_.join(name_);
    }

    void announce()
    {
        Message m;
        m.kind = Message::Announce;
        m.from = name_;
        m.pub = keys_.pub;
        ch_.send(m);
    }

    void listen()
    {
        for (const auto& m : ch_.recv(name_))
            ingest(m);
    }

    void send(const std::string& to, const std::string& plaintext)
    {
        auto it = peers_.find(to);
        if (it == peers_.end())
            throw std::runtime_error("unknown peer");
        Message m;
        m.kind = Message::Cipher;
        m.from = name_;
        m.to = to;
        m.ciphertext = it->second.encrypt(plaintext);
        ch_.send(m);
    }

    std::vector<std::string> receive()
    {
        std::vector<std::string> out;
        for (const auto& m : ch_.recv(name_))
        {
            ingest(m);
            if (m.kind == Message::Cipher && m.to == name_)
                out.push_back(keys_.decrypt(m.ciphertext));
        }
        return out;
    }
};

int main()
{
    Channel wire;
    Agent alice("alice", wire);
    Agent bob("bob", wire);
    alice.announce();
    bob.announce();
    alice.listen();
    bob.listen();
    alice.send("bob", "hello from Alice");
    auto inbox = bob.receive();
    assert(inbox.size() == 1 && inbox[0] == "hello from Alice");
    bob.send("alice", "hello from Bob");
    assert(alice.receive().at(0) == "hello from Bob");
    alice.send("bob", "hello from Alice");
    assert(bob.receive().at(0) == "hello from Alice");
    std::cout << "alice <-> bob: ok\n";
    return 0;
}
