# Algorithms

inf.cpp names each kernel after the algorithm it implements. They are first-class functions, not hidden internals.

## Operation → algorithm

| Operation | Algorithm | Notes |
|---|---|---|
| `+` / `-` | Schoolbook carry | Optimal. |
| `*` (native int) | `mul` — single-limb schoolbook | Falls back to big×big if the scalar does not fit in 32 bits. |
| `*` (big×big), ≤ 32 limbs | `schoolbook` | |
| `*` (big×big), ≤ 256 limbs | `karatsuba` | Split at half, three recursive products. |
| `*` (big×big), ≤ 4096 limbs | `ntt` | 16-bit digits, length-`2^k` NTT, three 31-bit primes, CRT. |
| `*` (big×big), larger or NTT too long | `schonhage_strassen` | Fermat-ring FFT, \(\mathbb{Z}/(2^N+1)\mathbb{Z}\). |
| `/` `%` (native int) | `div` / `mod` | Single-limb long division. |
| `/` `%` (big×big) | `AlgoD` | Knuth TAOCP Vol. 2 §4.3.1, Algorithm D, base 2³². |
| `/` `%` (both ≥ 64 limbs) | Newton reciprocal | Then one multiply and a short correction; falls back to AlgoD if the estimate is off. |
| `<<` `>>` | Limb + intra-limb bit shift | |
| `gcd` | Lehmer | Native leading digits + 2×2 matrix; Euclidean on a single limb. |
| `lcm` | Via `gcd` | `abs(a)/gcd * abs(b)`. |
| `modinv` | Extended Euclidean | Throws if not invertible. |
| `modexp` | Montgomery + binary exp | Odd moduli. Even moduli use `(r*b)%m`. |
| `sqrt` | Newton | Initial guess from bit length. |
| `pow` | Binary exponentiation | Negative exponents throw. |
| Primality | Baillie–PSW | Small primes, then Miller–Rabin base 2 + strong Lucas. |
| `prime` | Random odd walk | Candidates from XorShift64. |
| `rand` | XorShift64* | Fills 2³² limbs and sets the high bit so the bit length is exact. |
| `sieve` | Segmented trial | Fast path when the range fits in 64-bit. |

## Multiply cutoffs

```
limbs ≤ 32     → schoolbook
limbs ≤ 256    → karatsuba
limbs ≤ 4096   → ntt          (if transform length ≤ 2^23)
otherwise      → schonhage_strassen
```

`ntt` is a Cooley–Tukey **number-theoretic transform**. It splits each 32-bit limb into two 16-bit digits and transforms under

- 998244353 (g = 3)
- 897581057 (g = 3)
- 754974721 (g = 11)

If the length would exceed \(2^{23}\), it calls `schonhage_strassen`.

`schonhage_strassen` implements the practical Fermat-ring FFT (Schönhage & Strassen, *Schnelle Multiplikation großer Zahlen*, Computing 1971): split into \(M=2^m\) pieces, DFT over \(\mathbb{Z}/(2^N+1)\mathbb{Z}\) with twiddles that are shifts, pointwise multiply (via `operator*`, so smaller kernels recurse), invert, reconstruct.

The RSA **test** (`tests/rsa.cpp`) uses RSAES-OAEP with SHA-256 and CRT decryption (RFC 8017 PKCS #1 v2.2 §7.1 / B.2.1; SHA-256 from FIPS 180-4). Those helpers were implemented from the specs, not copied from another library. The library itself only supplies `modexp`, `prime`, and `modinv`.

## Representation

Little-endian `uint32_t` limbs, sign +1 / −1. Decimal appears only at I/O (chunks of 10⁹). Intermediates are `uint64_t`; there is no `__int128` requirement.
