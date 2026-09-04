# API

C++17, one header, standard library only.

```cpp
#include "inf.hpp"

int main() {
    inf::integer a = "12345678901234567890";
    inf::integer b = 42;
    std::cout << a * b << "\n";
    std::cout << inf::pow(a, 3) << "\n";
}
```

The type lives in `namespace inf`. Named algorithms (`schoolbook`, `karatsuba`, `ntt`, `schonhage_strassen`, `AlgoD`, `gcd`, …) are also in that namespace and are found by ADL.

## `inf::integer`

- Constructors: default (0), `const char*`, `std::string`, any integral type, copy/move
- Assignment from the same
- Arithmetic: `+ - * / %` and the compound forms, prefix/postfix `++` `--`, unary `-`
- Shifts: `<< >> <<= >>=` with an integral or `integer` count
- Comparisons: `< > <= >= == !=` against `integer` or integrals
- Stream `<<` / `>>`
- Cast to an arithmetic type (may overflow)
- `c_str()` — decimal `std::string`
- `integer::rand(bits)` — uniformly random **positive** integer with that many bits

## Named algorithms

| Function | Role |
|---|---|
| `add` `sub` `mul` | Building-block arithmetic (`mul` is × a 64-bit scalar) |
| `schoolbook` `karatsuba` `ntt` `schonhage_strassen` | Explicit multiply (`*` dispatches by size) |
| `AlgoD` | Returns `qr<integer>` `{q, r}` |
| `div` `mod` | Divide / remainder by a 64-bit scalar |
| `abs` `min` `max` | |
| `gcd` `lcm` `sqrt` `pow` | |
| `modinv` `modexp` | Modular inverse and exponentiation |
| `sieve(lo, hi)` | Primes in `[lo, hi]` |
| `log2` | Exact bit length |
| `log10` / `size` | Decimal digit count |
| `debug` | Prints sign and raw limbs |

## Primality

```cpp
inf::primality check(8);          // argument kept for API kinship; BPSW does not use it
bool p = check(n);
inf::integer q = inf::prime(64, check);  // 64-bit prime
```

## Hash

`std::hash<inf::integer>` is specialized (decimal string), so `integer` can be an `unordered_map` key.
