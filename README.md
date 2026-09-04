# inf.cpp

infinite integers in one C++ header with zero dependencies.

--- 

A high-efficiency, high-precision infinite-integer arithmetic library. Drop `inf.hpp` into any C++17 project — no GMP, no Boost, no build step.

```cpp
#include "inf.hpp"

int main() {
    inf::integer a = "12345678901234567890";
    inf::integer b = 42;
    std::cout << a * b << "\n";
}
```

## Algorithms

Schoolbook → Karatsuba → NTT → Schönhage–Strassen, Knuth Algorithm D, Lehmer GCD, Montgomery `modexp`, Baillie–PSW. Details: [docs/algorithms.md](docs/algorithms.md). API: [docs/api.md](docs/api.md).

## Layout

| Path | Role |
|---|---|
| [`inf.hpp`](inf.hpp) | The library |
| [`docs/`](docs/) | Algorithms and API |
| [`tests/`](tests/) | `./tests/build` (Alpine image) or `make -C tests unit` / `make -C tests rsa` (RSA-OAEP Alice/Bob) |

## License

MIT — see [LICENSE](LICENSE).
