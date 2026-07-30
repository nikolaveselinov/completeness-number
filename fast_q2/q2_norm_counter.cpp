#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * Exact binary Gray-code counter for the packed reduced-norm quadratic forms
 * produced by drinfeld_complete.core._quadratic_form_q2.
 *
 * Input (whitespace separated):
 *
 *   dimension max_degree field_degree word_count target_count stop_when_seen
 *   diagonal[dimension][word_count]
 *   cross[i][j][word_count] for 0 <= i < j < dimension
 *   target_polynomial_codes[target_count]
 *
 * Packed words and target codes may be decimal or 0x-prefixed hexadecimal.
 * Output:
 *
 *   SUMMARY iterations exhaustive seen target_count invalid seconds
 *   COUNT target_code multiplicity
 *
 * "exhaustive" is one precisely when all 2^dimension-1 nonzero vectors were
 * visited.  With stop_when_seen=1, a non-exhaustive run is nevertheless a
 * rigorous positivity certificate for every requested target: each COUNT is
 * accompanied by at least one explicitly encountered vector.
 */

namespace {

// Degree-nine characteristics require 20 coefficients over F_{2^18},
// i.e. 360 packed bits.  Eight words cover the degree-nine and degree-ten
// computations while keeping the fixed-size hot-loop representation compact.
constexpr unsigned MAX_WORDS = 8;
constexpr unsigned TABLE_CHUNK_BITS = 8;

struct Packed {
    std::array<std::uint64_t, MAX_WORDS> word{};
};

inline void xor_into(Packed& left, const Packed& right, unsigned words) {
    for (unsigned k = 0; k < words; ++k) {
        left.word[k] ^= right.word[k];
    }
}

inline Packed xored(Packed left, const Packed& right, unsigned words) {
    xor_into(left, right, words);
    return left;
}

std::uint64_t parse_u64(const std::string& token) {
    std::size_t used = 0;
    const auto value = std::stoull(token, &used, 0);
    if (used != token.size()) {
        throw std::runtime_error("invalid integer token: " + token);
    }
    return value;
}

Packed read_packed(unsigned words) {
    Packed value;
    for (unsigned k = 0; k < words; ++k) {
        std::string token;
        if (!(std::cin >> token)) {
            throw std::runtime_error("unexpected end of packed input");
        }
        value.word[k] = parse_u64(token);
    }
    return value;
}

inline std::uint64_t coefficient(
    const Packed& packed,
    unsigned index,
    unsigned field_degree,
    std::uint64_t field_mask
) {
    const unsigned bit = index * field_degree;
    const unsigned word = bit >> 6;
    const unsigned shift = bit & 63U;
    std::uint64_t value = packed.word[word] >> shift;
    if (shift + field_degree > 64U) {
        value |= packed.word[word + 1] << (64U - shift);
    }
    return value & field_mask;
}

/*
 * Normalize a coherent raw norm.  Its coefficients must all be either zero
 * or its leading coefficient.  Return the monic F_2[T] bit code, or UINT32_MAX
 * if the coherence invariant is violated.
 */
inline std::uint32_t normalized_code(
    const Packed& packed,
    unsigned max_degree,
    unsigned field_degree,
    std::uint64_t field_mask
) {
    int degree = static_cast<int>(max_degree);
    std::uint64_t leading = 0;
    while (degree >= 0) {
        leading = coefficient(
            packed, static_cast<unsigned>(degree), field_degree, field_mask
        );
        if (leading != 0) {
            break;
        }
        --degree;
    }
    if (degree < 0) {
        return UINT32_MAX;
    }

    std::uint32_t code = 0;
    for (int k = 0; k <= degree; ++k) {
        const auto value = coefficient(
            packed, static_cast<unsigned>(k), field_degree, field_mask
        );
        if (value == leading) {
            code |= std::uint32_t{1} << k;
        } else if (value != 0) {
            return UINT32_MAX;
        }
    }
    return code;
}

}  // namespace

int main() {
    try {
        unsigned dimension = 0;
        unsigned max_degree = 0;
        unsigned field_degree = 0;
        unsigned words = 0;
        unsigned target_count = 0;
        unsigned stop_when_seen = 0;
        if (!(std::cin >> dimension >> max_degree >> field_degree >> words
              >> target_count >> stop_when_seen)) {
            throw std::runtime_error("missing header");
        }
        if (dimension == 0 || dimension >= 63) {
            throw std::runtime_error("dimension must lie in [1,62]");
        }
        if (max_degree >= 31) {
            throw std::runtime_error("max_degree must be below 31");
        }
        if (field_degree == 0 || field_degree >= 64) {
            throw std::runtime_error("field_degree must lie in [1,63]");
        }
        if (words == 0 || words > MAX_WORDS) {
            throw std::runtime_error("unsupported packed word count");
        }
        if ((max_degree + 1) * field_degree > words * 64) {
            throw std::runtime_error("packed word count is too small");
        }

        std::vector<Packed> diagonal(dimension);
        for (auto& value : diagonal) {
            value = read_packed(words);
        }

        std::vector<std::vector<Packed>> cross(
            dimension, std::vector<Packed>(dimension)
        );
        for (unsigned i = 0; i < dimension; ++i) {
            for (unsigned j = i + 1; j < dimension; ++j) {
                cross[i][j] = read_packed(words);
                cross[j][i] = cross[i][j];
            }
        }

        const unsigned code_space = 1U << (max_degree + 1);
        std::vector<int> target_index(code_space, -1);
        std::vector<std::uint32_t> targets(target_count);
        for (unsigned k = 0; k < target_count; ++k) {
            std::string token;
            if (!(std::cin >> token)) {
                throw std::runtime_error("missing target polynomial");
            }
            const auto code = static_cast<std::uint32_t>(parse_u64(token));
            if (code >= code_space) {
                throw std::runtime_error("target code exceeds max_degree");
            }
            if (target_index[code] != -1) {
                throw std::runtime_error("duplicate target polynomial");
            }
            target_index[code] = static_cast<int>(k);
            targets[k] = code;
        }

        const unsigned chunks =
            (dimension + TABLE_CHUNK_BITS - 1) / TABLE_CHUNK_BITS;
        using Table = std::vector<Packed>;
        std::vector<std::vector<Table>> lookup(
            dimension, std::vector<Table>(chunks)
        );
        for (unsigned flip = 0; flip < dimension; ++flip) {
            for (unsigned chunk = 0; chunk < chunks; ++chunk) {
                const unsigned offset = chunk * TABLE_CHUNK_BITS;
                const unsigned width =
                    std::min(TABLE_CHUNK_BITS, dimension - offset);
                auto& table = lookup[flip][chunk];
                table.resize(1U << width);
                for (unsigned mask = 1; mask < table.size(); ++mask) {
                    const unsigned low = std::countr_zero(mask);
                    table[mask] = table[mask & (mask - 1)];
                    xor_into(
                        table[mask], cross[flip][offset + low], words
                    );
                }
            }
        }

        const std::uint64_t total = (std::uint64_t{1} << dimension) - 1;
        const std::uint64_t field_mask =
            (std::uint64_t{1} << field_degree) - 1;
        std::vector<std::uint64_t> counts(target_count, 0);
        unsigned seen = 0;
        std::uint64_t invalid = 0;
        std::uint64_t iterations = 0;
        std::uint64_t active = 0;
        Packed packed;

        const auto started = std::chrono::steady_clock::now();
        for (std::uint64_t n = 1; n <= total; ++n) {
            const unsigned flip = std::countr_zero(n);
            Packed delta = diagonal[flip];
            for (unsigned chunk = 0; chunk < chunks; ++chunk) {
                const unsigned offset = chunk * TABLE_CHUNK_BITS;
                const unsigned width =
                    std::min(TABLE_CHUNK_BITS, dimension - offset);
                const unsigned mask =
                    static_cast<unsigned>((active >> offset)
                                          & ((1U << width) - 1));
                xor_into(delta, lookup[flip][chunk][mask], words);
            }
            xor_into(packed, delta, words);
            active ^= std::uint64_t{1} << flip;
            ++iterations;

            const auto code = normalized_code(
                packed, max_degree, field_degree, field_mask
            );
            if (code == UINT32_MAX) {
                ++invalid;
                continue;
            }
            const int index = target_index[code];
            if (index >= 0) {
                if (counts[index] == 0) {
                    ++seen;
                }
                ++counts[index];
                if (stop_when_seen && seen == target_count) {
                    break;
                }
            }
        }
        const auto stopped = std::chrono::steady_clock::now();
        const double seconds =
            std::chrono::duration<double>(stopped - started).count();

        /*
         * A nonzero Hom vector may never have zero raw norm in this definite
         * setting, and coherence was verified while constructing the form.
         * Treat any invalid normalization as a hard failure.
         */
        if (invalid != 0) {
            throw std::runtime_error(
                "encountered " + std::to_string(invalid)
                + " incoherent or zero raw norms"
            );
        }

        std::cout << std::setprecision(17);
        std::cout << "SUMMARY " << iterations << ' '
                  << (iterations == total ? 1 : 0) << ' ' << seen << ' '
                  << target_count << ' ' << invalid << ' ' << seconds << '\n';
        for (unsigned k = 0; k < target_count; ++k) {
            std::cout << "COUNT " << targets[k] << ' ' << counts[k] << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "q2_norm_counter: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
