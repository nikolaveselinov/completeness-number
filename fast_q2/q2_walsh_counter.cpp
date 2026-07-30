#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * Exact Boolean-Walsh counter for an F_2-valued packed quadratic form.
 *
 * This is an isolated experimental executable.  In particular, it shares no
 * build target or mutable state with q2_norm_counter.
 *
 * Input (whitespace separated):
 *
 *   dimension max_degree field_degree word_count include_zero_vector
 *   diagonal[dimension][word_count]
 *   cross[i][j][word_count] for 0 <= i < j < dimension
 *
 * Packed words are little-endian 64-bit limbs and may be decimal or
 * 0x-prefixed hexadecimal.  Every extension-field coefficient block is
 * checked before counting: a block must be literally the integer 0 or 1.
 * This precondition is essential because the method treats the packed norm
 * as a vector-valued Boolean quadratic map.
 *
 * Output:
 *
 *   SUMMARY dimension max_degree character_sums nonzero_sums
 *           include_zero_vector seconds
 *   COUNT polynomial_code multiplicity
 *
 * One COUNT line is emitted for every code from zero through
 * 2^(max_degree+1)-1.
 */

namespace {

using SignedWide = __int128_t;
using UnsignedWide = __uint128_t;

constexpr unsigned MAX_DIMENSION = 62;

struct BooleanForm {
    std::uint64_t linear = 0;
    std::array<std::uint64_t, MAX_DIMENSION> rows{};
};

std::uint64_t parse_u64(const std::string& token) {
    std::size_t used = 0;
    const auto value = std::stoull(token, &used, 0);
    if (used != token.size()) {
        throw std::runtime_error("invalid integer token: " + token);
    }
    return value;
}

std::string decimal(UnsignedWide value) {
    if (value == 0) {
        return "0";
    }
    std::string result;
    while (value != 0) {
        const auto digit = static_cast<unsigned>(value % 10);
        result.push_back(static_cast<char>('0' + digit));
        value /= 10;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

bool bit_at(const std::vector<std::uint64_t>& words, std::size_t bit) {
    return ((words[bit >> 6] >> (bit & 63U)) & 1U) != 0;
}

std::vector<std::uint64_t> read_words(unsigned word_count) {
    std::vector<std::uint64_t> words(word_count);
    for (auto& word : words) {
        std::string token;
        if (!(std::cin >> token)) {
            throw std::runtime_error("unexpected end of packed input");
        }
        word = parse_u64(token);
    }
    return words;
}

template <typename AddMonomial>
void validate_and_transpose(
    const std::vector<std::uint64_t>& words,
    const std::string& label,
    unsigned max_degree,
    unsigned field_degree,
    AddMonomial add_monomial
) {
    const std::size_t packed_width =
        (static_cast<std::size_t>(max_degree) + 1) * field_degree;
    for (unsigned degree = 0; degree <= max_degree; ++degree) {
        const std::size_t start =
            static_cast<std::size_t>(degree) * field_degree;
        const bool constant_one = bit_at(words, start);
        for (unsigned offset = 1; offset < field_degree; ++offset) {
            if (bit_at(words, start + offset)) {
                throw std::runtime_error(
                    label + ", T^" + std::to_string(degree)
                    + " coefficient block is not literally 0 or 1"
                );
            }
        }
        if (constant_one) {
            add_monomial(degree);
        }
    }

    const std::size_t supplied_width =
        static_cast<std::size_t>(words.size()) * 64;
    for (std::size_t bit = packed_width; bit < supplied_width; ++bit) {
        if (bit_at(words, bit)) {
            throw std::runtime_error(
                label + " has nonzero bits beyond the requested packed width"
            );
        }
    }
}

SignedWide quadratic_character_sum(
    const BooleanForm& form,
    unsigned dimension
) {
    auto rows = form.rows;
    std::uint64_t linear = form.linear;
    std::uint64_t active =
        (std::uint64_t{1} << dimension) - std::uint64_t{1};
    bool constant = false;
    unsigned eliminated_pairs = 0;

    while (true) {
        unsigned left = dimension;
        unsigned right = dimension;
        std::uint64_t remaining = active;
        while (remaining != 0) {
            const unsigned candidate = std::countr_zero(remaining);
            const std::uint64_t neighbours = rows[candidate] & active;
            if (neighbours != 0) {
                left = candidate;
                right = std::countr_zero(neighbours);
                break;
            }
            remaining &= remaining - 1;
        }
        if (left == dimension) {
            break;
        }

        const std::uint64_t left_bit = std::uint64_t{1} << left;
        const std::uint64_t right_bit = std::uint64_t{1} << right;
        const bool left_constant = (linear & left_bit) != 0;
        const bool right_constant = (linear & right_bit) != 0;
        const std::uint64_t left_terms =
            rows[left] & active & ~right_bit;
        const std::uint64_t right_terms =
            rows[right] & active & ~left_bit;

        active &= ~(left_bit | right_bit);
        linear &= active;
        remaining = active;
        while (remaining != 0) {
            const unsigned index = std::countr_zero(remaining);
            rows[index] &= active;
            remaining &= remaining - 1;
        }

        constant ^= left_constant && right_constant;
        if (left_constant) {
            linear ^= right_terms;
        }
        if (right_constant) {
            linear ^= left_terms;
        }

        std::uint64_t left_remaining = left_terms;
        while (left_remaining != 0) {
            const unsigned left_index = std::countr_zero(left_remaining);
            const std::uint64_t left_index_bit =
                std::uint64_t{1} << left_index;
            std::uint64_t right_remaining = right_terms;
            while (right_remaining != 0) {
                const unsigned right_index =
                    std::countr_zero(right_remaining);
                const std::uint64_t right_index_bit =
                    std::uint64_t{1} << right_index;
                if (left_index == right_index) {
                    linear ^= left_index_bit;
                } else {
                    rows[left_index] ^= right_index_bit;
                    rows[right_index] ^= left_index_bit;
                }
                right_remaining &= right_remaining - 1;
            }
            left_remaining &= left_remaining - 1;
        }
        ++eliminated_pairs;
    }

    if ((linear & active) != 0) {
        return 0;
    }
    const unsigned exponent =
        eliminated_pairs + std::popcount(active);
    const SignedWide magnitude = SignedWide{1} << exponent;
    return constant ? -magnitude : magnitude;
}

void fwht(std::vector<SignedWide>& values) {
    for (std::size_t width = 1; width < values.size(); width *= 2) {
        const std::size_t step = 2 * width;
        for (std::size_t start = 0; start < values.size(); start += step) {
            for (std::size_t offset = 0; offset < width; ++offset) {
                const std::size_t left = start + offset;
                const std::size_t right = left + width;
                const SignedWide left_value = values[left];
                const SignedWide right_value = values[right];
                values[left] = left_value + right_value;
                values[right] = left_value - right_value;
            }
        }
    }
}

}  // namespace

int main() {
    try {
        unsigned dimension = 0;
        unsigned max_degree = 0;
        unsigned field_degree = 0;
        unsigned word_count = 0;
        unsigned include_zero = 0;
        if (!(std::cin >> dimension >> max_degree >> field_degree
              >> word_count >> include_zero)) {
            throw std::runtime_error("missing header");
        }
        if (dimension == 0 || dimension > MAX_DIMENSION) {
            throw std::runtime_error("dimension must lie in [1,62]");
        }
        if (max_degree >= 62) {
            throw std::runtime_error("max_degree must be below 62");
        }
        if (field_degree == 0) {
            throw std::runtime_error("field_degree must be positive");
        }
        if (include_zero > 1) {
            throw std::runtime_error("include_zero_vector must be 0 or 1");
        }
        const std::size_t output_dimension =
            static_cast<std::size_t>(max_degree) + 1;
        if (dimension + output_dimension >= 127) {
            throw std::runtime_error(
                "dimension plus output dimension exceeds exact accumulator"
            );
        }
        if (output_dimension >= std::numeric_limits<std::size_t>::digits) {
            throw std::runtime_error("output code space does not fit size_t");
        }
        if (output_dimension >
            std::numeric_limits<std::size_t>::max() / field_degree) {
            throw std::runtime_error("packed width overflow");
        }
        const std::size_t packed_width = output_dimension * field_degree;
        const std::size_t expected_words = (packed_width + 63) / 64;
        if (word_count == 0 || word_count != expected_words) {
            throw std::runtime_error(
                "word_count must exactly cover the requested packed width"
            );
        }

        std::vector<BooleanForm> components(output_dimension);
        for (unsigned index = 0; index < dimension; ++index) {
            const auto words = read_words(word_count);
            validate_and_transpose(
                words,
                "diagonal[" + std::to_string(index) + "]",
                max_degree,
                field_degree,
                [&](unsigned degree) {
                    components[degree].linear |=
                        std::uint64_t{1} << index;
                }
            );
        }
        for (unsigned left = 0; left < dimension; ++left) {
            for (unsigned right = left + 1; right < dimension; ++right) {
                const auto words = read_words(word_count);
                validate_and_transpose(
                    words,
                    "cross[" + std::to_string(left) + "]["
                        + std::to_string(right) + "]",
                    max_degree,
                    field_degree,
                    [&](unsigned degree) {
                        components[degree].rows[left] |=
                            std::uint64_t{1} << right;
                        components[degree].rows[right] |=
                            std::uint64_t{1} << left;
                    }
                );
            }
        }
        std::string trailing;
        if (std::cin >> trailing) {
            throw std::runtime_error(
                "unexpected trailing input token: " + trailing
            );
        }

        const auto started = std::chrono::steady_clock::now();
        const std::size_t output_count =
            std::size_t{1} << output_dimension;
        std::vector<SignedWide> spectrum(output_count);
        spectrum[0] = SignedWide{1} << dimension;

        BooleanForm scalar_form;
        std::size_t nonzero_character_sums = 1;
        for (std::size_t ordinal = 1; ordinal < output_count; ++ordinal) {
            const unsigned flip = std::countr_zero(ordinal);
            scalar_form.linear ^= components[flip].linear;
            for (unsigned row = 0; row < dimension; ++row) {
                scalar_form.rows[row] ^= components[flip].rows[row];
            }
            const SignedWide character =
                quadratic_character_sum(scalar_form, dimension);
            const std::size_t gray_code = ordinal ^ (ordinal >> 1);
            spectrum[gray_code] = character;
            nonzero_character_sums += character != 0;
        }

        fwht(spectrum);
        const SignedWide divisor =
            static_cast<SignedWide>(output_count);
        UnsignedWide total = 0;
        for (std::size_t code = 0; code < output_count; ++code) {
            if (spectrum[code] % divisor != 0) {
                throw std::runtime_error(
                    "Fourier inversion was nonintegral at code "
                    + std::to_string(code)
                );
            }
            spectrum[code] /= divisor;
            if (spectrum[code] < 0) {
                throw std::runtime_error(
                    "Fourier inversion was negative at code "
                    + std::to_string(code)
                );
            }
        }
        if (!include_zero) {
            if (spectrum[0] < 1) {
                throw std::runtime_error(
                    "zero input was missing from the zero fibre"
                );
            }
            --spectrum[0];
        }
        for (const SignedWide count : spectrum) {
            total += static_cast<UnsignedWide>(count);
        }
        const UnsignedWide expected_total =
            (UnsignedWide{1} << dimension) - (include_zero ? 0 : 1);
        if (total != expected_total) {
            throw std::runtime_error(
                "Walsh multiplicities have the wrong total"
            );
        }

        const double seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - started
            ).count();
        std::cout << "SUMMARY " << dimension << ' ' << max_degree << ' '
                  << output_count << ' ' << nonzero_character_sums << ' '
                  << include_zero << ' ' << std::setprecision(17)
                  << seconds << '\n';
        for (std::size_t code = 0; code < output_count; ++code) {
            std::cout << "COUNT " << code << ' '
                      << decimal(static_cast<UnsignedWide>(spectrum[code]))
                      << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 2;
    }
}
