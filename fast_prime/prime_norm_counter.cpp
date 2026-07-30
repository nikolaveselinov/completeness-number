#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * Exact reflected p-ary Gray counter for a quadratic form
 *
 *   Q(x) = sum_i D_i x_i^2 + sum_{i<j} C_ij x_i x_j
 *
 * with coefficients in F_p[T]. Input is whitespace separated:
 *
 *   p dimension max_degree target_count
 *   D_i[k]                              (all i, then all k)
 *   C_ij[k]                             (i<j, then all k)
 *   target polynomial codes
 *
 * A polynomial sum a_k T^k is encoded as sum a_k p^k. Output is:
 *
 *   SUMMARY iterations target_count invalid seconds
 *   COUNT target_code multiplicity
 *
 * Every one of the p^dimension-1 nonzero vectors is visited. Consequently a
 * zero COUNT is an exact absence certificate, not a failed witness search.
 */

namespace {

std::uint64_t checked_power(std::uint64_t base, unsigned exponent) {
    std::uint64_t value = 1;
    for (unsigned k = 0; k < exponent; ++k) {
        if (value > std::numeric_limits<std::uint64_t>::max() / base) {
            throw std::runtime_error("vector space exceeds uint64 range");
        }
        value *= base;
    }
    return value;
}

unsigned inverse_mod(unsigned value, unsigned prime) {
    for (unsigned candidate = 1; candidate < prime; ++candidate) {
        if ((value * candidate) % prime == 1) {
            return candidate;
        }
    }
    throw std::runtime_error("leading coefficient is not invertible");
}

inline void add_scaled(
    std::vector<unsigned>& destination,
    const std::vector<unsigned>& source,
    int scalar,
    unsigned prime
) {
    for (std::size_t k = 0; k < destination.size(); ++k) {
        int value = static_cast<int>(destination[k])
                    + scalar * static_cast<int>(source[k]);
        value %= static_cast<int>(prime);
        if (value < 0) {
            value += static_cast<int>(prime);
        }
        destination[k] = static_cast<unsigned>(value);
    }
}

std::uint64_t normalized_code(
    const std::vector<unsigned>& polynomial,
    unsigned prime
) {
    int degree = static_cast<int>(polynomial.size()) - 1;
    while (degree >= 0 && polynomial[static_cast<std::size_t>(degree)] == 0) {
        --degree;
    }
    if (degree < 0) {
        return std::numeric_limits<std::uint64_t>::max();
    }
    const unsigned inverse = inverse_mod(
        polynomial[static_cast<std::size_t>(degree)], prime
    );
    std::uint64_t code = 0;
    std::uint64_t place = 1;
    for (int k = 0; k <= degree; ++k) {
        const unsigned coefficient =
            (polynomial[static_cast<std::size_t>(k)] * inverse) % prime;
        code += place * coefficient;
        place *= prime;
    }
    return code;
}

std::vector<unsigned> read_polynomial(unsigned length, unsigned prime) {
    std::vector<unsigned> result(length);
    for (auto& coefficient : result) {
        if (!(std::cin >> coefficient)) {
            throw std::runtime_error("unexpected end of quadratic-form input");
        }
        if (coefficient >= prime) {
            throw std::runtime_error("quadratic-form coefficient outside F_p");
        }
    }
    return result;
}

}  // namespace

int main() {
    try {
        unsigned prime = 0;
        unsigned dimension = 0;
        unsigned max_degree = 0;
        unsigned target_count = 0;
        if (!(std::cin >> prime >> dimension >> max_degree >> target_count)) {
            throw std::runtime_error("missing header");
        }
        if (prime < 3 || (prime & 1U) == 0) {
            throw std::runtime_error("p must be an odd prime");
        }
        if (dimension == 0) {
            throw std::runtime_error("dimension must be positive");
        }
        const unsigned length = max_degree + 1;
        std::vector<std::vector<unsigned>> diagonal(dimension);
        for (auto& value : diagonal) {
            value = read_polynomial(length, prime);
        }
        std::vector<std::vector<std::vector<unsigned>>> cross(
            dimension,
            std::vector<std::vector<unsigned>>(
                dimension, std::vector<unsigned>(length, 0)
            )
        );
        for (unsigned i = 0; i < dimension; ++i) {
            for (unsigned j = i + 1; j < dimension; ++j) {
                cross[i][j] = read_polynomial(length, prime);
                cross[j][i] = cross[i][j];
            }
        }

        const std::uint64_t code_space = checked_power(prime, length);
        if (code_space > static_cast<std::uint64_t>(
                             std::numeric_limits<std::size_t>::max()
                         )) {
            throw std::runtime_error("target code space is too large");
        }
        std::vector<int> target_index(
            static_cast<std::size_t>(code_space), -1
        );
        std::vector<std::uint64_t> targets(target_count);
        for (unsigned index = 0; index < target_count; ++index) {
            std::uint64_t code = 0;
            if (!(std::cin >> code)) {
                throw std::runtime_error("missing target polynomial code");
            }
            if (code >= code_space) {
                throw std::runtime_error("target code exceeds max degree");
            }
            if (target_index[static_cast<std::size_t>(code)] != -1) {
                throw std::runtime_error("duplicate target polynomial");
            }
            target_index[static_cast<std::size_t>(code)] =
                static_cast<int>(index);
            targets[index] = code;
        }

        const std::uint64_t total = checked_power(prime, dimension) - 1;
        std::vector<std::uint64_t> counts(target_count, 0);
        std::vector<unsigned> coefficients(dimension, 0);
        std::vector<int> directions(dimension, 1);
        std::vector<unsigned> value(length, 0);
        std::vector<std::vector<unsigned>> gradients(
            dimension, std::vector<unsigned>(length, 0)
        );
        std::uint64_t invalid = 0;

        const auto started = std::chrono::steady_clock::now();
        for (std::uint64_t iteration = 0; iteration < total; ++iteration) {
            unsigned flip = 0;
            while (
                static_cast<int>(coefficients[flip]) + directions[flip] < 0
                || static_cast<int>(coefficients[flip]) + directions[flip]
                       >= static_cast<int>(prime)
            ) {
                directions[flip] = -directions[flip];
                ++flip;
                if (flip >= dimension) {
                    throw std::runtime_error("p-ary Gray traversal overflow");
                }
            }
            const int direction = directions[flip];

            add_scaled(value, gradients[flip], direction, prime);
            add_scaled(value, diagonal[flip], 1, prime);
            for (unsigned index = 0; index < dimension; ++index) {
                if (index != flip) {
                    add_scaled(
                        gradients[index],
                        cross[index][flip],
                        direction,
                        prime
                    );
                }
            }
            add_scaled(
                gradients[flip], diagonal[flip], 2 * direction, prime
            );
            coefficients[flip] = static_cast<unsigned>(
                static_cast<int>(coefficients[flip]) + direction
            );

            const std::uint64_t code = normalized_code(value, prime);
            if (code == std::numeric_limits<std::uint64_t>::max()) {
                ++invalid;
                continue;
            }
            const int index =
                target_index[static_cast<std::size_t>(code)];
            if (index >= 0) {
                ++counts[static_cast<std::size_t>(index)];
            }
        }
        const auto stopped = std::chrono::steady_clock::now();
        const double seconds =
            std::chrono::duration<double>(stopped - started).count();
        if (invalid != 0) {
            throw std::runtime_error(
                "encountered " + std::to_string(invalid)
                + " zero raw norms"
            );
        }

        std::cout << std::setprecision(17);
        std::cout << "SUMMARY " << total << ' ' << target_count << ' '
                  << invalid << ' ' << seconds << '\n';
        for (unsigned index = 0; index < target_count; ++index) {
            std::cout << "COUNT " << targets[index] << ' '
                      << counts[index] << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "prime_norm_counter: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
