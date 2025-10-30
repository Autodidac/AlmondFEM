#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace almond::fem::detail
{
    class DenseMatrix
    {
    public:
        DenseMatrix() = default;

        explicit DenseMatrix(std::size_t dimension)
            : m_dimension(dimension)
            , m_data(dimension * dimension, 0.0)
        {
        }

        void resize(std::size_t dimension)
        {
            m_dimension = dimension;
            m_data.assign(dimension * dimension, 0.0);
        }

        [[nodiscard]] std::size_t size() const noexcept { return m_dimension; }

        [[nodiscard]] double& operator()(std::size_t row, std::size_t column) noexcept
        {
            assert(row < m_dimension && column < m_dimension);
            return m_data[row * m_dimension + column];
        }

        [[nodiscard]] double operator()(std::size_t row, std::size_t column) const noexcept
        {
            assert(row < m_dimension && column < m_dimension);
            return m_data[row * m_dimension + column];
        }

        [[nodiscard]] const std::vector<double>& data() const noexcept { return m_data; }
        [[nodiscard]] std::vector<double>& data() noexcept { return m_data; }

    private:
        std::size_t m_dimension{0};
        std::vector<double> m_data{};
    };

    inline std::vector<double> solve(DenseMatrix matrix, std::vector<double> rhs, double pivot_tolerance)
    {
        const auto n = matrix.size();
        std::vector<double> solution(n, 0.0);

        for (std::size_t k = 0; k < n; ++k)
        {
            std::size_t pivot_row = k;
            double pivot_value = std::abs(matrix(k, k));
            for (std::size_t i = k + 1; i < n; ++i)
            {
                const double candidate = std::abs(matrix(i, k));
                if (candidate > pivot_value)
                {
                    pivot_value = candidate;
                    pivot_row = i;
                }
            }

            if (pivot_value <= pivot_tolerance)
            {
                throw std::runtime_error("Matrix is singular to working precision");
            }

            if (pivot_row != k)
            {
                for (std::size_t j = k; j < n; ++j)
                {
                    std::swap(matrix(k, j), matrix(pivot_row, j));
                }
                std::swap(rhs[k], rhs[pivot_row]);
            }

            const double pivot = matrix(k, k);
            for (std::size_t i = k + 1; i < n; ++i)
            {
                const double factor = matrix(i, k) / pivot;
                matrix(i, k) = 0.0;
                for (std::size_t j = k + 1; j < n; ++j)
                {
                    matrix(i, j) -= factor * matrix(k, j);
                }
                rhs[i] -= factor * rhs[k];
            }
        }

        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(n) - 1; i >= 0; --i)
        {
            double sum = rhs[static_cast<std::size_t>(i)];
            for (std::size_t j = static_cast<std::size_t>(i) + 1; j < n; ++j)
            {
                sum -= matrix(static_cast<std::size_t>(i), j) * solution[j];
            }
            solution[static_cast<std::size_t>(i)] = sum / matrix(static_cast<std::size_t>(i), static_cast<std::size_t>(i));
        }

        return solution;
    }
} // namespace almond::fem::detail

