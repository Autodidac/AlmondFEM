#pragma once

#include "dense_matrix.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace almond::fem::detail
{
    // Stores explicit triplets produced during assembly.
    struct CooMatrix
    {
        std::size_t dimension{0};
        std::vector<int> rows{};
        std::vector<int> cols{};
        std::vector<double> values{};
        std::vector<int> row_prefix{};
    };

    // Canonical sparse matrix for CPU solves. Owns CSR buffers and provides
    // utility helpers such as dense materialisation for debugging.
    class CsrMatrix
    {
    public:
        CsrMatrix() = default;

        explicit CsrMatrix(const CooMatrix& coo)
            : m_dimension(coo.dimension)
            , m_row_ptr(coo.row_prefix)
            , m_col_idx(coo.cols)
            , m_values(coo.values)
        {
            validate();
        }

        explicit CsrMatrix(CooMatrix&& coo)
            : m_dimension(coo.dimension)
            , m_row_ptr(std::move(coo.row_prefix))
            , m_col_idx(std::move(coo.cols))
            , m_values(std::move(coo.values))
        {
            validate();
        }

        [[nodiscard]] std::size_t dimension() const noexcept { return m_dimension; }
        [[nodiscard]] const std::vector<int>& row_ptr() const noexcept { return m_row_ptr; }
        [[nodiscard]] const std::vector<int>& col_idx() const noexcept { return m_col_idx; }
        [[nodiscard]] const std::vector<double>& values() const noexcept { return m_values; }

        [[nodiscard]] std::vector<int>& row_ptr() noexcept { return m_row_ptr; }
        [[nodiscard]] std::vector<int>& col_idx() noexcept { return m_col_idx; }
        [[nodiscard]] std::vector<double>& values() noexcept { return m_values; }

        [[nodiscard]] DenseMatrix to_dense() const
        {
            DenseMatrix dense(m_dimension);
            for (std::size_t row = 0; row < m_dimension; ++row)
            {
                const auto row_begin = static_cast<std::size_t>(m_row_ptr[row]);
                const auto row_end = static_cast<std::size_t>(m_row_ptr[row + 1]);
                for (std::size_t idx = row_begin; idx < row_end; ++idx)
                {
                    const auto column = static_cast<std::size_t>(m_col_idx[idx]);
                    dense(row, column) = m_values[idx];
                }
            }
            return dense;
        }

    private:
        void validate()
        {
            if (m_row_ptr.size() != m_dimension + 1)
            {
                throw std::invalid_argument("CSR row pointer length must equal dimension + 1");
            }
            if (m_col_idx.size() != m_values.size())
            {
                throw std::invalid_argument("CSR column and value arrays must have identical sizes");
            }
        }

        std::size_t m_dimension{0};
        std::vector<int> m_row_ptr{};
        std::vector<int> m_col_idx{};
        std::vector<double> m_values{};
    };

    // SELL-C-σ layout derived from CSR using locally sorted slices. Enables
    // vector-friendly traversal without mutating the canonical CSR storage.
    class SellCSigmaMatrix
    {
    public:
        SellCSigmaMatrix() = default;

        SellCSigmaMatrix(const CsrMatrix& csr, std::size_t chunk_size = 32)
        {
            build_from_csr(csr, chunk_size);
        }

        [[nodiscard]] std::size_t dimension() const noexcept { return m_dimension; }
        [[nodiscard]] std::size_t chunk_size() const noexcept { return m_chunk_size; }
        [[nodiscard]] const std::vector<int>& slice_ptr() const noexcept { return m_slice_ptr; }
        [[nodiscard]] const std::vector<int>& col_idx() const noexcept { return m_col_idx; }
        [[nodiscard]] const std::vector<double>& values() const noexcept { return m_values; }
        [[nodiscard]] const std::vector<int>& row_map() const noexcept { return m_row_map; }

    private:
        void build_from_csr(const CsrMatrix& csr, std::size_t chunk_size)
        {
            m_dimension = csr.dimension();
            m_chunk_size = std::max<std::size_t>(1, chunk_size);

            const auto& row_ptr = csr.row_ptr();
            const auto& col_idx = csr.col_idx();
            const auto& values = csr.values();

            const std::size_t slice_count = (m_dimension + m_chunk_size - 1) / m_chunk_size;
            m_slice_ptr.clear();
            m_slice_ptr.reserve(slice_count + 1);
            m_slice_ptr.push_back(0);

            m_row_map.clear();
            m_row_map.reserve(m_dimension);

            for (std::size_t slice = 0; slice < slice_count; ++slice)
            {
                const std::size_t row_begin = slice * m_chunk_size;
                const std::size_t row_end = std::min(m_dimension, row_begin + m_chunk_size);

                std::vector<std::pair<std::size_t, std::size_t>> rows;
                rows.reserve(row_end - row_begin);
                for (std::size_t row = row_begin; row < row_end; ++row)
                {
                    const auto start = static_cast<std::size_t>(row_ptr[row]);
                    const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
                    rows.emplace_back(row, end - start);
                }

                std::sort(rows.begin(), rows.end(), [](const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });

                const std::size_t max_nnz = rows.empty() ? 0 : rows.front().second;
                for (const auto& [row_index, nnz] : rows)
                {
                    (void)nnz;
                    m_row_map.push_back(static_cast<int>(row_index));
                }

                for (std::size_t entry = 0; entry < max_nnz; ++entry)
                {
                    for (const auto& [row_index, nnz] : rows)
                    {
                        const auto start = static_cast<std::size_t>(row_ptr[row_index]);
                        if (entry < nnz)
                        {
                            const auto idx = start + entry;
                            m_col_idx.push_back(col_idx[idx]);
                            m_values.push_back(values[idx]);
                        }
                        else
                        {
                            m_col_idx.push_back(-1);
                            m_values.push_back(0.0);
                        }
                    }
                }

                m_slice_ptr.push_back(static_cast<int>(m_col_idx.size()));
            }
        }

        std::size_t m_dimension{0};
        std::size_t m_chunk_size{32};
        std::vector<int> m_slice_ptr{};
        std::vector<int> m_col_idx{};
        std::vector<double> m_values{};
        std::vector<int> m_row_map{};
    };

    // Lightweight ELLPACK view (primarily for documentation/tests) created
    // from CSR data without reassembling element contributions.
    class EllMatrix
    {
    public:
        EllMatrix() = default;

        explicit EllMatrix(const CsrMatrix& csr)
        {
            build_from_csr(csr);
        }

        [[nodiscard]] std::size_t dimension() const noexcept { return m_dimension; }
        [[nodiscard]] std::size_t entries_per_row() const noexcept { return m_entries_per_row; }
        [[nodiscard]] const std::vector<int>& col_idx() const noexcept { return m_col_idx; }
        [[nodiscard]] const std::vector<double>& values() const noexcept { return m_values; }

    private:
        void build_from_csr(const CsrMatrix& csr)
        {
            m_dimension = csr.dimension();
            const auto& row_ptr = csr.row_ptr();
            const auto& col_idx = csr.col_idx();
            const auto& values = csr.values();

            m_entries_per_row = 0;
            for (std::size_t row = 0; row < m_dimension; ++row)
            {
                const auto start = static_cast<std::size_t>(row_ptr[row]);
                const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
                m_entries_per_row = std::max<std::size_t>(m_entries_per_row, end - start);
            }

            m_col_idx.assign(m_dimension * m_entries_per_row, -1);
            m_values.assign(m_dimension * m_entries_per_row, 0.0);

            for (std::size_t row = 0; row < m_dimension; ++row)
            {
                const auto start = static_cast<std::size_t>(row_ptr[row]);
                const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
                for (std::size_t entry = 0; entry < end - start; ++entry)
                {
                    const auto idx = start + entry;
                    const auto offset = row * m_entries_per_row + entry;
                    m_col_idx[offset] = col_idx[idx];
                    m_values[offset] = values[idx];
                }
            }
        }

        std::size_t m_dimension{0};
        std::size_t m_entries_per_row{0};
        std::vector<int> m_col_idx{};
        std::vector<double> m_values{};
    };

    inline void multiply(const CsrMatrix& csr, const std::vector<double>& x, std::vector<double>& result)
    {
        const auto n = csr.dimension();
        if (result.size() != n)
        {
            result.assign(n, 0.0);
        }
        else
        {
            std::fill(result.begin(), result.end(), 0.0);
        }

        const auto& row_ptr = csr.row_ptr();
        const auto& col_idx = csr.col_idx();
        const auto& values = csr.values();

        for (std::size_t row = 0; row < n; ++row)
        {
            const auto begin = static_cast<std::size_t>(row_ptr[row]);
            const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
            double sum = 0.0;
            for (std::size_t idx = begin; idx < end; ++idx)
            {
                const auto column = static_cast<std::size_t>(col_idx[idx]);
                sum += values[idx] * x[column];
            }
            result[row] = sum;
        }
    }

    inline std::vector<double> diagonal(const CsrMatrix& csr)
    {
        const auto n = csr.dimension();
        std::vector<double> diag(n, 0.0);

        const auto& row_ptr = csr.row_ptr();
        const auto& col_idx = csr.col_idx();
        const auto& values = csr.values();

        for (std::size_t row = 0; row < n; ++row)
        {
            const auto begin = static_cast<std::size_t>(row_ptr[row]);
            const auto end = static_cast<std::size_t>(row_ptr[row + 1]);
            for (std::size_t idx = begin; idx < end; ++idx)
            {
                if (static_cast<std::size_t>(col_idx[idx]) == row)
                {
                    diag[row] = values[idx];
                    break;
                }
            }
        }

        return diag;
    }
} // namespace almond::fem::detail
