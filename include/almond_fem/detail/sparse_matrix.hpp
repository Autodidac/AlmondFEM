#pragma once

#include "dense_matrix.hpp"

#include <cstddef>
#include <vector>

namespace almond::fem::detail
{
    struct CooMatrix
    {
        std::size_t dimension{0};
        std::vector<int> rows{};
        std::vector<int> cols{};
        std::vector<double> values{};
        std::vector<int> row_prefix{};
    };

    struct CsrMatrix
    {
        std::size_t dimension{0};
        std::vector<int> row_ptr{};
        std::vector<int> col_idx{};
        std::vector<double> values{};
    };

    inline CsrMatrix compress_to_csr(const CooMatrix& coo)
    {
        CsrMatrix csr;
        csr.dimension = coo.dimension;
        csr.row_ptr = coo.row_prefix;
        csr.col_idx = coo.cols;
        csr.values = coo.values;
        return csr;
    }

    inline DenseMatrix csr_to_dense(const CsrMatrix& csr)
    {
        DenseMatrix dense(csr.dimension);
        for (std::size_t row = 0; row < csr.dimension; ++row)
        {
            const auto row_begin = static_cast<std::size_t>(csr.row_ptr[row]);
            const auto row_end = static_cast<std::size_t>(csr.row_ptr[row + 1]);
            for (std::size_t idx = row_begin; idx < row_end; ++idx)
            {
                const auto column = static_cast<std::size_t>(csr.col_idx[idx]);
                dense(row, column) = csr.values[idx];
            }
        }
        return dense;
    }
} // namespace almond::fem::detail
