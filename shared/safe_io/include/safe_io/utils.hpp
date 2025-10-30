#pragma once

#include <fmt/core.h>
#include <iostream>
#include <utility>

namespace safe_io
{
    // Returns a reference to std::cout for safe, global access to output.
    std::ostream& out() noexcept;

    // A wrapper function for formatted output using {fmt}.
    // This function allows type-safe, efficient formatting with automatic newline.
    template <typename... Args>
    void print(fmt::format_string<Args...> fmt_str, Args&&... args)
    {
        auto message = fmt::vformat(fmt_str, fmt::make_format_args(std::forward<Args>(args)...));
        out() << message << '\n';
    }
} // namespace safe_io
