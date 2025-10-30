#pragma once

#include <fmt/core.h>
#include <iostream>
#include <utility>

namespace safe_io
{
    namespace detail
    {
        // Platform specific setup to keep console output UTF-8 friendly.
        void configure_console() noexcept;
    }

    // Returns a reference to std::cout for safe, global access to output.
    inline std::ostream& out() noexcept
    {
        detail::configure_console();
        return std::cout;
    }

    // A wrapper function for formatted output using {fmt}.
    // This function allows type-safe, efficient formatting with automatic newline.
    template <typename... Args>
    void print(fmt::format_string<Args...> fmt_str, Args&&... args)
    {
        out() << fmt::vformat(fmt_str, fmt::make_format_args(std::forward<Args>(args)...)) << '\n';
    }
} // namespace safe_io
