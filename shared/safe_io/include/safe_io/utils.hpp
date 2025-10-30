#pragma once

#include <fmt/format.h>   // println, format_string, format
#include <fmt/ostream.h>  // enable fmt for types with operator<<
#include <iostream>       // std::cout
#include <string>
#include <string_view>

namespace safe_io {
    namespace detail {
        void configure_console() noexcept;
    }

    // Optional: keep a stable ostream handle if other code uses it.
    inline std::ostream& out() noexcept {
        detail::configure_console();
        return std::cout;
    }

    // Format → std::string.
    template <class... Args>
    [[nodiscard]] inline std::string sformat(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        return fmt::format(fmt_str, std::forward<Args>(args)...);
    }

    // Print to stdout with newline. Perfect-forward args. No arg-store, no vformat.
    template <class... Args>
    inline void print(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        fmt::println(stdout, fmt_str, std::forward<Args>(args)...);
        std::fflush(stdout);
    }

    // Print to stderr with newline.
    template <class... Args>
    inline void eprint(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        fmt::println(stderr, fmt_str, std::forward<Args>(args)...);
        std::fflush(stderr);
    }

    // Print to stderr and abort.
    template <class... Args>
    [[noreturn]] inline void fatal(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        fmt::println(stderr, fmt_str, std::forward<Args>(args)...);
        std::fflush(stderr);
        std::terminate();
    }

} // namespace safe_io

// If you use fmt::join / ranges anywhere, include this in that TU:
//   #include <fmt/ranges.h>
