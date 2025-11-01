#pragma once

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <iostream>
#include <string>
#include <string_view>

namespace safe_io {
    namespace detail {
        void configure_console() noexcept;
    }

    inline std::ostream& out() noexcept {
        detail::configure_console();
        return std::cout;
    }

    // Format to string
    template <class... Args>
    [[nodiscard]] inline std::string sformat(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        return fmt::format(fmt_str, std::forward<Args>(args)...);
    }

    // Print to stdout
    template <class... Args>
    inline void print(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        fmt::print(stdout, fmt_str, std::forward<Args>(args)...);
        fmt::print(stdout, "\n");
        std::fflush(stdout);
    }

    // Print to stderr
    template <class... Args>
    inline void eprint(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
        fmt::print(stderr, "\n");
        std::fflush(stderr);
    }

    // Print and terminate
    template <class... Args>
    [[noreturn]] inline void fatal(fmt::format_string<Args...> fmt_str, Args&&... args) {
        detail::configure_console();
        fmt::print(stderr, fmt_str, std::forward<Args>(args)...);
        fmt::print(stderr, "\n");
        std::fflush(stderr);
        std::terminate();
    }

} // namespace safe_io
