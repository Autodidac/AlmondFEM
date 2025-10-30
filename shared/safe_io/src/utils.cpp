#include "safe_io/utils.hpp"

#include <cstdio>
#include <iostream>

#ifdef _WIN32
#include <Windows.h>
#include <io.h>
#endif

namespace safe_io {
    namespace detail {

        void configure_console() noexcept {
#ifdef _WIN32
            static const bool configured = [] {
                // Only tweak real terminals.
                if (!_isatty(_fileno(stdout)) && !_isatty(_fileno(stderr)))
                    return true;

                // UTF-8 code pages.
                ::SetConsoleOutputCP(CP_UTF8);
                ::SetConsoleCP(CP_UTF8);

                // Enable ANSI/VT sequences for colors.
                const HANDLE hOut = ::GetStdHandle(STD_OUTPUT_HANDLE);
                if (hOut != INVALID_HANDLE_VALUE) {
                    DWORD mode = 0;
                    if (::GetConsoleMode(hOut, &mode)) {
                        ::SetConsoleMode(hOut, mode | ENABLE_PROCESSED_OUTPUT | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
                    }
                }

                // Do NOT switch to _O_U8TEXT. We print narrow UTF-8 via fmt.
                return true;
                }();
            (void)configured;
#endif
        }

    } // namespace detail
} // namespace safe_io
