#include "safe_io/utils.hpp"

#include <cstdio>
#include <iostream>

#ifdef _WIN32
#    include <Windows.h>
#    include <fcntl.h>
#    include <io.h>
#endif

namespace safe_io
{
    namespace
    {
        void configure_console() noexcept
        {
#ifdef _WIN32
            static const bool configured = [] {
                if (!_isatty(_fileno(stdout)))
                {
                    return true;
                }

                ::SetConsoleOutputCP(CP_UTF8);
                ::SetConsoleCP(CP_UTF8);

                const HANDLE handle = ::GetStdHandle(STD_OUTPUT_HANDLE);
                if (handle != INVALID_HANDLE_VALUE)
                {
                    DWORD mode = 0;
                    if (::GetConsoleMode(handle, &mode))
                    {
                        ::SetConsoleMode(handle, mode | ENABLE_PROCESSED_OUTPUT | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
                    }
                }

                _setmode(_fileno(stdout), _O_U8TEXT);
                return true;
            }();
            (void)configured;
#endif
        }
    } // namespace

    std::ostream& out() noexcept
    {
        configure_console();
        return std::cout;
    }
} // namespace safe_io
