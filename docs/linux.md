# Linux build guide

This guide targets Debian/Ubuntu-like distributions but applies to most modern Linux
setups. It explains how to provision compilers, use the scripted build flow, and diagnose
common issues when compiling AlmondFEM.

## Required packages
- Build tools: `sudo apt install build-essential cmake ninja-build git`
- Optional LLVM toolchain: `sudo apt install clang lld`
- vcpkg (optional dependency manager):
  ```bash
  git clone https://github.com/microsoft/vcpkg.git ~/vcpkg
  ~/vcpkg/bootstrap-vcpkg.sh
  export VCPKG_ROOT="$HOME/vcpkg"
  export PATH="$VCPKG_ROOT:$PATH"
  ```
- Vulkan SDK (optional, for graphics experiments): install from [LunarG](https://vulkan.lunarg.com/sdk/home#linux)
  and export `VULKAN_SDK`/`VK_LAYER_PATH`.

Persist any environment variables in your shell profile (`~/.bashrc`, `~/.zshrc`, etc.).

## Scripted workflow
The root-level shell scripts wrap CMake with consistent arguments and directory layouts:

```bash
./cmake/configure.sh gcc Debug
./build.sh gcc Debug
./install.sh gcc Debug
./run.sh gcc Debug
```

Swap `gcc` for `clang` to exercise the LLVM toolchain. Choose `Release` or `RelWithDebInfo`
for optimised binaries. `install.sh` deposits outputs in `built/bin/<Compiler>-<Config>` so
you can run them outside of the build tree.

## IDE integrations
- **VS Code** – install the *C/C++*, *CMake Tools*, and *CMake* extensions. Select the
  desired kit and use the command palette or the scripts above.
- **CLion** – open the repository folder and select the preset (e.g. `gcc-debug`) when
  prompted. CLion reads `CMakePresets.json` automatically.

## Troubleshooting
- **Scripts not executable** – run `chmod +x cmake/*.sh *.sh`.
- **Ninja missing** – install via `sudo apt install ninja-build` or remove it from `PATH` to
  fall back to Unix Makefiles.
- **vcpkg not detected** – ensure `VCPKG_ROOT` or `VCPKG_INSTALLATION_ROOT` is exported before
  configuring.
- **Linker errors referencing Vulkan** – source `setup-env.sh` from the Vulkan SDK or unset
  Vulkan-related environment variables if you are not using the graphics demos.
