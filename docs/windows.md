# Windows build guide

This guide captures the tooling required to build and debug AlmondFEM on Windows. You can
use the scripted flow (via Git Bash or WSL), Visual Studio presets, or the standalone
solution depending on preference.

## Prerequisites
- **Visual Studio 2022** (or Build Tools) with the *Desktop development with C++* workload.
- **Ninja** (recommended) installed via [Chocolatey](https://chocolatey.org/) or scoop.
- **Git for Windows** to provide Git Bash and Unix-compatible shells for the scripts.
- **vcpkg** (optional dependency manager):
  ```powershell
  git clone https://github.com/microsoft/vcpkg.git C:\dev\vcpkg
  C:\dev\vcpkg\bootstrap-vcpkg.bat
  $env:VCPKG_ROOT = "C:\dev\vcpkg"
  $env:PATH = "$env:VCPKG_ROOT;$env:PATH"
  ```
- **Vulkan SDK** (optional) from [LunarG](https://vulkan.lunarg.com/sdk/home#windows). Add
  `%VULKAN_SDK%\Bin` and `%VULKAN_SDK%\Lib` to `PATH` if using the graphics experiments.

## Scripted workflow
Run the scripts from **x64 Native Tools Command Prompt for VS** or **Developer PowerShell**
so `cl.exe` and Ninja are available:

```powershell
bash ./cmake/configure.sh msvc Debug
bash ./build.sh msvc Debug
bash ./install.sh msvc Debug
bash ./run.sh msvc Debug
```

The scripts accept `msvc`, `gcc`, or `clang` to match your installed toolchains. Outputs are
written to `Bin/<Compiler>-<Config>` and `built/bin/<Compiler>-<Config>`.

## Visual Studio workflows
- Open `app1.sln` for the classic solution experience. Projects link against the same
  sources as the CMake targets.
- Alternatively, choose **File → Open → CMake...** and load `CMakeLists.txt`. Visual Studio
  reads `CMakePresets.json`, exposing the same presets as the scripts.

## VS Code workflows
1. Install the *C/C++*, *CMake Tools*, and *CMake* extensions.
2. Open the repository folder and select the `msvc` kit.
3. Use the *CMake: Configure* and *CMake: Build* commands or call the scripts from the
   integrated terminal.

## Troubleshooting
- **vcpkg not detected** – set `VCPKG_ROOT` or `VCPKG_INSTALLATION_ROOT` before running the
  configure script.
- **Ninja missing** – `choco install ninja` or rely on the bundled Visual Studio generator if
  Ninja is unavailable.
- **Link errors** – ensure you launched the x64 developer prompt so MSVC environment
  variables are initialised.
- **Vulkan loader issues** – confirm `%VULKAN_SDK%` is exported and `%VK_LAYER_PATH%` points to
  the SDK `Bin` directory.
