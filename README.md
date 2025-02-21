# Welcome to mdot

**mdot** is a simple C++ project for matrix multiplication, real, complex,
normal, triangular or hermitian.

You can start the project by using *gcc*
```bash
# entering the shell
nix-shell --arg clangSupport false
# building the project
nix-build . --arg clangSupport false
```
The project is built by default with *clang*.


Tip: one line code formatter for C/C++ projects

```bash
nixpkgs-fmt .

clang-format -i $(find . -path "./b*" -prune  -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")

cmake-format -i $(find . -path "./b*" -prune  -name "*.cmake" -o -name "CMakeLists.txt")
```
