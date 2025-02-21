# golden_cpp

This is a skeleton template for a C/C++ project.

## My development tools are
- nix :heart: (packaging from hell)
- clang-format (formatter)
- vscode (IDE) with 
  - ms-vscode.cmake-tools
  - ms-vscode.cpptools
- ctest or sometimes catch2 (unit testing)
- Markdown (documentation)

## Start developing

### Use a classic Nix overlay

#### Option 1: use the software

The [nokxpkgs](https://github.com/nokx5/nokxpkgs#add-nokxpkgs-to-your-nix-channel) channel and associate overlay can be imported with the `-I` command or by setting the `NIX_PATH` environment variable.

```bash
# you can avoid this export by adding nokxpkgs to your channels
export NIX_PATH=nixpkgs=https://github.com/nokx5/nokxpkgs/archive/main.tar.gz
```

The package can then be used.

```bash
nix-shell -p golden_cpp
$ exit
```

#### Option 2: create a `default.nix` to build and develop the software

Create a simple `default.nix` from the channel.

```bash
$ cat default.nix 
{ pkgs ? import (builtins.fetchTarball
  "https://github.com/nokx5/nokxpkgs/archive/main.tar.gz") { } }:
with pkgs; {
  inherit golden_cpp;
}
```

You can then enter the shell or buid the project.

```bash
# option a: develop the local software
nix-shell -A golden_cpp
$ exit

# option b: build the local software
nix-build -A golden_cpp
unlink result
```

### Use an experimental flake feature

```bash
nix develop
```



```bash
nix build .#golden_cpp
```











## Code Snippets

One line code formatter for C/C++ projects

```bash
nixfmt $(find -name "*.nix")

clang-format -i $(find . -path "./build*" -prune  -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")

cmake-format -i $(find . -path "./build*" -prune  -name "*.cmake" -o -name "CMakeLists.txt")
```
