# golden_cpp

This is a skeleton template for a C/C++ project.

## My development tools are
- nix :snowflake: (packaging from hell :heart:)
- clang-format (formatter)
- vscode (IDE) with 
  - ms-vscode.cmake-tools
  - ms-vscode.cpptools
- ctest - sometimes catch2 (unit testing)
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
nix-shell -p golden_cpp --command "cli_golden_cpp"
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

You can then enter the shell or build the project.

```bash
# option a: develop the local software
nix-shell -A golden_cpp
$ exit

# option b: build the local software
nix-build -A golden_cpp
unlink result
```

### Use the experimental flake feature

NOTE: this section requires the experimental `flake` and `nix-command` features. Please refer to the official documentation for nix flakes. We highly recommand you to have a look to nix flakes since the issue of channel pinning is locked in the `flake.lock` file.

#### Option 1: use the software

The package can be used easily with flakes.

```bash
nix shell github:nokx5/golden_cpp --command cli_golden_cpp
```

#### Option 2: build and develop the software

You can enter the shell or build the project with flakes in a very convenient way.

```bash
# option a: develop the local software
nix develop
$ exit

# option b: build the local software
nix build .#golden_cpp
unlink result
```

## Code Snippets

One line code formatter for C/C++ projects

```bash
nixfmt $(find -name "*.nix")

clang-format -i $(find . -path "./build*" -prune  -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")

cmake-format -i $(find . -path "./build*" -prune  -name "*.cmake" -o -name "CMakeLists.txt")
```
