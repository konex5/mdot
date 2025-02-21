# Welcome to the golden cpp template

[![CI-linux](https://github.com/nokx5/golden-cpp/workflows/CI-linux/badge.svg)](https://github.com/nokx5/golden-cpp/actions/workflows/ci-linux.yml) [![CI-linux](https://github.com/nokx5/golden-cpp/workflows/CI-darwin/badge.svg)](https://github.com/nokx5/golden-cpp/actions/workflows/ci-darwin.yml) [![doc](https://github.com/nokx5/golden-cpp/workflows/doc-api/badge.svg)](https://nokx5.github.io/golden-cpp) [![MIT License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/nokx5/golden-cpp/blob/master/LICENSE)

This is a skeleton template for a C/C++ project. Please find all the documentation [here](https://nokx5.github.io/golden-cpp) and the source code [here](https://github.com/nokx5/golden-cpp).

## My development tools are
- nix :snowflake: (packaging from hell :heart:)
- clang-format (formatter)
- vscode (IDE) with
  - ms-vscode.cmake-tools
  - ms-vscode.cpptools
- ctest - or catch2 (unit testing)
- markdown (documentation)

## Use the classic Nix commands

### Use the software (without git clone)

The [nokxpkgs](https://github.com/nokx5/nokxpkgs#add-nokxpkgs-to-your-nix-channel) channel and associate overlay can be imported with the `-I` command or by setting the `NIX_PATH` environment variable.

```bash
nix-shell -I nixpkgs=https://github.com/nokx5/nokxpkgs/archive/main.tar.gz -p golden-cpp --command cli_golden
```

### Develop the software

Start by cloning the [git repository](https://github.com/nokx5/golden-cpp) locally and enter it.

#### Option 1: Develop the software (minimal requirements)

You can develop or build the local software easily with the minimal requirements.

```bash
# option a: develop with a local shell
nix-shell --expr '{pkgs ? import <nixpkgs> {} }: with pkgs; callPackage ./derivation.nix {src = ./.; }'

# option b: build the local project
nix-build --expr '{pkgs ? import <nixpkgs> {} }: with pkgs; callPackage ./derivation.nix {src = ./.; }'
```

Note that you can write the nix expression directly to the `default.nix` file to avoid typing `--expr` each time.

 #### Option 2: Develop the software (supercharged :artificial_satellite:)

You can enter the supercharged environment for development.

```bash
nix-shell shell.nix
```

## Use the experimental flake feature

**NOTE:** This section requires the experimental `flake` and `nix-command` features. Please refer to the official documentation for nix flakes. The advantage of using nix flakes is that you avoid channel pinning issues.

After Nix was installed, update to the unstable feature with:

```bash
nix-env -f '<nixpkgs>' -iA nixUnstable
```

And enable experimental features with:

```bash
mkdir -p ~/.config/nix
echo 'experimental-features = nix-command flakes' >> ~/.config/nix/nix.conf
```

### Use the software (without git clone)

```
nix shell github:nokx5/golden-cpp --command cli_golden
```

### Develop the software

Start by cloning the [git repository](https://github.com/nokx5/golden-cpp) locally and enter it.

#### Option 1: Develop the software

```bash
# option a: develop with a local shell
nix develop .#golden-cpp
# or
nix-shell . -A packages.x86_64-linux.golden-cpp

# option b: build the local project
nix build .#golden-cpp
# or
nix-build . -A packages.x86_64-linux.golden-cpp
```

#### Option 2: Develop the software (supercharged :artificial_satellite:)

You can enter the supercharged development environment.

```bash
nix develop
# or
nix-shell . -A devShell.x86_64-linux
```

## Code Snippets

One line code formatter for C/C++ projects

```bash
nixpkgs-fmt .

clang-format -i $(find . -path "./build*" -prune  -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")

cmake-format -i $(find . -path "./build*" -prune  -name "*.cmake" -o -name "CMakeLists.txt")
```
