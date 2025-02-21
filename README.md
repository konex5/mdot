# Welcome to the golden cpp template

This is a skeleton template for a C/C++ project. Please find all the documentation [here](https://nokx5.github.io/golden_cpp) and the source code [here](https://github.com/nokx5/golden_cpp).

## My development tools are
- nix :snowflake: (packaging from hell :heart:)
- clang-format (formatter)
- vscode (IDE) with
  - ms-vscode.cmake-tools
  - ms-vscode.cpptools
- ctest - sometimes catch2 (unit testing)
- Markdown (documentation)

## Use the classic Nix commands

### Use the software (without git clone)

The [nokxpkgs](https://github.com/nokx5/nokxpkgs#add-nokxpkgs-to-your-nix-channel) channel and associate overlay can be imported with the `-I` command or by setting the `NIX_PATH` environment variable.

```bash
nix-shell -I nixpkgs=https://github.com/nokx5/nokxpkgs/archive/main.tar.gz -p golden_cpp --command cli_golden
```

### Develop the software

Start by cloning the [git repository](https://github.com/nokx5/golden_cpp) locally and enter it.

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
nix-shell derivation-shell.nix
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
nix shell github:nokx5/golden_cpp --command cli_golden
```

### Develop the software

Start by cloning the [git repository](https://github.com/nokx5/golden_cpp) locally and enter it.

#### Option 1: Develop the software

```bash
# option a: develop with a local shell
nix run .#cli_golden

# option b: build the local project
nix build .#golden_cpp
```

#### Option 2: Develop the software (supercharged :artificial_satellite:)

You can enter the development supercharged environment.

```bash
nix develop
```

## Code Snippets

One line code formatter for C/C++ projects

```bash
nixfmt $(find -name "*.nix")

clang-format -i $(find . -path "./build*" -prune  -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")

cmake-format -i $(find . -path "./build*" -prune  -name "*.cmake" -o -name "CMakeLists.txt")
```
