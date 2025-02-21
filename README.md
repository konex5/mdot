# golden_cpp

This is a skeleton for a cpp project template (called a golden project).

## My development tools are
- nix <3
- clang-format (formatting)
- vscode (IDE) with 
  - ms-vscode.cmake-tools
  - ms-vscode.cpptools
- ctest/?catch2? (tests)
- Markdown (docs)

## nix

1. Enter a nixpkgs shell (deprecated)
   ```bash
   nix-shell --pure default-shell.nix --arg clangSupport true
   ```

2. Enter [nokxpkgs](https://github.com/nokx5/nokxpkgs#add-nokxpkgs-to-your-nix-channel) workflow
    ```bash
    # you can avoid this export by adding nokxpkgs to your channels
    export NIX_PATH=nokxpkgs=https://github.com/nokx5/nokxpkgs/archive/main.tar.gz
    
    # option 1: develop the local software
    nix-shell -A dev
    $ exit
    
    # option 2: build the local software
    nix-build -A dev
    unlink result*
    
    # option 3: use your local changes
    nix-shell -I nixpkgs=$PWD -p dev
    $ exit
    ```

## format
One line code formatter for pure cpp projects

```bash
clang-format -i $(find . -path "./build*" -prune  -name "*.c" -o -name "*.cpp" -o -name "*.h" -o -name "*.hpp")

cmake-format -i $(find . -path "./build*" -prune  -name "*.cmake" -o -name "CMakeLists.txt")
```
