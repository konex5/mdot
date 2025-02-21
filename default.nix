{ pkgs ? import <nokxpkgs> { }, clangSupport ? true, cudaSupport ? false }:

with pkgs;

let
  golden_cpp = pkgs.golden_cpp.overrideAttrs (oldAttrs: rec {
    name = "local_dev";
    src = ./.; # pkgs.nix-gitignore.gitignoreSource [ ".git" ] ./.;
  });
  dev = pkgs.golden_cpp.overrideAttrs (oldAttrs: rec {
    # buildInputs = oldAttrs.buildInputs ++ [];
    nativeBuildInputs = oldAttrs.nativeBuildInputs
      ++ [ cacert cmakeCurses emacs-nox gdb git pkg-config ]
      ++ lib.optional (!clangSupport) [ gcc gcovr lcov ]; # vscode valgrind
    CC = if clangSupport then "clang" else "gcc";
    CXX = if clangSupport then "clang++" else "g++";
    shellHook = ''
      export SSL_CERT_FILE=${cacert}/etc/ssl/certs/ca-bundle.crt
    '';
  });
in pkgs // { inherit golden_cpp dev; }
