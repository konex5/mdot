{ pkgs ? import <nixpkgs> { config.allowUnfree = true; }, clangSupport ? false, cudaSupport ? false }:

with pkgs;
assert hostPlatform.isx86_64;

let
  mkCustomShell = mkShell.override { stdenv = if clangSupport then clangStdenv else gccStdenv; };
  dbg = if clangSupport then lldb else gdb;

  vscodeExt = vscode-with-extensions.override {
    vscodeExtensions = with vscode-extensions; [ bbenoist.nix eamodio.gitlens ms-vscode.cpptools ] ++
      vscode-utils.extensionsFromVscodeMarketplace [
        {
          name = "cmake";
          publisher = "twxs";
          version = "0.0.17";
          sha256 = "CFiva1AO/oHpszbpd7lLtDzbv1Yi55yQOQPP/kCTH4Y=";
        }
        {
          name = "cmake-tools";
          publisher = "ms-vscode";
          version = "1.11.26";
          sha256 = "wnR+C8h78JoCAR99i+pOX3v9wQxBjC9UJr+tou1RzMk=";
        }
        {
          name = "emacs-mcx";
          publisher = "tuttieee";
          version = "0.31.0";
          sha256 = "McSWrOSYM3sMtZt48iStiUvfAXURGk16CHKfBHKj5Zk=";
        }
      ];
  };

in
mkCustomShell {
  buildInputs = [
    # stdenv.cc.cc.lib
    # libcxxabi
    boost
    spdlog
    tbb
    # zlib
  ] ++ lib.optionals (hostPlatform.isLinux) [ glibcLocales ];

  nativeBuildInputs = [ cmake gnumake ninja ] ++ [
    bashCompletion
    bashInteractive
    cacert
    clang-tools
    cmake-format
    cmakeCurses
    dbg
    # cppcheck
    emacs-nox
    # fmt
    git
    # llvm
    nixpkgs-fmt
    pkg-config
  ] ++ lib.optionals (hostPlatform.isLinux) [
    #typora hugo
    vscodeExt
  ];

  LANG = "en_US.UTF-8";

  shellHook = ''
    mkdir -p $(pwd)/.trash_config
    export HOME=$(pwd)/.trash_config
    export SSL_CERT_FILE=${cacert}/etc/ssl/certs/ca-bundle.crt
  '';
}
