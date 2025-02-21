{ pkgs ? import <nixpkgs> { config.allowUnfree = true; }
, clangSupport ? false
, cudaSupport ? false
}:

with pkgs;
assert hostPlatform.isx86_64;

let
  vscodeExt = vscode-with-extensions.override {
    vscodeExtensions = with vscode-extensions;
      [ bbenoist.Nix eamodio.gitlens ms-vscode.cpptools ]
      ++ vscode-utils.extensionsFromVscodeMarketplace [
        {
          name = "cmake-tools";
          publisher = "ms-vscode";
          version = "1.7.3";
          sha256 = "6UJSJETKHTx1YOvDugQO194m60Rv3UWDS8UXW6aXOko=";
        }
        {
          name = "emacs-mcx";
          publisher = "tuttieee";
          version = "0.31.0";
          sha256 = "McSWrOSYM3sMtZt48iStiUvfAXURGk16CHKfBHKj5Zk=";
        }
        {
          name = "cmake";
          publisher = "twxs";
          version = "0.0.17";
          sha256 = "CFiva1AO/oHpszbpd7lLtDzbv1Yi55yQOQPP/kCTH4Y=";
        }
      ];
  };
  mkCustomShell = mkShell.override {
    stdenv = if clangSupport then clangStdenv else gccStdenv;
  };
in
mkCustomShell {
  nativeBuildInputs = [ catch2 cmake gnumake ninja ] ++ [
    bashCompletion
    cacert
    cppcheck
    clang-tools
    cmake-format
    cmakeCurses
    cmake-language-server
    fmt
    gdb
    git
    gnumake
    nixpkgs-fmt
    pkg-config
    emacs-nox
  ] ++ lib.optionals (hostPlatform.isLinux) [ typora vscodeExt ] ++ [ hugo ];
  buildInputs = [
    # libcxxabi
    boost17x
    spdlog
    tbb
    zlib
  ] ++ lib.optionals (hostPlatform.isLinux) [ glibcLocales ];

  shellHook = ''
    export SSL_CERT_FILE=${cacert}/etc/ssl/certs/ca-bundle.crt
  '';

}
