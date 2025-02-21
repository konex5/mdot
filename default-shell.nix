{ pkgs ? import <nixpkgs> { }, clangSupport ? true, cudaSupport ? false }:

with pkgs;

let stdenv = if clangSupport then clangStdenv else gccStdenv;
in (mkShell.override { inherit stdenv; }) rec {

  nativeBuildInputs = [ catch2 cmake gdb gnumake ninja pkg-config ] # valgrind
    ++ lib.optional (!clangSupport) [ gcovr lcov ]; # [ sssd cacert ]

  buildInputs = [
    # blas hdf5 tbb
    # tbb
    boost
    # eigen
  ] # ++ lib.optionals cudaSupport [ cudatoolkit.cc ] # cudatoolkit nvidia_x11
    ++ lib.optional stdenv.hostPlatform.isDarwin fixDarwinDylibNames;

  # shellHook = ''
  #   export SSL_CERT_FILE=${cacert}/etc/ssl/certs/ca-bundle.crt
  #   export PATH=/nix/var/nix/profiles/per-user/$USER/tools-dev/bin/:$PATH
  # '';
}
