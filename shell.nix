{ pkgs ? import
    (
      builtins.fetchTarball {
        url = "https://github.com/NixOS/nixpkgs/archive/c00959877fb06b09468562518b408acda886c79e.tar.gz";
        sha256 = "sha256:02anj9mbzy45bszlmv7k42mb5y7yj2lxc5vpbxgd3f5qljn5ih7y";
      }
    )
    { }
, clangSupport ? true
}:

with pkgs;

assert hostPlatform.isx86_64;

let
  mkCustomShell = mkShell.override { stdenv = if clangSupport then clangStdenv else gccStdenv; };
  cc = if clangSupport then clangStdenv.cc else gccStdenv.cc;
  dbg = if clangSupport then lldb else gdb;
  together = symlinkJoin { name = "ccAndDbg"; paths = [ cc dbg ]; postBuild = "echo links added"; };
  name = if clangSupport then "clang++" else "g++";
in
mkCustomShell {
  buildInputs = [ boost17x hdf5-cpp lapack tbb ];
  nativeBuildInputs = [ together cmakeMinimal gnumake ninja ] ++ lib.optional stdenv.hostPlatform.isDarwin fixDarwinDylibNames;
  hardeningDisable = [ "all" ];
  shellHook = ''
    mkdir -p $(pwd)/.trash_config
    export HOME=$(pwd)/.trash_config
    export CPP=${together}/bin/${name}
  '';

}
