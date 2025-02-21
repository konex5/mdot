{
  description = "A basic flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils }:
    utils.lib.eachSystem [ "x86_64-linux" ] (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
        };
	
        golden_cpp = (with pkgs;
	  stdenv.mkDerivation rec {
            pname = "golden_cpp";
            version = "0.0.1";
          
            src = self;          
            buildInputs = [ boost ];
          
            nativeBuildInputs = [ catch2 cmake gnumake ninja ];
          
            cmakeFlags = [
              "-DCMAKE_BUILD_TYPE=Release"
              "-DPROJECT_TESTS=On"
              "-DPROJECT_SANDBOX=OFF"
            ];
            hardeningEnable = [ "format" "fortify" "pic" ];
            ninjaFlags = [ "-v" ];
            makeFlags = [ "VERBOSE=1" ];
          
            enableParallelBuilding = true;
            enableParallelChecking = true;
          
            doCheck = true;
          });
        dev = (with pkgs;
          golden_cpp.overrideAttrs (oldAttrs: rec {
            nativeBuildInputs = oldAttrs.nativeBuildInputs ++ [
              cacert
              clang-tools
              cmakeCurses
              gdb
              pkg-config
              hugo
              emacs-nox
              vscode
              nixfmt
            ];
            shellHook = ''
              export SSL_CERT_FILE=${cacert}/etc/ssl/certs/ca-bundle.crt
            '';
          }));

      in {
        packages = { inherit golden_cpp; };
        defaultPackage = self.packages.${system}.golden_cpp;
        devShell = dev;
      });
}
