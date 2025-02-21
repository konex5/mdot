{
  description = "A simple C/C++ flake";

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
          let stdEnv = clangStdenv; # gccStdenv # clangStdenv
          in stdEnv.mkDerivation rec {
            pname = "golden_cpp";
            version = "0.0.1";
            src = self;
            buildInputs = [ boost17x ];
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
            nativeBuildInputs = let
              vscodeExt = vscode-with-extensions.override {
                vscodeExtensions = with vscode-extensions;
                  [ bbenoist.Nix ms-vscode.cpptools ]
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
            in [
              # stdenv.cc.cc
              # libcxxabi	      
              bashCompletion
              cacert
              clang-tools
              cmake-format
              cmakeCurses
              gdb
              pkg-config
              hugo
              emacs-nox
              nixfmt
              vscodeExt
              typora
            ] ++ oldAttrs.nativeBuildInputs;
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
