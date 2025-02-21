{
  description = "A simple C/C++ flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { self, nixpkgs }:

    let
      forCustomSystems = custom: f: nixpkgs.lib.genAttrs custom (system: f system);
      supportedSystems = [ "x86_64-linux" "i686-linux" "aarch64-linux" "darwin" ];
      forAllSystems = forCustomSystems supportedSystems;

      nixpkgsFor = forAllSystems (system:
        import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          overlays = [ self.overlay ];
        }
      );
      #pkgs = nixpkgsFor.${"x86_64-linux"};
      
      repoName = "golden-cpp";
      version = nixpkgsFor.${"x86_64-linux"}.golden-cpp.version;
    in
    {
      overlay = final: prev: {
        golden-cpp = prev.callPackage ./derivation.nix {
          src = self;
          stdenv = final.gccStdenv;
        };
        golden-cpp-clang = prev.callPackage ./derivation.nix {
          src = self;
          stdenv = final.clangStdenv;
        };
        fullDev = prev.callPackage ./shell.nix { pkgs = final; clangSupport = false; };
      };

      hydraJobs = forAllSystems (system: {

        build = self.packages.${system}.golden-cpp;

        tarball =
          nixpkgsFor.${system}.releaseTools.sourceTarball rec {
            name = "${repoName}-tarball";
            inherit version;
            src = self;
            postDist = ''
              cp README.md $out/
            '';
          };

        dockerImage = nixpkgsFor.${system}.dockerTools.buildLayeredImage {
          name = "${repoName}-docker";
          tag = "flake";
          created = "now";
          contents = [ self.defaultApp ];
          config = {
            Cmd = [ "cli_golden" ];
            # Env = [ "CMDLINE=ENABLED" ];
            # ExposedPorts = { "8000" = { }; };
          };
        };


        coverage =
          nixpkgsFor.${system}.releaseTools.coverageAnalysis {
            name = "${repoName}-coverage";
            src = self.hydraJobs.tarball;
            lcovFilter = [ "*/tests/*" ];
          };

        release = nixpkgsFor.${system}.releaseTools.aggregate
          {
            name = "${repoName}-${self.hydraJobs.tarball.version}";
            constituents =
              [
                self.hydraJobs.tarball
                self.hydraJobs.build.x86_64-linux
                self.hydraJobs.build.i686-linux
              ];
            meta.description = "Release golden-cpp";
          };

      });

      packages = forAllSystems (system:
        with nixpkgsFor.${system}; {
          inherit golden-cpp golden-cpp-clang fullDev;
        });

      defaultPackage = forAllSystems (system:
        self.packages.${system}.golden-cpp);

      apps = forAllSystems (system: {
        golden-cpp = {
          type = "app";
          program = "${self.packages.${system}.golden-cpp}/bin/cli_golden";
        };
        golden-cpp-clang = {
          type = "app";
          program = "${self.packages.${system}.golden-cpp-clang}/bin/cli_golden";
        };
      }
      );

      defaultApp = forAllSystems (system: self.apps.${system}.golden-cpp);

      templates = forAllSystems (system: {
        golden-cpp = {
          description = "C/C++ template";
          path = "${self.packages.${system}.golden-cpp}";
        };
        golden-cpp-clang = {
          description = "C/C++ template";
          path = "${self.packages.${system}.golden-cpp-clang}";
        };
      });

      defaultTemplate = self.templates.golden-cpp;
    };
}
