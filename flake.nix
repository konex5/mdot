{
  description = "golden-cpp - A simple C/C++ flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  };

  outputs = { self, nixpkgs }:

    let
      forCustomSystems = custom: f: nixpkgs.lib.genAttrs custom (system: f system);
      allSystems = [ "x86_64-linux" "i686-linux" "aarch64-linux" "x86_64-darwin" ];
      devSystems = [ "x86_64-linux" "x86_64-darwin" ];
      forAllSystems = forCustomSystems allSystems;
      forDevSystems = forCustomSystems devSystems;

      nixpkgsFor = forAllSystems (system:
        import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          overlays = [ self.overlay ];
        }
      );

      repoName = "golden-cpp";
      repoVersion = nixpkgsFor.x86_64-linux.golden-cpp.version;
      repoDescription = "golden-cpp - A simple C/C++ flake";
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
      };

      devShell = forDevSystems (system:
        let pkgs = nixpkgsFor.${system}; in pkgs.callPackage ./shell.nix { clangSupport = false; }
      );

      hydraJobs = forDevSystems (system: {

        build = nixpkgsFor.${system}.golden-cpp;
        build-clang = nixpkgsFor.${system}.golden-cpp-clang;

        docker = with nixpkgsFor.${system}; dockerTools.buildLayeredImage {
          name = "${repoName}-docker-${repoVersion}";
          tag = "latest";
          created = "now";
          contents = [ golden-cpp ];
          config = {
            Cmd = [ "cli_golden" ];
            # Env = [ "CMDLINE=ENABLED" ];
            # ExposedPorts = { "8000" = { }; };
          };
        };

        # deb = with nixpkgsFor.x86_64-linux; releaseTools.debBuild {
        #   name = "${repoName}-debian";
        #   diskImage = vmTools.diskImageFuns.debian8x86_64 {};
        #   src = golden-cpp.src;
        #   # buildInputs = [];
        # };

        # rpm = with nixpkgsFor.x86_64-linux; releaseTools.rpmBuild {
        #   name = "${repoName}-redhat";
        #   diskImage = vmTools.diskImageFuns.centos7x86_64 {};
        #   src = golden-cpp.src;
        #   # buildInputs = [];
        # };

        # tarball =
        #   nixpkgsFor.${system}.releaseTools.sourceTarball {
        #     name = "${repoName} - tarball";
        #     src = "./.";
        #   };

        # clang-analysis =
        #   nixpkgsFor.${system}.releaseTools.clangAnalysis {
        #     name = "${repoName}-clang-analysis";
        #     src = self;
        #   };

        # coverage =
        #   nixpkgsFor.${system}.releaseTools.coverageAnalysis {
        #     name = "${repoName}-coverage";
        #     src = self.hydraJobs.tarball;
        #     #lcovFilter = [ "*/tests/*" ];
        #   };

        release = nixpkgsFor.${system}.releaseTools.aggregate
          {
            name = "${repoName}-release-${repoVersion}";
            constituents =
              [
                self.hydraJobs.${system}.build
                #self.hydraJobs.${system}.build-clang                
                #self.hydraJobs.${system}.docker
                #self.hydraJobs.x86_64-linux.deb
                #self.hydraJobs.x86_64-linux.rpm
                #self.hydraJobs.${system}.tarball
                #self.hydraJobs.${system}.coverage
              ];
            meta.description = "hydraJobs: ${repoDescription}";
          };

      });

      packages = forAllSystems (system:
        with nixpkgsFor.${system}; {
          inherit golden-cpp golden-cpp-clang;
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

      templates = {
        golden-cpp = {
          description = "template - ${repoDescription}";
          path = ./.;
        };
      };

      defaultTemplate = self.templates.golden-cpp;
    };
}
