{
  description = "A simple C/C++ flake";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils }:
    utils.lib.eachSystem [ "x86_64-linux" ] (system:
      let
        overlay = pkgs-self: pkgs-super: {
          project_gcc = pkgs-super.callPackage ./derivation.nix {
            src = self;
            stdenv = pkgs-self.gccStdenv;
          };
          project_clang = pkgs-super.callPackage ./derivation.nix {
            src = self;
            stdenv = pkgs-self.clangStdenv;
          };
          project_dev = pkgs-super.callPackage ./derivation-shell.nix { };
        };
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ overlay ];
        };
      in {
        packages = { golden-cpp = pkgs.project_gcc; };
        defaultPackage = self.packages.${system}.golden-cpp;

        apps = {
          cli_golden = {
            type = "app";
            program = "${self.defaultPackage.${system}}/bin/cli_golden";
          };
          cli_silver = {
            type = "app";
            program = "${self.defaultPackage.${system}}/bin/cli_silver";
          };
        };

        defaultApp = self.apps.${system}.cli_golden;

        devShell = pkgs.project_dev;
      });
}
