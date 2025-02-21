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
          golden-cpp = pkgs-super.callPackage ./derivation.nix {
            src = self;
            stdenv = pkgs-self.gccStdenv;
          };
          golden-cpp-clang = pkgs-super.callPackage ./derivation.nix {
            src = self;
            stdenv = pkgs-self.clangStdenv;
          };
          shell-dev = pkgs-super.callPackage ./shell.nix { pkgs = pkgs-self; };
        };
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ overlay ];
        };
      in {
        packages = with pkgs; {
          inherit golden-cpp golden-cpp-clang shell-dev;
        };
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

        devShell = self.packages.${system}.shell-dev;
      });
}
