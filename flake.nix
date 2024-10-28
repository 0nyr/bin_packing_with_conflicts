{
  description = "Julia environment";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    flake-utils.url = github:numtide/flake-utils;
  };

  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        # Override the Nix package set to allow unfree packages
        pkgs = import nixpkgs {
          system = system; 
          config.allowUnfree = true; 
        };
        # WARN: Nix packaging system doesn't support all packages, so rely on Julia package manager instead.
        # Use Julia in REPL mode, then package mode and install packages that way.
        julia = pkgs.julia-bin.overrideDerivation (oldAttrs: { doInstallCheck = false; });
      in
      {
        # development environment
        devShells.default = pkgs.mkShell {
          packages = [
            julia
            pkgs.electron # for PlotlyJS, used as backend
            pkgs.gurobi
          ];

          shellHook = ''
            export JULIA_NUM_THREADS="auto"
            export JULIA_PROJECT="turing"
            export JULIA_BINDIR=${julia}/bin
            export JULIA_EDITOR="code"
            export LD_LIBRARY_PATH=${pkgs.gurobi}/lib:$LD_LIBRARY_PATH
            export GUROBI_HOME=${pkgs.gurobi}
            export GUROBI_VERSION=$(basename $(ls -d ${pkgs.gurobi}) | sed 's/.*-\([0-9]\+\)\.\([0-9]\+\).*/\1\2/')
            export ELECTRON_PATH=${pkgs.electron}/libexec/electron/electron
            echo "Nix shell loaded."
          '';
        };
      }
    );
}