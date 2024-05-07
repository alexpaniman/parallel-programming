{
  description = "Parallel programming";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {};
    in
  {

    packages.${system}.default = pkgs.clangStdenv.mkDerivation {
      name = "parallel-programming";
      src = ./.;

      nativeBuildInputs = with pkgs; [
        mpi
        clang-tools
      ];

      buildInputs = with pkgs; [];
    };

  };
}
