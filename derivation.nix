{ stdenv
, boost
, cmakeMinimal
, fixDarwinDylibNames
, hdf5-cpp
, gnumake
, lapack
, lib
, nix-gitignore
, ninja
, tbb
, src ? ./.
, version ? "0.0.2"
}:

stdenv.mkDerivation {
  pname = "mdot";
  inherit version;
  inherit src;
  buildInputs = [ boost hdf5-cpp lapack tbb ];
  nativeBuildInputs = [ cmakeMinimal gnumake ninja ] ++ lib.optional stdenv.hostPlatform.isDarwin fixDarwinDylibNames;
  cmakeFlags = [
    "-DCMAKE_BUILD_TYPE=Release"
    "-DPROJECT_TESTS=On"
    "-DPROJECT_SANDBOX=OFF"
    "--no-warn-unused-cli"
  ];
  hardeningEnable = [ "format" "fortify" "pic" ];
  ninjaFlags = [ "-v" ];
  makeFlags = [ "VERBOSE=1" ];
  enableParallelBuilding = true;
  enableParallelChecking = true;
  doCheck = true;
}
