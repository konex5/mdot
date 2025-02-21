{ stdenv, src, boost17x, cmakeMinimal, fixDarwinDylibNames, gnumake, lapack, lib, ninja, tbb }:

# lapack spdlog

stdenv.mkDerivation rec {
  pname = "mdot";
  version = "0.0.2";
  name = "${pname}-${version}"
  inherit src;

  buildInputs = [ boost17x lapack tbb ];
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

  # outputs = [ "doc" "out" ]; # "bin" "dev"  "lib"
  # preInstall = ''
  #   mkdir -p $doc/html
  #   make doc-build -C $src/docs DESTINATION=$TMP/tmp-doc HUGO_OPTS="--baseURL $doc/html/"
  #   mv $TMP/tmp-doc/public/* $doc/html/
  # '';
}
