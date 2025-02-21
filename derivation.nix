{ stdenv, src, boost17x, catch2, cmake, cacert, curl, gnumake, hugo, ninja
, pandoc }:

stdenv.mkDerivation rec {
  pname = "golden_cpp";
  version = "0.0.1";
  inherit src;

  buildInputs = [ boost17x ];
  nativeBuildInputs = [ catch2 cmake gnumake ninja ] ++ [ cacert curl hugo ];
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

  outputs = [ "doc" "out" ]; # "bin" "dev"  "lib"
  preInstall = ''
    mkdir -p $doc/html
    make doc-build -C $src/docs DESTINATION=$TMP/tmp-doc HUGO_OPTS="--baseURL $doc/html/"
    mv $TMP/tmp-doc/public/* $doc/html/
  '';
}
