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
  ];
  hardeningEnable = [ "format" "fortify" "pic" ];
  ninjaFlags = [ "-v" ];
  makeFlags = [ "VERBOSE=1" ];
  enableParallelBuilding = true;
  enableParallelChecking = true;
  doCheck = true;

  # outputs = [ "doc" "out" ]; # "bin" "dev"  "lib"
  # preInstall = ''
  #   export DESTINATION=$doc/tmp/target-doc
  #   export THEMEDIR=$doc/tmp/theme    
  #   # build doc
  #   make doc-build -C $src/docs DESTINATION=$DESTINATION THEMEDIR=$THEMEDIR HUGO_OPTS="--baseURL $doc/share/doc/html/"
  #   # install doc and clean
  #   mv $DESTINATION/public $doc/html
  #   rm -rf $doc/tmp
  # '';
}
