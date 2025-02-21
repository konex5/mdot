{ stdenv, src, boost17x, catch2, cmake, curl, gnumake, hugo, ninja, pandoc }:

stdenv.mkDerivation rec {
  pname = "golden_cpp";
  version = "0.0.1";
  inherit src;

  buildInputs = [ boost17x ];
  nativeBuildInputs = [ catch2 cmake gnumake ninja ] ++ [ curl hugo ];
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


  outputs = [ "bin" "dev" "doc" "lib" "out"];
  preInstall = ''
    mkdir -p $bin $dev $doc $lib $out
    mkdir -p $out/share/doc/
    # install your doc here
    # make -C $src/docs DESTINATION=$out/share/doc/target-doc
'';
}
