{ stdenv, src, boost17x, catch2, cmake, gnumake, ninja, pandoc }:

stdenv.mkDerivation rec {
  pname = "golden_cpp";
  version = "0.0.1";
  inherit src;

  buildInputs = [ boost17x ];
  nativeBuildInputs = [ catch2 cmake gnumake ninja ] ++ [pandoc];
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
    for docinputs in $src/docs/*.md; do 
        pandoc $docinputs -f markdown -t html5 -H $src/docs/simple.css -s -o $out/share/doc/$(basename $docinputs .md).html --template $src/docs/template.html --toc --toc-depth=2 --lua-filter=$src/docs/links-to-html.lua
    done
'';
}
