#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w /io/wheelhouse/
    fi
}


# Install headers required for the extension
curl -sSL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz | tar -xvz
curl -sSL  https://github.com/Tencent/rapidjson/archive/refs/tags/v1.1.0.tar.gz | tar -xvz
curl -sSL https://boostorg.jfrog.io/artifactory/main/release/1.80.0/source/boost_1_80_0.tar.gz | tar -xvz
mkdir include_dir 
cp -r boost_*/boost include_dir
cp -r eigen-*/Eigen include_dir
cp -r rapidjson-*/include/rapidjson include_dir

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r /io/dev-requirements.txt
    "${PYBIN}/pip" wheel /io/ --no-deps -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    repair_wheel "$whl"
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install pyket -I ./include_dir --no-index --f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/nosetests" pyket)
done