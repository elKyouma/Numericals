#!/bin/bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
ninja -C build
./build/tests/bin/numericals_tests

