#!/bin/bash
cmake -B build
ninja -C build
./build/tests/bin/numericals_tests

