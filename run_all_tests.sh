#!/bin/bash
cmake -B build
ninja -D build
./build/tests/bin/numericals_tests

