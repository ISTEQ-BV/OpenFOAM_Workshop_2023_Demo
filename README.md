# OpenFOAM_Workshop_2023_Demo
Demo code accompanying the talk "Implementing memory locality optimizations in OpemFOAM based code"
by Ilya Popov and Dmitry Astakhov, ISTEQ BV, 2023

This is demo not suitable for real world use yet.

## Installation

1. Clone code repository

        git clone https://github.com/ISTEQ-BV/OpenFOAM_Workshop_2023_Demo.git
        cd OpenFOAM_Workshop_2023_Demo

2. Checkout git submodules

        git submodule update --init

3. Activate your OpenFOAM environment

        . ~/OpenFOAM/OpenFOAM-v2212/etc/bashrc

4. Compile code

        cd src
        wmake
        cd ../

5. Run the nechmark in the test case directory

        cd test_case
        ./run.sh
