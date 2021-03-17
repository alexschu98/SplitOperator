# SplitOperator

[![C/C++ Build](https://github.com/alexschu98/SplitOperator/actions/workflows/c-cpp-build.yml/badge.svg)](https://github.com/alexschu98/SplitOperator/actions/workflows/c-cpp-build.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## About

A 1D, 2D and 3D implementation of the split operator algorithm for a particle in a harmonic potential. Based on the implementation from the [Arcane Algorithm Archive](https://www.algorithm-archive.org/contents/split-operator_method/split-operator_method.html). The extensions to 2D and 3D were done using the [Eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page).

## Build
This project uses [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) so make sure your compiler can find the needed header files. Please note the 3d version requires the Eigen Tensor class from Eigen-unsupported. 

```zsh
mkdir bin
make
cd bin
./splitop1d
./splitop2d
./splitop3d
```

## License

[MIT License](https://opensource.org/licenses/MIT)
