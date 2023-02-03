# DGFEM1D

One dimensional DG and FEM solver

## Prerequisites

* cmake
  ```sh
  sudo apt install cmake
  ```
* Eigen
  [Eigen](https://gitlab.com/libeigen/eigen)
  
## Installation

1. Clone the repo
   ```sh
   git clone git@github.com:erdemcaliskan/DGFEM1D.git
   ```
   
2. Build the project
    ```sh
    mkdir build && cd build
    cmake ../.
    cmake --build .
    ```