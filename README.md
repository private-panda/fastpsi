### The source codes of "Efficient Private Set Intersection by Utilizing Oblivious Transfer Extension" in AsiaCCS'25.

## 1 The dependency: https://github.com/Visa-Research/volepsi

## 2 How to run this code?
### 1). Put the source code in volepsi and add 'add_subdirectory(fastpsi)' in CMakeLists.txt of volePSI
### 2). Compile volepsi as the 'README.md' of volePSI
### 3). `cd out/build/linux/fastpsi`
### 4). Running examples

for PSI

./fastpsi -role 0 -nnx 16 -nny 16 -v & ./fastpsi -role 1 -nnx 16 -nny 16 -v

./fastpsi -role 0 -nnx 20 -nny 20 -v & ./fastpsi -role 1 -nnx 20 -nny 20 -v

./fastpsi -role 0 -nnx 24 -nny 24 -v & ./fastpsi -role 1 -nnx 24 -nny 24 -v

for OPPRF

./fastpsi -role 0 -nnx 16 -nny 16 -v -theta 1 & ./fastpsi -role 1 -nnx 16 -nny 16 -v -theta 1

./fastpsi -role 0 -nnx 16 -nny 16 -v -theta 2 & ./fastpsi -role 1 -nnx 16 -nny 16 -v -theta 2

./fastpsi -role 0 -nnx 16 -nny 16 -v -theta 4 & ./fastpsi -role 1 -nnx 16 -nny 16 -v -theta 4

./fastpsi -role 0 -nnx 16 -nny 16 -v -theta 8 & ./fastpsi -role 1 -nnx 16 -nny 16 -v -theta 8
