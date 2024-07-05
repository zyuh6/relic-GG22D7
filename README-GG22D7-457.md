# GG22D7-457

### Algorithms

We implemented the pairing computation related to pairing-based protocols on "GG22D7-457" based on the famous [RELIC cryptographic toolkit](https://github.com/relic-toolkit/relic), including optimal ate pairing and super-optimal ate pairing.  

### Requirements

The build process requires the [CMake](https://cmake.org/) cross-platform build system. The [GMP](https://gmplib.org/) library is also needed in our benchmarks.

### Build instructions

Instructions for building the library can be found in the [Wiki](https://github.com/relic-toolkit/relic/wiki/Building).

### main functions
  
The main source code of our algorithms are distributed in different folders.  The main functions are:
* pp_map_oatep_k22 (fp22_t r, const ep11_t Q, const ep_t P): computes the optimal ate pairing of two given points $Q\in  \mathbb{G}_2$ and $P\in \mathbb{G}_1$.
* pp_map_soatep_k22_2iso (fp22_t r, const ep11_t Q, const ep_t P): computes the super-optimal ate pairing of two given points $Q\in  \mathbb{G}_2$ and $P\in \mathbb{G}_1$. 
* pp_exp_k22 (fp22_t c, fp22_t a): computes c = (a^(p^22 - 1)/r)^s, used for the final exponentiation in pairing computation. 
* fp22_exp_cyc_x (fp22_t c, const fp22_t a) : computes the x-th power of a cyclotomic 22-degree extension field element using addition chain, where x=-779523 is the seed. 

### finite field arithmitcs
     multiplications in $F_{p^11}$ : /src/fpx/relic_fp11_mul.c
     squairing in $F_{p^11}$:        /src/fpx/relic_fp11_sqr.c
     inversion in $F_{p^11}$:        line 580 in /src/fpx/relic_fpx_inv.c 
    
     multiplications in $F_{p^22}$ : /src/fpx/relic_fp22_mul.c
     squairing in $F_{p^22}$:        /src/fpx/relic_fp22_sqr.c
     inversion in $F_{p^22}$:        line 755 in /src/fpx/relic_fpx_inv.c 

### Testings, benckmarks and comparisons
* Testings and benckmarks: Function testing and benckmarking on GG22D7-P457 can be done by performing the following commands：

    1. mkdir build && cd build 
    2. ../preset/x64-pbc-gg22d7-457.sh ../
    3. make
    4. cd bin 
    5. ./test_gg22d7
    6. ./bench_pp_gg22d7

* 5 and 6 are to check that our implementation is correct and to obtain the clock cycles of pairing computation on GG22D7-P457, respectively. 
 
 *  For other curves, you can perform the following commands:
   1. mkdir build && cd build 
   2. ../preset/ < preset >.sh ./
   3. make
   4. cd bin 
   5. ./bench_pc
