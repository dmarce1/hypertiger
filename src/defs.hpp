/*
 * types.hpp
 *
 *  Created on: May 26, 2015
 *      Author: dmarce1
 */

#include <hpx/config.hpp>

//#define OCTOTIGER_RESTART_LOAD_SEQ

//#define OCTOTIGER_USE_NODE_CACHE

#ifdef OCTOTIGER_HAVE_GRAV_PAR
# define USE_GRAV_PAR
#endif

#define icoarse(i,j,k) ((k) * (NX * NX / 4) + (j) * (NX / 2) + (i))

#ifdef OCTOTIGER_HAVE_SILO
# define DO_OUTPUT
#endif

//#define OCTOTIGER_HAVE_RADIATION

#ifdef OCTOTIGER_HAVE_RADIATION
#define RADIATION
#endif

#ifndef TYPES444_HPP_

//#define CWD
#define BIBI
//#define OLD_SCF

//#define WD_EOS

#define EXPERIMENT
#ifdef RADIATION
#define NRF 4
#else
#define NRF 0
#define NRADF 0
#endif

#define NPF 5

#define abort_error() printf( "Error in %s on line %i\n", __FILE__, __LINE__); abort()

#define USE_SIMD

//#define USE_DRIVING
#define DRIVING_RATE 0.01
#define DRIVING_TIME 1.00

#define ZCORRECT

#define USE_PPM
//#define USE_MINMOD

#if !defined(OCTOTIGER_FORCEINLINE)
#   if defined(__NVCC__) || defined(__CUDACC__)
#       define OCTOTIGER_FORCEINLINE inline
#   elif defined(_MSC_VER)
#       define OCTOTIGER_FORCEINLINE __forceinline
#   elif defined(__GNUC__)
#       define OCTOTIGER_FORCEINLINE inline __attribute__ ((__always_inline__))
#   else
#       define OCTOTIGER_FORCEINLINE inline
#   endif
#endif

#include "real.hpp"
typedef long long int integer;
#define TYPES444_HPP_

typedef unsigned char byte;

enum gsolve_type {
	RHO, DRHODT
};

#include <array>
#include <iostream>

#define USE_ROTATING_FRAME
//#define OUTPUT_FREQ (100.0)

//#define USE_SPHERICAL
constexpr integer M_POLES = 3;
constexpr integer L_POLES = M_POLES;

//#define GRID_SIZE real(2.0)

constexpr real DEFAULT_OMEGA = 0.0;

//const integer MAX_LEVEL = 5;

enum boundary_type {
	OUTFLOW, REFLECT
};

constexpr integer NDIM = 3;

constexpr integer NSPECIES = 5;

constexpr integer INX = 8;
constexpr integer BW = 2;
constexpr integer NX = INX + 2 * BW;
constexpr integer N3 = NX * NX * NX;
constexpr integer DN[NDIM] = { 1, NX, NX * NX };

constexpr integer NDIR = 27;

constexpr integer XDIM = 0;
constexpr integer YDIM = 1;
constexpr integer ZDIM = 2;

constexpr integer FXM = 0;
constexpr integer FXP = 1;
constexpr integer FYM = 2;
constexpr integer FYP = 3;
constexpr integer FZM = 4;
constexpr integer FZP = 5;

constexpr integer NFACE = 2 * NDIM;
constexpr integer NVERTEX = 8;
constexpr integer NCHILD = 8;

constexpr real ZERO = real(0);
constexpr real ONE = real(1);
constexpr real TWO = real(2);
constexpr real THREE = real(3);
constexpr real FOUR = real(4);

constexpr real HALF = real(real(1) / real(2));
constexpr real SIXTH = real(real(1) / real(6));
constexpr real TWELFTH = real(real(1) / real(12));

constexpr real cfl = real(real(2) / real(15));
constexpr integer NRK = 2;
constexpr real rk_beta[2] = { ONE, HALF };

template<typename T>
constexpr inline T sqr(T const& val) {
	return val * val;
}

template<typename T>
constexpr inline T cube(T const& val) {
	return val * val * val;
}

template<typename T>
constexpr inline T average(T const& s1, T const& s2) {
	return 0.5 * (s1 + s2);
}
;

template<typename T>
inline void inplace_average(T& s1, T& s2) {
	s1 = s2 = average(s1, s2);
}
;

template<typename T>
std::size_t write(std::ostream& strm, T && t) {
	typedef typename std::decay<T>::type output_type;
	strm.write(reinterpret_cast<char const*>(&t), sizeof(output_type));
	return sizeof(output_type);
}

template<typename T>
std::size_t write(std::ostream& strm, T* t, std::size_t size) {
	strm.write(reinterpret_cast<char const*>(t), sizeof(T) * size);
	return sizeof(T) * size;
}

template<typename T>
std::size_t read(std::istream& strm, T & t) {
	typedef typename std::decay<T>::type input_type;
	strm.read(reinterpret_cast<char*>(&t), sizeof(input_type));
	return sizeof(input_type);
}

template<typename T>
std::size_t read(std::istream& strm, T* t, std::size_t size) {
	strm.read(reinterpret_cast<char*>(t), sizeof(T) * size);
	return sizeof(T) * size;
}

#endif /* TYPES_HPP_ */
