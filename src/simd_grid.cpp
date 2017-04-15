/*
 * simd.cpp
 *
 *  Created on: Apr 12, 2017
 *      Author: dminacore
 */

#include "simd_grid.hpp"
#include "geometry.hpp"
#if 0
__attribute__((constructor))
static void unit_test() {
	printf("simd_grid test\n");
	constexpr integer NF = 3;
	simd_grid<NF> U;
	for (integer i = 0; i != NX; ++i) {
		for (integer j = 0; j != NX; ++j) {
			for (integer k = 0; k != NX; ++k) {
				U(i, j, k, 0) = k;
				U(i, j, k, 1) = i;
				U(i, j, k, 2) = j;
			}
		}
	}
	/*	for (integer i = 1; i != NX / 2 - 1; ++i) {
	 for (integer j = 1; j != NX / 2 - 1; ++j) {
	 for (integer k = 1; k != NX / 2 - 1; ++k) {
	 integer nf = 0;
	 auto uxp = U.get_shift(ZDIM, +1, nf + NF * icoarse(i, j, k));
	 auto uxm = U.get_shift(ZDIM, -1, nf + NF * icoarse(i, j, k));
	 auto u = U(nf + NF * icoarse(i, j, k));
	 for (integer j0 = 0; j0 != 2; ++j0) {
	 for (integer k0 = 0; k0 != 2; ++k0) {
	 const geo::octant c0( { j0, k0, 0 });
	 const geo::octant c1( { j0, k0, 1 });
	 printf("\n%e %e\n", uxp[c0], uxp[c1]);
	 printf("%e %e\n", u[c0], u[c1]);
	 printf("%e %e\n", uxm[c0], uxm[c1]);
	 }
	 }
	 }
	 }*/
	for (integer j = 1; j != NX / 2 - 1; ++j) {
		for (integer k = 1; k != NX / 2 - 1; ++k) {
			for (integer i = 1; i != NX / 2 - 1; ++i) {
				integer nf = 1;
				auto uxp = U.get_shift(XDIM, +1, nf + NF * icoarse(i, j, k));
				auto uxm = U.get_shift(XDIM, -1, nf + NF * icoarse(i, j, k));
				auto u = U(nf + NF * icoarse(i, j, k));
				for (integer j0 = 0; j0 != 2; ++j0) {
					for (integer k0 = 0; k0 != 2; ++k0) {
						const geo::octant c0( { 0, j0, k0 });
						const geo::octant c1( { 1, j0, k0 });
						printf("\n%e %e\n", uxp[c0], uxp[c1]);
						printf("%e %e\n", u[c0], u[c1]);
						printf("%e %e\n", uxm[c0], uxm[c1]);
					}
				}
			}
		}
	}
	for (integer k = 1; k != NX / 2 - 1; ++k) {
		for (integer i = 1; i != NX / 2 - 1; ++i) {
			for (integer j = 1; j != NX / 2 - 1; ++j) {
				integer nf = 2;
				auto uxp = U.get_shift(YDIM, +1, nf + NF * icoarse(i, j, k));
				auto uxm = U.get_shift(YDIM, -1, nf + NF * icoarse(i, j, k));
				auto u = U(nf + NF * icoarse(i, j, k));
				for (integer j0 = 0; j0 != 2; ++j0) {
					for (integer k0 = 0; k0 != 2; ++k0) {
						const geo::octant c0( { j0, 0, k0 });
						const geo::octant c1( { j0, 1, k0 });
						printf("\n%e %e\n", uxp[c0], uxp[c1]);
						printf("%e %e\n", u[c0], u[c1]);
						printf("%e %e\n", uxm[c0], uxm[c1]);
					}
				}
			}
		}
	}
	exit(0);
}

#endif
