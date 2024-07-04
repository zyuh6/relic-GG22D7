/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2019 RELIC Authors
 *
 * This file is part of RELIC. RELIC is legal property of its developers,
 * whose names are not listed here. Please refer to the COPYRIGHT file
 * for contact information.
 *
 * RELIC is free software; you can redistribute it and/or modify it under the
 * terms of the version 2.1 (or later) of the GNU Lesser General Public License
 * as published by the Free Software Foundation; or version 2.0 of the Apache
 * License as published by the Apache Software Foundation. See the LICENSE files
 * for more details.
 *
 * RELIC is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the LICENSE files for more details.
 *
 * You should have received a copy of the GNU Lesser General Public or the
 * Apache License along with RELIC. If not, see <https://www.gnu.org/licenses/>
 * or <https://www.apache.org/licenses/>.
 */

/**
 * @file
 *
 * Implementation of Miller doubling for curves of embedding degree 12.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_pp.h"
#include "relic_fp_low.h"
#include "relic_fpx_low.h"
#include "relic_util.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/
void pp_dbl_k22(fp22_t L, ep11_t RR, const ep11_t T, const ep_t P) {
	fp11_t YY, ZZ, beta, t0, t1;
	ep11_t R;

	fp11_null(YY);
	fp11_null(ZZ);
	fp11_null(beta);
	fp11_null(t0);
	fp11_null(t1);
	ep11_null(R);

	RLC_TRY {
		fp11_new(YY);
		fp11_new(ZZ);
		fp11_new(beta);
		fp11_new(t0);
		fp11_new(t1);
		ep11_new(R);

		ep11_shared_dbl(R, YY, ZZ, beta, T);
		for(int i=0;i<11;i++) fp_mul(t1[i],ZZ[i],P->x);
		fp11_mul_art(t0,T->x);
		fp11_sub(t1,t0,t1);
		fp11_mul(t1,beta,t1);
		fp11_mul_art(t0,YY);
		fp11_add(t0,t0,t0);
		fp11_sub(L[1],t1,t0);
		fp11_mul(L[0],ZZ,R->z);
		for(int i=0;i<11;i++) fp_mul(L[0][i],L[0][i],P->y);
		
		ep11_copy(RR,R);

	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(YY);
		fp11_free(ZZ);
		fp11_free(beta);
		fp11_free(t0);
		fp11_free(t1);
		ep11_free(R);
	}
}

void pp_dbl1_k22(fp22_t L, ep11_t RR, const ep11_t T, const ep_t P) {
	fp11_t YY, ZZ, beta, t0, t1;
	ep11_t R;

	fp11_null(YY);
	fp11_null(ZZ);
	fp11_null(beta);
	fp11_null(t0);
	fp11_null(t1);
	ep11_null(R);

	RLC_TRY {
		fp11_new(YY);
		fp11_new(ZZ);
		fp11_new(beta);
		fp11_new(t0);
		fp11_new(t1);
		ep11_new(R);
		ep11_shared_dbl1(R, YY, beta, T);
		fp11_mul_art(t0,T->x);
		fp_sub(t0[0],t0[0],P->x);
		fp11_mul(t1,beta,t0);
		fp11_mul_art(t0,YY);
		fp11_add(t0,t0,t0);
		fp11_sub(L[1],t1,t0);
		for(int i=0;i<11;i++) fp_mul(L[0][i],(R->z)[i],P->y);

		ep11_copy(RR,R);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(YY);
		fp11_free(ZZ);
		fp11_free(beta);
		fp11_free(t0);
		fp11_free(t1);
		ep11_free(R);
	}
}

void pp_dbl_k22_2iso(fp22_t L1, fp22_t L2, ep11_t RR, const ep11_t T, const ep_t P1, const ep_t P2) {
	fp11_t YY, ZZ, beta, B, C, D, E;
	ep11_t R;

	fp11_null(YY);
	fp11_null(ZZ);
	fp11_null(beta);
	fp11_null(B);
	fp11_null(C);
	fp11_null(D);
	fp11_null(E);
	ep11_null(R);

	RLC_TRY {
		fp11_new(YY);
		fp11_new(ZZ);
		fp11_new(beta);
		fp11_new(B);
		fp11_new(C);
		fp11_new(D);
		fp11_new(E);
		ep11_new(R);

		ep11_shared_dbl(R, YY, ZZ, beta, T);
		fp11_add(B,YY,YY);
		fp11_mul(C,beta,ZZ);
		fp11_mul(D,beta,T->x);
		fp11_mul(E,R->z,ZZ);
		for(int i=0;i<11;i++) fp_mul(L2[0][i],P2->y,E[i]);
		for(int i=0;i<11;i++) fp_mul(L1[0][i],P1->y,E[i]);
		fp11_sub(D,D,B);
		fp11_mul_art(D,D);
		for(int i=0;i<11;i++) fp_mul(B[i],C[i],P2->x);
		fp11_sub(L2[1],D,B);
		for(int i=0;i<11;i++) fp_mul(B[i],C[i],P1->x);
		fp11_sub(L1[1],D,B);

		ep11_copy(RR,R);

	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(YY);
		fp11_free(ZZ);
		fp11_free(beta);
		fp11_free(B);
		fp11_free(C);
		fp11_free(D);
		fp11_free(E);
		fp11_free(R);
	}
}

void pp_dbl1_k22_2iso(fp22_t L1, fp22_t L2, ep11_t RR, const ep11_t T, const ep_t P1, const ep_t P2) {
	fp11_t YY, ZZ, beta, B, C, D, E;
	ep11_t R;

	fp11_null(YY);
	fp11_null(ZZ);
	fp11_null(beta);
	fp11_null(B);
	fp11_null(C);
	fp11_null(D);
	fp11_null(E);
	ep11_null(R);

	RLC_TRY {
		fp11_new(YY);
		fp11_new(ZZ);
		fp11_new(beta);
		fp11_new(B);
		fp11_new(C);
		fp11_new(D);
		fp11_new(E);
		ep11_new(R);

		ep11_shared_dbl1(R, YY, beta, T);
		fp11_add(B,YY,YY);
		fp11_mul(D,beta,T->x);
		for(int i=0;i<11;i++) fp_mul(L2[0][i],P2->y,(R->z)[i]);
		for(int i=0;i<11;i++) fp_mul(L1[0][i],P1->y,(R->z)[i]);
		fp11_sub(D,D,B);
		fp11_mul_art(D,D);
		for(int i=0;i<11;i++) fp_mul(B[i],beta[i],P2->x);
		fp11_sub(L2[1],D,B);
		for(int i=0;i<11;i++) fp_mul(B[i],beta[i],P1->x);
		fp11_sub(L1[1],D,B);

		ep11_copy(RR,R);

	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(YY);
		fp11_free(ZZ);
		fp11_free(beta);
		fp11_free(B);
		fp11_free(C);
		fp11_free(D);
		fp11_free(E);
		fp11_free(R);
	}
}