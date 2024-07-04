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
void pp_add_k22(fp22_t L, ep11_t RR, const ep11_t T, const ep11_t Q, const ep_t P) {
	fp11_t alpha, ZTQ, t0;
	dv11_t u0, u1, u2;
	ep11_t R;
	fp11_null(alpha);
	fp11_null(ZTQ);
	fp11_null(t0);
	dv11_null(u0);
	dv11_null(u1);
	dv11_null(u2);
	ep11_null(R);

	RLC_TRY {
		fp11_new(alpha);
		fp11_new(ZTQ);
		fp11_new(t0);
		dv11_new(u0);
		dv11_new(u1);
		dv11_new(u2);
		ep11_new(R);

		ep11_shared_add(R, alpha, ZTQ, T, Q);
		fp11_mul_art(t0, Q->y);
		fp11_muln_low(u0, ZTQ, t0);
		fp11_mul_art(t0, Q->x);
		fp_sub(t0[0],t0[0],P->x);
		fp11_muln_low(u1,t0,alpha);
		fp11_subc_low(u0,u1,u0);
		fp11_rdc(L[1], u0);
		for(int i=0;i<11;i++) fp_mul(L[0][i], ZTQ[i], P->y);

		ep11_copy(RR,R);
	
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(alpha);
		fp11_free(ZTQ);
		fp11_free(t0);
		dv11_free(u0);
		dv11_free(u1);
		dv11_free(u2);
		ep11_free(R);
	}
}

void pp_add_k22_2iso(fp22_t L1, fp22_t L2, ep11_t RR, const ep11_t T, const ep11_t Q, const ep_t P1, const ep_t P2) {
	fp11_t alpha, ZTQ, t0;
	dv11_t u0, u1, u2;
	ep11_t R;
	fp11_null(alpha);
	fp11_null(ZTQ);
	fp11_null(t0);
	dv11_null(u0);
	dv11_null(u1);
	dv11_null(u2);
	ep11_null(R);

	RLC_TRY {
		fp11_new(alpha);
		fp11_new(ZTQ);
		fp11_new(t0);
		dv11_new(u0);
		dv11_new(u1);
		dv11_new(u2);
		ep11_new(R);

		ep11_shared_add(R, alpha, ZTQ, T, Q);
		fp11_mul_art(t0, Q->y);
		fp11_muln_low(u0, ZTQ, t0);
		fp11_mul_art(t0, Q->x);
		fp11_muln_low(u1, alpha, t0);
		fp11_subc_low(u0, u1, u0);
		for(int i=0;i<11;i++) fp_muln_low(u2[i], alpha[i], P2->x);
		fp11_subc_low(u2,u0,u2);
		fp11_rdc(L2[1], u2);
		for(int i=0;i<11;i++) fp_muln_low(u2[i], alpha[i], P1->x);
		fp11_subc_low(u2,u0,u2);
		fp11_rdc(L1[1], u2);
		for(int i=0;i<11;i++) fp_mul(L2[0][i], ZTQ[i], P2->y);
		for(int i=0;i<11;i++) fp_mul(L1[0][i], ZTQ[i], P1->y);

		ep11_copy(RR,R);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(alpha);
		fp11_free(ZTQ);
		fp11_free(t0);
		dv11_free(u0);
		dv11_free(u1);
		dv11_free(u2);
		ep11_free(R);
	}
}

