/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2022 RELIC Authors
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
 * Implementation of addition on prime elliptic curves over an undecic
 * extension field.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

/**
 * Adds two points represented in affine coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param r					- the result.
 * @param s					- the resulting slope.
 * @param p					- the first point to add.
 * @param q					- the second point to add.
 */
static void ep11_add_basic_imp(ep11_t r, fp11_t s, const ep11_t p, const ep11_t q) {
	fp11_t t0, t1, t2;

	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);

	RLC_TRY {
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);

		/* t0 = x2 - x1. */
		fp11_sub(t0, q->x, p->x);
		/* t1 = y2 - y1. */
		fp11_sub(t1, q->y, p->y);

		/* If t0 is zero. */
		if (fp11_is_zero(t0)) {
			if (fp11_is_zero(t1)) {
				/* If t1 is zero, q = p, should have doubled. */
				ep11_dbl_slp_basic(r, s, p);
			} else {
				/* If t1 is not zero and t0 is zero, q = -p and r = infty. */
				ep11_set_infty(r);
			}
		} else {
			/* t2 = 1/(x2 - x1). */
			fp11_inv(t2, t0);
			/* t2 = lambda = (y2 - y1)/(x2 - x1). */
			fp11_mul(t2, t1, t2);

			/* x3 = lambda^2 - x2 - x1. */
			fp11_sqr(t1, t2);
			fp11_sub(t0, t1, p->x);
			fp11_sub(t0, t0, q->x);

			/* y3 = lambda * (x1 - x3) - y1. */
			fp11_sub(t1, p->x, t0);
			fp11_mul(t1, t2, t1);
			fp11_sub(r->y, t1, p->y);

			fp11_copy(r->x, t0);
			fp11_copy(r->z, p->z);

			if (s != NULL) {
				fp11_copy(s, t2);
			}

			r->coord = BASIC;
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
	}
}

#endif /* EP_ADD == BASIC */

/*============================================================================*/
	/* Public definitions                                                         */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

void ep11_add_basic(ep11_t r, const ep11_t p, const ep11_t q) {
	if (ep11_is_infty(p)) {
		ep11_copy(r, q);
		return;
	}

	if (ep11_is_infty(q)) {
		ep11_copy(r, p);
		return;
	}

	ep11_add_basic_imp(r, NULL, p, q);
}

void ep11_add_slp_basic(ep11_t r, fp11_t s, const ep11_t p, const ep11_t q) {
	if (ep11_is_infty(p)) {
		ep11_copy(r, q);
		return;
	}

	if (ep11_is_infty(q)) {
		ep11_copy(r, p);
		return;
	}

	ep11_add_basic_imp(r, s, p, q);
}

#endif

void ep11_sub(ep11_t r, const ep11_t p, const ep11_t q) {
	ep11_t t;

	ep11_null(t);

	if (p == q) {
		ep11_set_infty(r);
		return;
	}

	RLC_TRY {
		ep11_new(t);

		ep11_neg(t, q);
		ep11_add(r, p, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep11_free(t);
	}
}

void ep11_shared_add(ep11_t RR, fp11_t alpha, fp11_t ZTQ, const ep11_t T, const ep11_t Q){
	fp11_t ZZ, U, S, H, I, J, W, V, temp;
	dv11_t t0, t1;
	ep11_t R;
	fp11_null(ZZ);
    fp11_null(U);
    fp11_null(S);
    fp11_null(H);
    fp11_null(I);
    fp11_null(J);
    fp11_null(W);
    fp11_null(V);
	fp11_null(temp);
	dv11_null(t0);
	dv11_null(t1);
	ep11_null(R);
	RLC_TRY {
		fp11_new(ZZ);
		fp11_new(U);
		fp11_new(S);
		fp11_new(H);
		fp11_new(I);
		fp11_new(J);
		fp11_new(W);
		fp11_new(V);
		fp11_new(temp);
		dv11_new(t0);
		dv11_new(t1);
		ep11_new(R);

		fp11_sqr(ZZ, T->z);
		fp11_mul(U,Q->x,ZZ);
		fp11_mul(temp,T->z,ZZ);
		fp11_mul(S,Q->y,temp);
		fp11_sub(H,U,T->x);
		fp11_add(temp,H,H);
		fp11_sqr(I,temp);
		fp11_mul(J,H,I);
		fp11_sub(alpha,S,T->y);
		fp11_add(W,alpha,alpha);
		fp11_mul(V,T->x,I);
		fp11_sqr(temp,W);
		fp11_sub(R->x,temp,J);
		fp11_sub(R->x,R->x,V);
		fp11_sub(R->x,R->x,V);
		fp11_sub(temp,V,R->x);
		fp11_muln_low(t0,W,temp);
		fp11_add(temp,T->y,T->y);
		fp11_muln_low(t1,temp,J);
		fp11_subc_low(t0,t0,t1);
		fp11_rdc(R->y,t0);
		fp11_mul(ZTQ,T->z,H);
		fp11_add(R->z,ZTQ,ZTQ);
		
		R->coord = PROJC;

		ep11_copy(RR,R);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(ZZ);
		fp11_free(U);
		fp11_free(S);
		fp11_free(H);
		fp11_free(I);
		fp11_free(J);
		fp11_free(W);
		fp11_free(V);
		fp11_free(temp);
		dv11_free(t0);
		dv11_free(t1);
		ep11_free(R);
	}
}

