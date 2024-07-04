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
 * Implementation of doubling on elliptic prime curves over an undecic
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
 * Doubles a point represented in affine coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[out] s			- the resulting slope.
 * @param[in] p				- the point to double.
 */
static void ep11_dbl_basic_imp(ep11_t r, fp11_t s, const ep11_t p) {
	fp11_t t0, t1, t2;

	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);

	RLC_TRY {
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);

		/* t0 = 1/(2 * y1). */
		fp11_dbl(t0, p->y);
		fp11_inv(t0, t0);

		/* t1 = 3 * x1^2 + a. */
		fp11_sqr(t1, p->x);
		fp11_copy(t2, t1);
		fp11_dbl(t1, t1);
		fp11_add(t1, t1, t2);

		ep11_curve_get_a(t2);
		fp11_add(t1, t1, t2);

		/* t1 = (3 * x1^2 + a)/(2 * y1). */
		fp11_mul(t1, t1, t0);

		if (s != NULL) {
			fp11_copy(s, t1);
		}

		/* t2 = t1^2. */
		fp11_sqr(t2, t1);

		/* x3 = t1^2 - 2 * x1. */
		fp11_dbl(t0, p->x);
		fp11_sub(t0, t2, t0);

		/* y3 = t1 * (x1 - x3) - y1. */
		fp11_sub(t2, p->x, t0);
		fp11_mul(t1, t1, t2);

		fp11_sub(r->y, t1, p->y);

		fp11_copy(r->x, t0);
		fp11_copy(r->z, p->z);

		r->coord = BASIC;
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

void ep11_dbl_basic(ep11_t r, const ep11_t p) {
	if (ep11_is_infty(p)) {
		ep11_set_infty(r);
		return;
	}

	ep11_dbl_basic_imp(r, NULL, p);
}

void ep11_dbl_slp_basic(ep11_t r, fp11_t s, const ep11_t p) {
	if (ep11_is_infty(p)) {
		ep11_set_infty(r);
		return;
	}

	ep11_dbl_basic_imp(r, s, p);
}

#endif

void ep11_shared_dbl(ep11_t RR, fp11_t YY, fp11_t ZZ, fp11_t beta, const ep11_t T) {
	fp11_t S, t0, t1, t2;
	dv11_t tt0, tt1;
	ep11_t R;
	fp11_null(S);
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	dv11_null(tt0);
	dv11_null(tt1);
	ep11_null(R);

	RLC_TRY {
		fp11_new(S);
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		dv11_new(tt0);
		dv11_new(tt1);
		ep11_new(R);

		fp11_sqr(YY, T->y);
		fp11_sqr(ZZ, T->z);
		fp11_mul(S, T->x, YY);
		fp11_mul_invu(t0, ZZ);
		fp11_sqr(t2, T->x);
		fp11_sub(t1, T->x, t0);
		fp11_add(t0, T->x, t0);
		fp11_mul(t0, t0, t1);
		fp11_add(beta, t0, t0);
		fp11_add(beta, beta, t0);
		fp11_sqr(t0, beta);
		fp11_add(t1, S, S);
		fp11_add(t2, t1, t1);
		fp11_add(t1, t2, t2);
		fp11_sub(R->x, t0, t1);
		fp11_sub(t0, t2, R->x);
		fp11_muln_low(tt0, beta, t0);
		fp11_add(t0,YY,YY);
		fp11_sqrn_low(tt1,t0);
		fp11_addd_low(tt1,tt1,tt1);
		fp11_subc_low(tt0,tt0,tt1);
		fp11_rdc(R->y, tt0);
		fp11_add(t0,T->y,T->z);
		fp11_sqr(R->z, t0);
		fp11_sub(R->z,R->z,YY);
		fp11_sub(R->z,R->z,ZZ);

		ep11_copy(RR,R);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(S);
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		dv11_free(tt0);
		dv11_free(tt1);
		ep11_free(R);
	}
}

void ep11_shared_dbl1(ep11_t RR, fp11_t YY, fp11_t beta, const ep11_t T) {
	fp_t h;
	fp11_t S, t0, t1, t2;
	dv11_t tt0, tt1;
	ep11_t R;
	fp_null(h);
	fp11_null(S);
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	dv11_null(tt0);
	dv11_null(tt1);
	ep11_null(R);

	RLC_TRY {
		fp_new(h);
		fp11_new(S);
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		dv11_new(tt0);
		dv11_new(tt1);
		ep11_new(R);

		fp11_sqr(YY, T->y);
		fp11_mul(S, T->x, YY);
		for(int i=1;i<10;i++) fp_copy(t1[i],(T->x)[i]);
		for(int i=1;i<10;i++) fp_copy(t0[i],(T->x)[i]);
		ep11_curve_get_half(h);
		fp_sub(t1[10],(T->x)[10],h);
		fp_add_dig(t1[0],(T->x)[0],1);
		fp_add(t0[10],(T->x)[10],h);
		fp_sub_dig(t0[0],(T->x)[0],1);
		fp11_mul(t0, t0, t1);
		fp11_add(beta, t0, t0);
		fp11_add(beta, beta, t0);
		fp11_sqr(t0, beta);
		fp11_add(t1, S, S);
		fp11_add(t2, t1, t1);
		fp11_add(t1, t2, t2);
		fp11_sub(R->x, t0, t1);
		fp11_sub(t0, t2, R->x);
		fp11_muln_low(tt0, beta, t0);
		fp11_add(t0,YY,YY);
		fp11_sqrn_low(tt1,t0);
		fp11_addd_low(tt1,tt1,tt1);
		fp11_subc_low(tt0,tt0,tt1);
		fp11_rdc(R->y, tt0);
		fp11_add(R->z,T->y,T->y);
	
		ep11_copy(RR,R);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(h);
		fp11_free(S);
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		dv11_free(tt0);
		dv11_free(tt1);
		ep11_free(R);
	}
}
