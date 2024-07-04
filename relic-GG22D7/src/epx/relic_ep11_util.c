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
 * Implementation of comparison for points on prime elliptic curves over
 * an undecic extension field.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

int ep11_is_infty(const ep11_t p) {
	return (fp11_is_zero(p->z) == 1);
}

void ep11_set_infty(ep11_t p) {
	fp11_zero(p->x);
	fp11_zero(p->y);
	fp11_zero(p->z);
	p->coord = BASIC;
}

void ep11_copy(ep11_t r, const ep11_t p) {
	fp11_copy(r->x, p->x);
	fp11_copy(r->y, p->y);
	fp11_copy(r->z, p->z);
	r->coord = p->coord;
}

void ep11_rand(ep11_t p) {
	bn_t n, k;

	bn_null(k);
	bn_null(n);

	RLC_TRY {
		bn_new(k);
		bn_new(n);

		//ep11_curve_get_ord(n);
		//bn_rand_mod(k, n);
		ep11_curve_get_gen(p);
		//ep11_mul_gen(p, k);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(k);
		bn_free(n);
	}
}

void ep11_blind(ep11_t r, const ep11_t p) {
	fp11_t rand;

	fp11_null(rand);

	RLC_TRY {
		fp11_new(rand);
		fp11_rand(rand);
#if EP_ADD == BASIC
		(void)rand;
		ep11_copy(r, p);
#else
		fp11_mul(r->z, p->z, rand);
		fp11_mul(r->y, p->y, rand);
		fp11_sqr(rand, rand);
		fp11_mul(r->x, r->x, rand);
		fp11_mul(r->y, r->y, rand);
		r->coord = EP_ADD;
#endif
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(rand);
	}
}

void ep11_rhs(fp11_t rhs, const ep11_t p) {
	fp11_t t0, t1;

	fp11_null(t0);
	fp11_null(t1);

	RLC_TRY {
		fp11_new(t0);
		fp11_new(t1);
		fp11_sqr(t0, p->x);                
		fp11_mul(t0, t0, p->x);				
		ep11_curve_get_b(t1);
		fp11_add(t0, t0, t1);
		fp11_copy(rhs, t0);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t1);
	}
}


int ep11_on_curve(const ep11_t p) {
	ep11_t t;
	int r = 0;

	ep11_null(t);

	RLC_TRY {
		ep11_new(t);

		ep11_norm(t, p);

		ep11_rhs(t->x, t);
		fp11_sqr(t->y, t->y);

		r = (fp11_cmp(t->x, t->y) == RLC_EQ) || ep11_is_infty(p);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		ep11_free(t);
	}
	return r;
}
void ep11_print(const ep11_t p) {
	fp11_print(p->x);
	fp11_print(p->y);
	fp11_print(p->z);
}
