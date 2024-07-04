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
 * Implementation of point multiplication on prime elliptic curves over
 * an undecic extension field.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep11_mul_basic(ep11_t r, const ep11_t p, const bn_t k) {
	ep11_t t;
	int8_t u, naf[2 * RLC_FP_BITS + 1];
	size_t l;

	ep11_null(t);

	if (bn_is_zero(k) || ep11_is_infty(p)) {
		ep11_set_infty(r);
		return;
	}

	RLC_TRY {
		ep11_new(t);

		l = 2 * RLC_FP_BITS + 1;
		bn_rec_naf(naf, &l, k, 2);

		ep11_set_infty(t);
		for (int i = l - 1; i >= 0; i--) {
			ep11_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep11_add(t, t, p);
			} else if (u < 0) {
				ep11_sub(t, t, p);
			}
		}

		ep11_norm(r, t);
		if (bn_sign(k) == RLC_NEG) {
			ep11_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep11_free(t);
	}
}

void ep11_mul_dig(ep11_t r, const ep11_t p, const dig_t k) {
	ep11_t t;
	bn_t _k;
	int8_t u, naf[RLC_DIG + 1];
	size_t l;

	ep11_null(t);
	bn_null(_k);

	if (k == 0 || ep11_is_infty(p)) {
		ep11_set_infty(r);
		return;
	}

	RLC_TRY {
		ep11_new(t);
		bn_new(_k);

		bn_set_dig(_k, k);

		l = RLC_DIG + 1;
		bn_rec_naf(naf, &l, _k, 2);

		ep11_copy(t, p);
		for (int i = l - 2; i >= 0; i--) {
			ep11_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep11_add(t, t, p);
			} else if (u < 0) {
				ep11_sub(t, t, p);
			}
		}

		ep11_norm(r, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep11_free(t);
		bn_free(_k);
	}
}
