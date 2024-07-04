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
 * Implementation of point normalization on prime elliptic curves over an undecic
 * extension field.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_ADD == PROJC || !defined(STRIP)

/**
 * Normalizes a point represented in projective coordinates.
 *
 * @param r			- the result.
 * @param p			- the point to normalize.
 */
static void ep11_norm_imp(ep11_t r, const ep11_t p, int inverted) {
	if (p->coord != BASIC) {
		fp11_t t0, t1;

		fp11_null(t0);
		fp11_null(t1);

		RLC_TRY {

			fp11_new(t0);
			fp11_new(t1);

			if (inverted) {
				fp11_copy(t1, p->z);
			} else {
				fp11_inv(t1, p->z);
			}
			fp11_sqr(t0, t1);
			fp11_mul(r->x, p->x, t0);
			fp11_mul(t0, t0, t1);
			fp11_mul(r->y, p->y, t0);
			fp11_set_dig(r->z, 1);
		}
		RLC_CATCH_ANY {
			RLC_THROW(ERR_CAUGHT);
		}
		RLC_FINALLY {
			fp11_free(t0);
			fp11_free(t1);
		}
	}

	r->coord = BASIC;
}

#endif /* EP_ADD == PROJC */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep11_norm(ep11_t r, const ep11_t p) {
	if (ep11_is_infty(p)) {
		ep11_set_infty(r);
		return;
	}

	if (p->coord == BASIC) {
		/* If the point is represented in affine coordinates, we just copy it. */
		ep11_copy(r, p);
	}
#if EP_ADD == PROJC || !defined(STRIP)
	ep11_norm_imp(r, p, 0);
#endif
}