/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2023 RELIC Authors
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
 * Implementation of doubling on elliptic prime curves over an octic extension
 * field.
 *
 * @ingroup epx
 */

#include "relic_core.h"
#include "relic_ep_dbl_tmpl.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

/**
 * Doubles a point represented in affine coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param[out] r			- the result.
 * @param[out] s			- the slope.
 * @param[in] p				- the point to double.
 */
TMPL_DBL_BASIC_IMP(ep8, fp8);

#endif /* EP_ADD == BASIC */

#if EP_ADD == PROJC || !defined(STRIP)

/**
 * Doubles a point represented in projective coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param r					- the result.
 * @param p					- the point to double.
 */
TMPL_DBL_PROJC_IMP(ep8, fp8);

#endif /* EP_ADD == PROJC */

#if EP_ADD == JACOB || !defined(STRIP)

/**
 * Doubles a point represented in Jacobian coordinates on an ordinary prime
 * elliptic curve.
 *
 * @param r					- the result.
 * @param p					- the point to double.
 */
TMPL_DBL_JACOB_IMP(ep8, fp8);

#endif /* EP_ADD == JACOB */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

#if EP_ADD == BASIC || !defined(STRIP)

void ep8_dbl_basic(ep8_t r, const ep8_t p) {
	if (ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}
	ep8_dbl_basic_imp(r, NULL, p);
}

void ep8_dbl_slp_basic(ep8_t r, fp8_t s, const ep8_t p) {
	if (ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

	ep8_dbl_basic_imp(r, s, p);
}

#endif

#if EP_ADD == PROJC || !defined(STRIP)

void ep8_dbl_projc(ep8_t r, const ep8_t p) {
	if (ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

	ep8_dbl_projc_imp(r, p);
}

#endif

#if EP_ADD == JACOB || !defined(STRIP)

void ep8_dbl_jacob(ep8_t r, const ep8_t p) {
	if (ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

	ep8_dbl_jacob_imp(r, p);
}

#endif
