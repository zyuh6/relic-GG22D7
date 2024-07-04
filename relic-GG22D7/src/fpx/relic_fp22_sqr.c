/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2012 RELIC Authors
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
 * Implementation of squaring in a 22-degree extension of a prime field.
 *
 * @ingroup fpx
 */

#include "relic_core.h"
#include "relic_fp_low.h"
#include "relic_fpx_low.h"


/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/
void fp22_sqr_basic(fp22_t c, const fp22_t a) {
	fp11_t t0, t1;

	fp11_null(t0);
	fp11_null(t1);

	RLC_TRY {
		fp11_new(t0);
		fp11_new(t1);

		fp11_add(t0, a[0], a[1]);
		fp11_mul_art(t1, a[1]);
		fp11_add(t1, a[0], t1);
		fp11_mul(t0, t0, t1);
		fp11_mul(c[1], a[0], a[1]);
		fp11_sub(c[0], t0, c[1]);
		fp11_mul_art(t1, c[1]);
		fp11_sub(c[0], c[0], t1);
		fp11_dbl(c[1], c[1]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t1);
	}
}


void fp22_sqr_lazyr(fp22_t c, const fp22_t a){
	fp11_t t0, t1;

	fp11_null(t0);
	fp11_null(t1);

	RLC_TRY {
		fp11_new(t0);
		fp11_new(t1);

		fp11_addn_low(t0, a[0], a[1]);
		fp11_mul_art(t1, a[1]);
		fp11_addn_low(t1, a[0], t1);
		fp11_mul(t0, t0, t1);
		fp11_mul(c[1], a[0], a[1]);
		fp11_sub(c[0], t0, c[1]);
		fp11_mul_art(t1, c[1]);
		fp11_sub(c[0], c[0], t1);
		fp11_dbl(c[1], c[1]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t1);
	}
}

void fp22_sqr_cyc(fp22_t c, const fp22_t a) {
    fp11_mul(c[1], a[0], a[1]);
    fp11_dbl(c[1], c[1]);
    fp11_sqr(c[0], a[0]);
    fp11_dbl(c[0], c[0]);
    fp_sub_dig(c[0][0], c[0][0], 1);
}
