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
 * Implementation of multiplication in a 22-degree extension of a prime field.
 *
 * @ingroup fpx
 */

#include "relic_core.h"
#include "relic_fp_low.h"
#include "relic_fpx_low.h"


/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/
void fp22_muln_low(dv22_t c, const fp22_t a, const fp22_t b) {
    dv22_t t;
    fp22_t d;
    dv22_null(t);
    fp22_null(d);
    RLC_TRY {
        dv22_new(t);
        fp22_new(d);
        fp11_muln_low(t[0], a[0], b[0]);
        fp11_muln_low(t[1], a[1], b[1]);
        fp11_addn_low(d[0], a[0], a[1]);
        fp11_addn_low(d[1], b[0], b[1]);
        fp11_muln_low(c[1], d[0], d[1]);
        fp11_subc_low(c[1], c[1], t[0]);
        fp11_subc_low(c[1], c[1], t[1]);
        fp11_mul_nor_low(c[0], t[1]);
		fp11_addd_low(c[0],c[0],t[0]);
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        dv22_free(t);
        fp22_free(d);
    }
}

void fp22_mul_lazyr(fp22_t d, const fp22_t a, const fp22_t b) {
    dv22_t c;
    dv22_null(c);
	RLC_TRY{
		dv22_new(c);
		fp22_muln_low(c,a,b);
		fp22_rdc(d,c);
    } RLC_CATCH_ANY {
        RLC_THROW(ERR_CAUGHT);
    } RLC_FINALLY {
        dv22_free(c);
    }
}
