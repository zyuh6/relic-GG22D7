/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (C) 20011-2019 RELIC Authors
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
 * Implementation of the final exponentiation for pairings over prime curves.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_pp.h"
#include "relic_util.h"

void pp_exp_k22(fp22_t c, fp22_t a) {
 	fp22_t t0, t1, f, g, m, m2;
	fp22_null(t0);
	fp22_null(t1);
	fp22_null(f);
	fp22_null(g);
	fp22_null(m2);
	fp22_null(m);
	RLC_TRY {
		fp22_new(t0);
		fp22_new(t1);
		fp22_new(f);
		fp22_new(g);
		fp22_new(m2);
		fp22_new(m);

		if(fp22_is_zero(a)) {
			fp22_zero(c);
			return;
		}

		/*easy part: m=a^{(p^11-1)(p+1)}= a^(p+1).conj / a^(p+1)*/
		fp22_frb(t0, a, 1);
		fp22_mul(t0, t0, a);
		fp22_inv_cyc(m, t0);
		fp22_inv(t0,t0);
		fp22_mul(m,m,t0);

		/* hard part*/
		fp22_sqr(m2,m);
		fp22_exp_cyc_x(a,m);
		fp22_exp_cyc_x(t0,a);
		fp22_mul(t1,m2,t0);
		fp22_inv_cyc(t0,a);
		fp22_mul(a,t1,t0);

		fp22_copy(f,a);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_mul(f,f,a);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_inv_cyc(t0,a);
		fp22_mul(f,f,t0);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_sqr(t0,a);
		fp22_mul(t0,t0,a);
		fp22_inv_cyc(t0,t0);
		fp22_mul(f,f,t0);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_inv_cyc(t0,a);
		fp22_mul(f,f,t0);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_mul(f,f,t0);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_inv_cyc(t1,a);
		fp22_mul(t0,t0,t1);
		fp22_mul(f,f,t0);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(f,f,1);
		fp22_sqr(t0,a);
		fp22_mul(t0,t0,a);
		fp22_inv_cyc(t0,t0);
		fp22_mul(f,f,t0);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(g,a);
		fp22_sqr(g,g);
		fp22_sqr(g,g);
		fp22_sqr(g,g);
		fp22_mul(g,a,g);
		fp22_frb(g,g,9);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(f,f);
		fp22_frb(t0,g,1);
		fp22_mul(f,f,t0);
		fp22_sqr(f,f);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_frb(t0,t0,9);
		fp22_mul(f,f,t0);
		fp22_sqr(f,f);

		fp22_exp_cyc_x(a,a);
		fp22_sqr(t0,m2);
		fp22_mul(t0,m,t0);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,m);
		fp22_mul(a,a,t0);

		fp22_sqr(g,a);
		fp22_mul(g,a,g);
		fp22_sqr(g,g);
		fp22_sqr(g,g);
		fp22_sqr(g,g);
		fp22_inv_cyc(g,g);
		fp22_mul(g,g,a);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_mul(g,g,t0);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_mul(g,g,t0);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_sqr(t0,a);
		fp22_mul(t0,t0,a);
		fp22_mul(g,g,t0);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_sqr(t0,t0);
		fp22_inv_cyc(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_mul(g,g,t0);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_sqr(t0,a);
		fp22_sqr(t0,t0);
		fp22_mul(t0,t0,a);
		fp22_inv_cyc(t0,t0);
		fp22_mul(g,g,t0);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_mul(g,g,a);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_sqr(t0,a);
		fp22_mul(t0,t0,a);
		fp22_mul(g,g,t0);

		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_mul(g,g,a);

		fp22_mul(f,f,g);

		fp22_exp_cyc_x(a,a);
		fp22_copy(g,a);
		fp22_exp_cyc_x(a,a);
		fp22_frb(g,g,1);
		fp22_mul(g,g,a);
		fp22_frb(g,g,9);
		fp22_mul(f,f,g);

		fp22_copy(c,f);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp22_free(t0);
		fp22_free(t1);
		fp22_free(f);
		fp22_free(g);	
		fp22_free(m2);	
		fp22_free(m);	
	}

}



