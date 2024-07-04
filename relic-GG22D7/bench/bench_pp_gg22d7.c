/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2010 RELIC Authors
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
 * Benchmarks for pairings defined over prime elliptic curves.
 *
 * @ingroup bench
 */

#include <stdio.h>

#include "relic.h"
#include "relic_bench.h"

static void oatepairing(void){
	ep_t P;
	ep11_t Q;
	fp22_t r;

	ep_null(P);
	ep11_null(Q);
	fp22_null(r);
	ep_new(P);
	ep11_new(Q);
	fp22_new(r);

	BENCH_RUN("pp_map_oatep_k22"){
		ep_rand(P);
		ep11_rand(Q);
		BENCH_ADD(pp_map_oatep_k22(r,Q,P));
	} BENCH_END;
	ep_free(P);
	ep11_free(Q);
	fp22_free(r);

}

static void oatepairing_2iso(void){
	ep_t P;
	ep11_t Q;
	fp22_t r;
	ep_null(P);
	ep11_null(Q);
	fp22_null(r);
	ep_new(P);
	ep11_new(Q);
	fp22_new(r);


	BENCH_RUN("pp_map_soatep_k22_2iso"){
		ep_rand(P);
		ep11_rand(Q);
		BENCH_ADD(pp_map_soatep_k22_2iso(r,Q,P));
	} BENCH_END;
	ep_free(P);
	ep11_free(Q);
	fp22_free(r);
}

static void ppexpk22(void) {
	fp22_t a,c;
	fp22_null(a);
	fp22_null(c);
	fp22_new(a);
	fp22_new(c);

	BENCH_RUN("final_exp"){
		fp22_rand(a);
		BENCH_ADD(pp_exp_k22(c, a));
	} BENCH_END;
	fp22_free(a);
	fp22_free(c);
}

int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}

	conf_print();

	util_banner("Benchmarks for the PP module:", 0);

	// if (ep_param_set_any_pairf() != RLC_OK) {
	// 	RLC_THROW(ERR_NO_CURVE);
	// 	core_clean();
	// 	return 0;
	// }

	if (ep_param_set_any() != RLC_OK) {
		RLC_THROW(ERR_NO_CURVE);
		core_clean();
		return 0;
	}

	ep_param_print();
	
	util_banner("final exponentiation and pairing:", 1);
	ppexpk22();
	oatepairing_2iso();
	oatepairing();

	core_clean();
	return 0;
}
