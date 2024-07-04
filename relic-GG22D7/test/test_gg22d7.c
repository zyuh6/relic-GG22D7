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
 * Tests for extensions defined over prime fields.
 *
 * @ingroup test
 */

#include "relic.h"
#include "relic_test.h"


static int exponentiation22(void) {
	int code = RLC_ERR;
	fp22_t a, b, c;
	bn_t d;

	fp22_null(a);
	fp22_null(b);
	fp22_null(c);
	bn_null(d);

	RLC_TRY {
		fp22_new(a);
		fp22_new(b);
		fp22_new(c);
		bn_new(d);
        
		TEST_CASE("frobenius and exponentiation are consistent") {
			fp22_rand(a);
			fp22_frb(b, a, 1);
			d->sign = RLC_POS;
			d->used = RLC_FP_DIGS;
			dv_copy(d->dp, fp_prime_get(), RLC_FP_DIGS);
			fp22_exp(c, a, d);
			TEST_ASSERT(fp22_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
	}
	RLC_CATCH_ANY {
		util_print("FATAL ERROR!\n");
		RLC_ERROR(end);
	}
	code = RLC_OK;
  end:
	fp22_free(a);
	fp22_free(b);
	fp22_free(c);
	bn_free(d);
	return code;
}

static int shared_add_dbl22(void) {
	int code = RLC_ERR;
	ep_t P, PP;
	ep11_t a, b, c, aa;
	fp22_t L1, L2;

	ep_null(P);
	ep_null(PP);
	ep11_null(a);
	ep11_null(b);
	ep11_null(c);
	ep11_null(aa);
	fp22_null(L1);
	fp22_null(L2);

	RLC_TRY{
		ep_new(P);
		ep_new(PP);
		ep11_new(a);
		ep11_new(b);
		ep11_new(c);
		ep11_new(aa);
		fp22_new(L1);
		fp22_new(L2);
		TEST_CASE("pp_dbl_k22_2iso is correct") {
			ep_curve_get_gen(P);
			ep_dbl_basic(PP,P);
			ep11_rand(a);
			ep11_dbl_basic(c,a);
			pp_dbl_k22_2iso(L1,L2,b,a,PP,P);
			
			TEST_ASSERT(ep11_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("pp_dbl1_k22_2iso is correct") {
			ep_curve_get_gen(P);
			ep11_rand(a);
			ep11_dbl_basic(c,a);
			pp_dbl1_k22_2iso(L1,L2,b,a,PP,P);
			
			TEST_ASSERT(ep11_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("pp_add_k22_2iso is correct") {
			ep_curve_get_gen(P);
			ep11_rand(a);
			ep11_dbl_basic(aa,a);
			ep_dbl_basic(PP,P);
			pp_add_k22_2iso(L2,L1,b,aa,a,PP,P);
			ep11_mul_dig(c,a,3);
			TEST_ASSERT(ep11_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("pp_dbl_k22 is correct") {
			ep_curve_get_gen(P);
			ep11_rand(a);
			ep11_dbl_basic(c,a);
			pp_dbl_k22(L1,b,a,P);
			
			TEST_ASSERT(ep11_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("pp_dbl1_k22 is correct") {
			ep_curve_get_gen(P);
			ep11_rand(a);
			ep11_dbl_basic(c,a);
			pp_dbl1_k22(L1,b,a,P);
			
			TEST_ASSERT(ep11_cmp(b, c) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("pp_add_k22 is correct") {
			ep_curve_get_gen(P);
			ep11_rand(a);
			ep11_dbl_basic(aa,a);
			pp_add_k22(L1,b,aa,a,P);
			ep11_mul_dig(c,a,3);
			TEST_ASSERT(ep11_cmp(b, c) == RLC_EQ, end);
		} TEST_END;
		
	}
	RLC_CATCH_ANY {
		util_print("FATAL ERROR!\n");
		RLC_ERROR(end);
	}
	code = RLC_OK;
  end:
	ep_free(P);
	ep_free(PP);
	ep11_free(a);
	ep11_free(b);
	ep11_free(c);
	ep11_free(aa);
	fp22_free(L1);
	fp22_free(L2);
	return code;
}

static int pairing22(void) {
	int code = RLC_ERR;
	ep11_t a;
	bn_t d, n;
	fp22_t r1,r2;
	ep_t P;

	ep11_null(a);
	bn_null(d);
	bn_null(n);
	fp22_null(r1);
	fp22_null(r2);
	ep_null(P);

	RLC_TRY {
		ep11_new(a);
		bn_new(d);
		bn_new(n);
		fp22_new(r1);
		fp22_new(r2);
		ep_new(P);

		TEST_CASE("pp_map_soatep_k22_2iso is bilinear") {
			ep11_rand(a);
			ep_rand(P);
			
			pp_map_soatep_k22_2iso(r1,a,P);
			bn_rand(d,RLC_POS,5);
			bn_rand(n,RLC_POS,6);
			ep11_mul_basic(a,a,d);
			ep_mul_basic(P,P,n);
			pp_map_soatep_k22_2iso(r2,a,P);

			bn_mul_basic(d,d,n);
			fp22_exp(r1,r1,d);

			TEST_ASSERT(fp22_cmp(r1, r2) == RLC_EQ, end);
		} TEST_END;

		TEST_CASE("pp_map_oatep_k22 is bilinear") {
			ep11_rand(a);
			ep_rand(P);
			
			pp_map_oatep_k22(r1,a,P);
			bn_rand(d,RLC_POS,3);
			bn_rand(n,RLC_POS,3);
			ep11_mul_basic(a,a,d);
			ep_mul_basic(P,P,n);
			pp_map_oatep_k22(r2,a,P);

			bn_mul_basic(d,d,n);
			fp22_exp(r1,r1,d);

			TEST_ASSERT(fp22_cmp(r1, r2) == RLC_EQ, end);
		} TEST_END;
		
	}
	RLC_CATCH_ANY {
		util_print("FATAL ERROR!\n");
		RLC_ERROR(end);
	}
	code = RLC_OK;
  end:
  	ep11_free(a);
	bn_free(d);
	bn_free(n);
	fp22_free(r1);
	fp22_free(r2);
	ep_free(P);
	return code;
}


int main(void) {
	if (core_init() != RLC_OK) {
		core_clean();
		return 1;
	}
	util_banner("Tests for the PP module", 0);

	if (ep_param_set_any() != RLC_OK) {
			core_clean();
			return 0;
	}
	fp_param_print();
	ep_param_print();

	util_banner("Arithmetic:", 1);

	if (exponentiation22() != RLC_OK) {
		core_clean();
		return 1;
	}

	if (shared_add_dbl22() != RLC_OK) {
		core_clean();
		return 1;
	}

	if (pairing22() != RLC_OK) {
		core_clean();
		return 1;
	}

	util_banner("All tests have passed.\n", 0);

	core_clean();
	return 0;
}
