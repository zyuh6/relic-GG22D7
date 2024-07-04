/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2019 RELIC Authors
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
 * Implementation of configuration for prime field extensions.
 *
 * @ingroup fpx
 */

#include "relic_core.h"
#include "relic_fpx_low.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

int fp2_field_get_qnr() {
	/* Override some of the results when cubic non-residue is also needed. */
#if FP_PRIME == 1150
	return 32;
#elif FP_PRIME == 158 || FP_PRIME == 256
	return 4;
#elif FP_PRIME == 446 && !defined(FP_QNRES)
	return 16;
#else
	return core_get()->qnr2;
#endif
}

int fp3_field_get_cnr() {
#if FP_PRIME == 638
	if (fp_param_get() == K18_638) {
		return 8;
	} else {
		return 3;
	}
#elif FP_PRIME == 768
	return -4;
#endif
	return core_get()->cnr3;
}

void fp2_field_init(void) {
	bn_t e;
	fp2_t t0, t1;
	ctx_t *ctx = core_get();

	bn_null(e);
	fp2_null(t0);
	fp2_null(t1);

	RLC_TRY {
		bn_new(e);
		fp2_new(t0);
		fp2_new(t1);

		/* Start by finding a quadratic/cubic non-residue. */
#ifdef FP_QNRES
		ctx->qnr2 = 1;
#else
		/* First start with u as quadratic non-residue. */
		ctx->qnr2 = 0;
		fp_zero(t0[0]);
		fp_set_dig(t0[1], 1);
		/* If it does not work, attempt (u + 1), otherwise double. */
		/* We cannot used QR test here due to Frobenius constants below. */
		if (fp2_srt(t1, t0)) {
			ctx->qnr2 = 1;
			fp_set_dig(t0[0], ctx->qnr2);
			while (fp2_srt(t1, t0) && util_bits_dig(ctx->qnr2) < RLC_DIG - 1) {
				/* Pick a power of 2 for efficiency. */
				ctx->qnr2 *= 2;
				fp_set_dig(t0[0], ctx->qnr2);
			}
		}
#endif /* !FP_QNRES */

		/* Compute QNR^(p - 1)/6 and consecutive powers. */
		fp2_set_dig(t1, 1);
		fp2_mul_nor(t0, t1);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_sub_dig(e, e, 1);
		bn_div_dig(e, e, 6);
		fp2_exp(t0, t0, e);

		fp_copy(ctx->fp2_p1[0][0], t0[0]);
		fp_copy(ctx->fp2_p1[0][1], t0[1]);
		fp2_sqr(t1, t0);
		fp_copy(ctx->fp2_p1[1][0], t1[0]);
		fp_copy(ctx->fp2_p1[1][1], t1[1]);
		fp2_mul(t1, t1, t0);
		fp_copy(ctx->fp2_p1[2][0], t1[0]);
		fp_copy(ctx->fp2_p1[2][1], t1[1]);
		fp2_sqr(t1, t0);
		fp2_sqr(t1, t1);
		fp_copy(ctx->fp2_p1[3][0], t1[0]);
		fp_copy(ctx->fp2_p1[3][1], t1[1]);
		fp2_mul(t1, t1, t0);
		fp_copy(ctx->fp2_p1[4][0], t1[0]);
		fp_copy(ctx->fp2_p1[4][1], t1[1]);

		/* Compute QNR^(p - (p mod 4))/4. */
		fp2_set_dig(t1, 1);
		fp2_mul_nor(t0, t1);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 4);
		fp2_exp(t0, t0, e);
		fp_copy(ctx->fp2_p2[0][0], t0[0]);
		fp_copy(ctx->fp2_p2[0][1], t0[1]);

		/* Compute QNR^(p - (p mod 8))/8. */
		fp2_set_dig(t1, 1);
		fp2_mul_nor(t0, t1);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 8);
		fp2_exp(t0, t0, e);
		fp_copy(ctx->fp2_p2[1][0], t0[0]);
		fp_copy(ctx->fp2_p2[1][1], t0[1]);

		/* Compute QNR^(p - (p mod 12))/12. */
		fp2_set_dig(t1, 1);
		fp2_mul_nor(t0, t1);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 12);
		fp2_exp(t0, t0, e);
		fp_copy(ctx->fp2_p2[2][0], t0[0]);
		fp_copy(ctx->fp2_p2[2][1], t0[1]);

		/* Compute QNR^(p - (p mod 24))/24. */
		fp2_set_dig(t1, 1);
		fp2_mul_nor(t0, t1);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 24);
		fp2_exp(t0, t0, e);
		fp_copy(ctx->fp2_p2[3][0], t0[0]);
		fp_copy(ctx->fp2_p2[3][1], t0[1]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		bn_free(e);
		fp2_free(t0);
		fp2_free(t1);
	}
}

void fp3_field_init(void) {
	bn_t e;
	fp3_t t0, t1;
	ctx_t *ctx = core_get();

	bn_null(e);
	fp3_null(t0);
	fp3_null(t1);

	RLC_TRY {
		bn_new(e);
		fp3_new(t0);
		fp3_new(t1);

		/* Start by trying a trivial quadratic non-residue. */
		ctx->cnr3 = 0;
		fp_zero(t0[0]);
		fp_set_dig(t0[1], 1);
		fp_zero(t0[2]);
		/* If it does not work, attempt (u + 1), otherwise double. */
		/* This code will fail if p \neq 1 mod 8 because square root in Fp^3
		 * relies on Frobenius. Must implement explicit test for those cases. */
		if (fp3_srt(t1, t0)) {
			ctx->cnr3 = 1;
			fp_set_dig(t0[0], ctx->cnr3);
			while (fp3_srt(t1, t0) && util_bits_dig(ctx->qnr2) < RLC_DIG - 1) {
				/* Pick a power of 2 for efficiency. */
				ctx->cnr3 *= 2;
				fp_set_dig(t0[0], ctx->cnr3);
			}
		}

		/* Compute t0 = u^((p - (p mod 3))/3). */
		if (fp_prime_get_cnr() < 0) {
			fp_set_dig(ctx->fp3_p0[0], -fp_prime_get_cnr());
			fp_neg(ctx->fp3_p0[0], ctx->fp3_p0[0]);
		} else {
			fp_set_dig(ctx->fp3_p0[0], fp_prime_get_cnr());
		}
		bn_read_raw(e, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 3);
		fp_exp(ctx->fp3_p0[0], ctx->fp3_p0[0], e);
		fp_sqr(ctx->fp3_p0[1], ctx->fp3_p0[0]);

		/* Compute t0 = u^((p - (p mod 6))/6). */
		fp3_set_dig(t1, 1);
		fp3_mul_nor(t0, t1);
		bn_read_raw(e, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 6);
		fp3_exp(t0, t0, e);
		if (fp3_field_get_cnr() == 0) {
			/* Look for a non-trivial subfield element.. */
			ctx->frb3[0] = 0;
			while (ctx->frb3[0] < 3 && fp_is_zero(t0[ctx->frb3[0]++]));
			/* Fill rest of table with powers of constant. */
			fp_copy(ctx->fp3_p1[0][0], t0[--ctx->frb3[0] % 3]);
			fp3_sqr(t1, t0);
			fp_copy(ctx->fp3_p1[1][0], t1[(2 * ctx->frb3[0]) % 3]);
			fp3_mul(t1, t1, t0);
			fp_copy(ctx->fp3_p1[2][0], t1[(3 * ctx->frb3[0]) % 3]);
			fp3_sqr(t1, t0);
			fp3_sqr(t1, t1);
			fp_copy(ctx->fp3_p1[3][0], t1[(4 * ctx->frb3[0]) % 3]);
			fp3_mul(t1, t1, t0);
			fp_copy(ctx->fp3_p1[4][0], t1[(5 * ctx->frb3[0]) % 3]);
		} else {
			fp_copy(ctx->fp3_p1[0][0], t0[0]);
			fp_copy(ctx->fp3_p1[0][1], t0[1]);
			fp_copy(ctx->fp3_p1[0][2], t0[2]);
			fp3_sqr(t1, t0);
			fp_copy(ctx->fp3_p1[1][0], t1[0]);
			fp_copy(ctx->fp3_p1[1][1], t1[1]);
			fp_copy(ctx->fp3_p1[1][2], t1[2]);
			fp3_mul(t1, t1, t0);
			fp_copy(ctx->fp3_p1[2][0], t1[0]);
			fp_copy(ctx->fp3_p1[2][1], t1[1]);
			fp_copy(ctx->fp3_p1[2][2], t1[2]);
			fp3_sqr(t1, t0);
			fp3_sqr(t1, t1);
			fp_copy(ctx->fp3_p1[3][0], t1[0]);
			fp_copy(ctx->fp3_p1[3][1], t1[1]);
			fp_copy(ctx->fp3_p1[3][2], t1[2]);
			fp3_mul(t1, t1, t0);
			fp_copy(ctx->fp3_p1[4][0], t1[0]);
			fp_copy(ctx->fp3_p1[4][1], t1[1]);
			fp_copy(ctx->fp3_p1[4][2], t1[2]);
		}

		/* Compute t0 = u^((p - (p mod 9))/9). */
		fp3_set_dig(t1, 1);
		fp3_mul_nor(t0, t1);
		bn_read_raw(e, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 9);
		fp3_exp(t0, t0, e);
		/* Look for a non-trivial subfield element.. */
		if (fp3_field_get_cnr() == 0) {
			ctx->frb3[1] = 0;
			while (ctx->frb3[1] < 3 && fp_is_zero(t0[ctx->frb3[1]++]));
			fp_copy(ctx->fp3_p2[0][0], t0[--ctx->frb3[1]]);
		} else {
			fp_copy(ctx->fp3_p2[0][0], t0[0]);
			fp_copy(ctx->fp3_p2[0][1], t0[1]);
			fp_copy(ctx->fp3_p2[0][2], t0[2]);
		}

		/* Compute t0 = u^((p - (p mod 18))/18). */
		fp3_set_dig(t1, 1);
		fp3_mul_nor(t0, t1);
		bn_read_raw(e, fp_prime_get(), RLC_FP_DIGS);
		bn_div_dig(e, e, 18);
		fp3_exp(t0, t0, e);
		/* Look for a non-trivial subfield element.. */
		if (fp3_field_get_cnr() == 0) {
			ctx->frb3[2] = 0;
			while (ctx->frb3[2] < 3 && fp_is_zero(t0[ctx->frb3[2]++]));
			fp_copy(ctx->fp3_p2[1][0], t0[--ctx->frb3[2]]);
		} else {
			fp_copy(ctx->fp3_p2[1][0], t0[0]);
			fp_copy(ctx->fp3_p2[1][1], t0[1]);
			fp_copy(ctx->fp3_p2[1][2], t0[2]);
		}
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		bn_free(e);
		fp3_free(t0);
		fp3_free(t1);
	}
}

void fp4_field_init() {
	bn_t e;
	fp4_t t0;
	ctx_t *ctx = core_get();

	bn_null(e);
	fp4_null(t0);

	RLC_TRY {
		bn_new(e);
		fp4_new(t0);

		fp4_set_dig(t0, 1);
		fp4_mul_art(t0, t0);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_sub_dig(e, e, 1);
		bn_div_dig(e, e, 6);
		fp4_exp(t0, t0, e);
		if (fp2_is_zero(t0[1])) {
			ctx->frb4 = 0;
			fp_copy(ctx->fp4_p1[0], t0[0][0]);
			fp_copy(ctx->fp4_p1[1], t0[0][1]);
		} else {
			ctx->frb4 = 1;
			fp_copy(ctx->fp4_p1[0], t0[1][0]);
			fp_copy(ctx->fp4_p1[1], t0[1][1]);
		}
	} RLC_CATCH_ANY {
	    RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		bn_free(e);
		fp4_free(t0);
	}
}

void fp8_field_init() {
	bn_t e;
	fp8_t t0;
	ctx_t *ctx = core_get();

	bn_null(e);
	fp8_null(t0);

	RLC_TRY {
		bn_new(e);
		fp8_new(t0);

		fp8_set_dig(t0, 1);
		fp8_mul_art(t0, t0);
		e->used = RLC_FP_DIGS;
		dv_copy(e->dp, fp_prime_get(), RLC_FP_DIGS);
		bn_sub_dig(e, e, 1);
		bn_div_dig(e, e, 6);
		fp8_exp(t0, t0, e);
		if (fp4_is_zero(t0[1])) {
			ctx->frb8 = 0;
			fp_copy(ctx->fp8_p1[0], t0[0][0][0]);
			fp_copy(ctx->fp8_p1[1], t0[0][0][1]);
		} else {
			ctx->frb8 = 1;
			fp_copy(ctx->fp8_p1[0], t0[1][1][0]);
			fp_copy(ctx->fp8_p1[1], t0[1][1][1]);
		}
	} RLC_CATCH_ANY {
	    RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		bn_free(e);
		fp8_free(t0);
	}
}

/*The below are all about fp11 and fp22 for GG22D7_457*/
//For the Re of Fp22, that is, Fp11
#define GG22D7_1_1 "3b44a0f63bd7ac44a8a603d9a594a6f845c2eb798ec0a1e326e6758ac5b21871843ff0954ceddbebe7b10166013a45bf978d32d5251eb054aa"
#define GG22D7_1_2 "14d2f86417b372064f89f8600b56ccd92fa8318ff21da722bef9e4ae3b622a7a92a122e828d392472867318f0bd3e2c188fc75efbe22ae2851"
#define GG22D7_1_3 "4511cf31705eb0da34ee90e61c6d8b0646426fac2ca292f7ad128ff7f61481d0580be4da458212fad94fdc80db6caddc910c19e72a34e5001d"
#define GG22D7_1_4 "b2fa4630a2bf2000afd8cb3f5d7626b0fee99dc7188dfc89c43eb5a317e0dc35c74a4804e4a8f8f7be316279c682fcb876049cdc3877247923"
#define GG22D7_1_5 "fca92cf46a771e0b552d6c7221ad876d3d9ee67ce4357400a1b01d319fecad9d513fbaf92a26c795ace8923b83f9409cf286ca9e0d1e7da6f3"
#define GG22D7_1_6 "1664d8138030ed9392f9d35247e8a249000239e90579eb7fda946f6421407a76c279ead11fb16d3941b49850bd72d32512deaafa14faa5593ab"
#define GG22D7_1_7 "112c1795ec6e469ee9e02d5774ae7e204402e650ebf452517867dd29abbfc4e05df9a025592f7c1078de8b1a2643edcc4abc8ca112759ef0dd4"
#define GG22D7_1_8 "6d94cacb346b5dbcca81955b0cd25982910bb70aff939e707daeedbd90ae9c8b79ce72c4b87f21be1cc2c4b2d0b806b97fd7c7612e17fd39a8"
#define GG22D7_1_9 "19e062fd5da8ca8a0789b3d76941fd42f749883e3da320cf834ad3a8d23528c992afe7ead3395d832f17b8c1aa60c362c51fb4859e59a07841c"
#define GG22D7_1_10 "12b106642732993a63934cf0ad7384be0a59254236b8f46bf488348ba2f746d519ef5e5fd042f8e98fc17fef3c5a5f680f3a9bccc0582792558"
#define GG22D7_1_11 "1a16c3df64b7c6e364d1afce15f10810331420be7d9696ee4adb8a8ff137cc322e85d1b5dc990c55e3dde2d345a1557660ed86614e83881fe00"
#define GG22D7_2_1 "1e0c4a41fe17182678a1acdb30b15a0b720872d7aa99819c6a2ddc34dcaa7ad05f74ecab0b3c48676dce562d741901508486bbffb17e8d6291"
#define GG22D7_2_2 "2aae95ea113afcc194957a9755640a17b615f5f2793fcc415fc4230285b4709d8848b1b1124f383691347e4d0bacfff6d090c064786eab9ef5"
#define GG22D7_2_3 "8469a614ebfac42b52dae1d0d05c035821a64d8a2031f6c5075f5ac349b42678c6dcf148bc5645f9f43c3d4e5c77c68d656ea5a2044a5de34e"
#define GG22D7_2_4 "9acf8e6831013d19d6628b7ba7e0957ffd9f22aa63b9aa70552d1f5e0a161c46a3d2eed37c0983ab414098c0f2389676cf09f6b45869fff8bc"
#define GG22D7_2_5 "1bdd88ed6389a785f87edb7e94caae67cfd31ff3d48abac1d1f9d02525008e7c71ac09def11e7e2a4cc5a03a989a911dbad549fac38992d1e6a"
#define GG22D7_2_6 "1721ccdfbfbe0a06de9464c39cbb969f0aa2ee16f856cacc95263f9c259eb6d85bd5d393430d5d7282863ce0d3893c16fc5fa6c7e5edf58f6a2"
#define GG22D7_2_7 "1990bd4d8edc47f883b8d497e262fa622fb648185c08a750498562ef07416b5ad67a4968bd90b255676f7352d581ff0265b15b325453541bf1"
#define GG22D7_2_8 "13717d84ec25561a34dca6bbbaf466e00628363839e2ec911782dab0ac15e91352e32e14af0f92195e6e5917ae7a854b2e434cbfa16caf66cae"
#define GG22D7_2_9 "4efce3014b558c21f9aa6818ddf604ad0d1438459b5e98afbc651bf15eec724df78de68fa3adb9dc6ef7a01294a2e1bf3b7ad1a8d55b435179"
#define GG22D7_2_10 "726531c1e8828f72e79c7c7e2b71eaa9051c7d9d48ca0cdccdd6c635dd389d5f9338a46068266922f08f86505a4da6375e4a91189596c4e950"
#define GG22D7_2_11 "843cb940288538241d34e266dbd528a3e3c5d6e79e45043750e8d4a6168f60bf29b3b4373f93566b6e5aaab54e637baba9313efd18b78db15b"
#define GG22D7_3_1 "3987b9b715a10d73bdd3d0966dbfd5e42819aad438a345c5bbc9b6fbb52271096fde12960419f663137bcfcdd94b919560206219e70920164b"
#define GG22D7_3_2 "6f8164aa18820b9305af7c52e642964c59ab524dcb147c1f5cd72323521a4fabe72d86d5de4f3c94a9547f7492ff359b2f80bab24d18004b16"
#define GG22D7_3_3 "62854dcfe4ae717bee9c47dd3a7d8b895f36c5a388cfd84d2e2d6068fd44e5d083504d72b0e8ee815c792ecc6b2321960d37d68da1ad460e03"
#define GG22D7_3_4 "16356c41be35cd9389ce69c1f340a0bda099c0465b84013525446d07b3b58c509160ea6f49a9f3205f185322eae8d4d4ee2c6afd3e2dad7f12c"
#define GG22D7_3_5 "b64077ef9f53c939a6f1f75fe3b79ea0aa4137f0b5b610852c4c6577ae70e888473c6d01e30df82e0d66a742c546a55049eb87c9457190c4fd"
#define GG22D7_3_6 "17264b62cfe5835075d2c64ff8dd467809915299b949731d3304fc884dcabaeb62bc4c5671001ddbf52290544fbc202a47bc578d0c900df5576"
#define GG22D7_3_7 "5b509fc24811844c1d9463b3d5ab3591ac9b852d22ca4606fa57fc2f59b4c14b7377c73c092cf0bef92e6db44b61569fcfccfcc7372902d3d6"
#define GG22D7_3_8 "1baf09ad25a9b657f6b5caa773a1a175a94f8d28bff80ab5f934c3b43c45e4835ad7b5cced5b837ed533572fb72506420f0146976335d5c3c8"
#define GG22D7_3_9 "e81cf067233e82fcfe7aaa07f305adf591efa4960f16a7415b60144f53cbefd99edd54b3e1cbd84bca2d6347a621874f4735109203c524b519"
#define GG22D7_3_10 "188c1c9b14103a0c868e025845cff838a8bc2e7559aebba76af27c1962131048d157e7779063bc2558e3d9162a46eff5bd9cd5d128ed42c4891"
#define GG22D7_3_11 "f2e2b49adcb10ff3db876c6cfe876757bfc7986b90b3ec7c7c637bb29019de75eba5504b81d8f0504fc3d8bcc9f85e54130b5283b7c579074b"
#define GG22D7_4_1 "10d8ff5409977cf0a1ad6681761615561c7fb682046a326fc898c20f1d1733ed050d7d2f79c068b7e1cbb0f583cf0058f82bbda5fc285500bf6"
#define GG22D7_4_2 "68fcc7cb0926875bc5f03ed878f1d93c2b241c0184787a7c4a42bf77fe731d8c2862f646d898127710cef185d19546d1e2780e9cc07884c55d"
#define GG22D7_4_3 "17246427f747786136f065bdca9189d9a75833ed9ef2de2d69695968f703025e898bc596715f448c4a6e965d268cd2c722273b86a510bc4975b"
#define GG22D7_4_4 "e493c46c0d5576c97e349cc095d0fb5e9462b5c1a7eefcac3471c73d117fb7b7dc8b8cc70652b3992e5a05588e8f06648464f55dfe3878ef6b"
#define GG22D7_4_5 "a5021539d6a87444b9c50beccfcf7d1884e3eec06a774ec1b374575b834bcec04aed071af13cfd1f3ffa0a3c1b25a36e6087b4a7db6711e316"
#define GG22D7_4_6 "f8030b214596933ccd93b4f948da3d80bb026a91c71a9712d336d66451a7fcc7ee2f3807083506af8d2dadaf37b2f194e69a81bf49d4e6d493"
#define GG22D7_4_7 "1912904ccf8962526be87ee976f0f123575064dd87726b74542c4a68a66470e267afd917ab415189f19dd9c4a10aa63fbe9cd521a8ad92b2739"
#define GG22D7_4_8 "16af86f3e343987c6a55c91eada7bded5efc644a9426d75e014a3ebc9e6aa8f906be8fd31b7d3609138f5bb8aa6f2eafefa8adca270f3ad6b73"
#define GG22D7_4_9 "cb3acd4437a485d524f965cf9555d254b34ca961931d883a5ded4effba235769b4a478ebf2fb7942e4611e478e4280b6da82572b99bf5a41ed"
#define GG22D7_4_10 "4e857f42137e360b2eefd8d749add338ec6c1552a1174a207004131820f195b68c2ab9c4b4776b139e532d3e1f61663cc3c6601a2f4a9135f1"
#define GG22D7_4_11 "1ca13824aef531ee187a761c1231c505cb18eab8f8b21f268cbc76a1112b2a226a95943abbc5a15c98ea7c927425a5159a3cf9a0fb5976502be"
#define GG22D7_5_1 "fa3b1dd42f61ca6855446a70e18a7a538618cacaa37ef00e9d52817843703a2afee2dd06df95d066e0c41e7bd435cd06b51aaa4a29c8a9a9e0"
#define GG22D7_5_2 "9395c77e68f2124f89f5e218e7d4a7bf672f8abb1b3f69b5421fc58fee304b3163d4f540b5ec70777e3331e31c776d3a29ebee4328b24d69a1"
#define GG22D7_5_3 "14c67496f48a2239f595bac1ead67319ac8a5937932e31b6c318309e024c48995bd300d639da15a11a4131d02c8fc229f59933cffa6b62e9777"
#define GG22D7_5_4 "17be32449f372b79db51b81d823866c5dd41f21948ffffab0d9eb4139a0a863bbf91c7c31ac09d983b19efbac03576727742c1af24a833462ea"
#define GG22D7_5_5 "79953ce4b5cd125ec8689b362473cd543984702581a8fe760ec51746ad459e0a6b69e12f57f83c0ec4e02c9eaffbcbdfa6b3005c29d1d42c52"
#define GG22D7_5_6 "1273c6c0cb71b5263267744a83c5e8cfed6f20041624546e70d7d1eae1e2c64539a33d6f29315b41b3f869eaddda05b448b845490813e394c13"
#define GG22D7_5_7 "de9098bc6b15e1e3b328d9315c505887950a5151e58fffa7d3a1fa0270ce5f433847b7bd11acab314bcc0f1232652b326c90d01bd222aefbaf"
#define GG22D7_5_8 "131747bddcab7fbd9c6df0712c6d85f595cf940ccdb67095b70eb66fb37925f2b54bd1131552de17295cc57b1b0dfa36a240198bcb117eeac7"
#define GG22D7_5_9 "1b54a982e6f4b7fb0c08acd128a626f580a05be35ed73b0b26905f1cffc5cc477fe214a70715cc346148f5cb76e00bfd97d5d9c5ec777ba57b3"
#define GG22D7_5_10 "a6f95206e008f21f24dc8da312b4c708a42d99234d730607983b4677bae9e90a3278d75737f655aecf95240562f88ba01fa4b24da07185c1e8"
#define GG22D7_5_11 "136b86b4ea2842335b5473690cfb579460798aa38dd4d9853121979e8db21615fd99e6549a1b71010158bdec98bebbf37aa25559b15f0bf214"
#define GG22D7_6_1 "1ce4c693a676cd931466f8fac322f742b1f33b6cf5b63f5d354df295a416be6c185274ff18b861c4edc2bf85f721c224cb76042b24a83290f8"
#define GG22D7_6_2 "129c912476185709795e1c392a7fe6ad9d98885e61f5314a5eea0e3203bed9d44bbb469905603a06bd1c70a2c3c27cbc52877def49683ee11df"
#define GG22D7_6_3 "14afd379e41d5a265cd4fb7ad7462e8b221bed1abb26dc8595745134f02bb4ad8c349eff0e15f49bdca308d6d4049f1f8bb0b5c686efd05f191"
#define GG22D7_6_4 "9a9bd148a106534d1073b328360b789267e8becf889790e64e4c83d424d91727605d51920a3fdd2564b64ee275026eb267056dad0cacc12385"
#define GG22D7_6_5 "106838be2b2645f7e560499878a84f5e91ec2cc9866b71578dc94e31fabac228627c22fb61c73fd19dbe5ae48e8fb7c47416c5fcb60cedb6770"
#define GG22D7_6_6 "118e3b802838773ce311d4fbec43f0340941e004bbfb60482f0d7b9d30656a59d13083b06931df28775e43f8660f362c3efea5449f68df8001"
#define GG22D7_6_7 "f4ca041f1ee50c25bfe409a5265508caa6af202d4c1444ca4278cf9d07ae3eae145b3197c89a0f3840d6e14c872bac24bdce39149dc7a2f6a5"
#define GG22D7_6_8 "55277bcd811d7f2601774567a968b464a52ac23439108dfcdc2cb94f23044992dd0a060f12458991500f463682ec224e536a3c4f32f3de5526"
#define GG22D7_6_9 "15e9f8f01d12f9e01b6c4aef819b006f0e7709d5c8d44ac6cbb320a0319c14bc8a3c4752eefd82600125d1b645c23c579c5f162544b34d53136"
#define GG22D7_6_10 "3fd6a9306b3dcd2f5c6a3752cec1a345881b9b8f7c64bdf4eadf1741b893ea3de8460029c38cb9750b935cf90b7a4823208268ce0fc0b13d9b"
#define GG22D7_6_11 "164fdb64fb3f934a3e2d1887a3fd3674369ba9119e409e97a20022aa0a4f3e899e17faf050a628be86b2eb800468df3a66b94414345160d8e5a"
#define GG22D7_7_1 "4a4d070620153e741d89a6c326e6ef84ab836752376c6d9ecaa91669d923aef4d22fb76ac5894e6ce1ab9065f9d3646eaf4345761d06157003"
#define GG22D7_7_2 "84a4ed47ec587b5a71b26535e627c19109598d93dc23c7ad9caab693aaea2dcbf1cd9569627a10c8636e73b206404dea2a010e9c1c1aa68d6d"
#define GG22D7_7_3 "108a59c40d3b3b6f84e9ef7ebd5260f1bac060d8ee10c91aed780e807c4abb62fd77976e2680c7588ee5451a63f57cf1c4b8e9f77c388bdc9d3"
#define GG22D7_7_4 "acdf120ecf3d7858ebfa04bd3e157c29d6d018c18b23c2eac897cf1762028d6d777e88f932e5c263ca9e9a7e8ee227938a107024f7af44447f"
#define GG22D7_7_5 "f821b2fc0dc2948c77ebfa07e435d4613ca1fe43e0553ea0ed2ebb596458da17964dbde9926e124174ad032a2bace2d18702dd3d36f61acd48"
#define GG22D7_7_6 "1d6198983fab6a275866ca54478a1cad21a8496a5aa818716df7bb9668b9ffea02f2cafe623fb14f94be5d29e6408e70a0d21d4d6f1c3d93466"
#define GG22D7_7_7 "b9a16313f8277140241ce8b4e5494d26de96a035c93bde7f1476b4bc1d68990b97de3899314968754c0fef2c2791dc699bcca4bc98e0be30df"
#define GG22D7_7_8 "12656d0fb1a02b2021ee9556d3385bfd39d03378146b7c7e0aa8e297a4235ca50aed3103680a6ec057a8637fe2b3a9fb1b927a6b1d90873cfc2"
#define GG22D7_7_9 "1739c039f457ad5efc7af913485a44d950cf1f0dfa1807c7f57c4305e68e693fb4ad8d915efa23eeff5d1a6fc41ec1431e269169e949b33d9c8"
#define GG22D7_7_10 "12f327e2eb0a2b41d60285cfee10d54500636040906263d4cb06b468650f88769c4b24268f0fcfdf729e9aef190ee8902b0125348e3fd4092aa"
#define GG22D7_7_11 "148b50ee41716cba6ed57328f3990b3d4ca1e90fdeb6dacc0e05eda9801e5bc92f15e687e988efca295d32c3d108f4cb2d86dfe8f11e3a179c4"
#define GG22D7_8_1 "b1d4545892dbb6185e29d7cd498733c62d1a2ef900f118cb612278b48ce072f7056d19d69ffa402df2ad31688c5d921ab562fa9c961a92e58"
#define GG22D7_8_2 "d1adc866ce473cc733c71948054c772be2bc68a36154b9550ed8ea09538f5e4c68a14567414fd418bb376cabaac31c20dac6517c50cc7dda2e"
#define GG22D7_8_3 "1202665936c57b92543d4c0cae5347f43a98296c5e3ae717574ff58804894c1b3555fef7e68455df7b72554480efd66fc0ebea628a313d04956"
#define GG22D7_8_4 "3ae2441230cf9ce659e3e8dbe85d73ec054249bb9803aa5aab3af9b141e59bbb52c29731aae8d5fbabdaa9d5914879d7e6059ad09402859776"
#define GG22D7_8_5 "179bf8afe0af1f84f47321a4cf69dcefffcf96ef11864a30dc8123ed72c470790bd357fddd0425f21127fcc65cc87335c6dcc6c6123db391c2"
#define GG22D7_8_6 "5693e47bc22f1acd30117f0d783a6ca642138556dd11cfe481223a5a64f7a475ddcb17ff8967cb4322258ada590bbfe45580ec3766c06633f6"
#define GG22D7_8_7 "15d54471d96085652752134be74f15f6f70116054554387c8b0de7d6d3ae4106ea314f082e0a1bbbb81ebc5b8251e4d7d09dcc269ed922b20d6"
#define GG22D7_8_8 "1a0f7ae77505af023b0e9a128d4ef79e49dc6b174eb9c5be28815549a83e647f018fb5ba6e93bea6d0d436a872698ec4cdb52516c175a72c8bb"
#define GG22D7_8_9 "e8b5babf262a48910c416bc85faaf0905c2f7a8ed13de05f841385ae8d3b1c006027c22a6da753460dd8130324a36f42696f210f6a768ab084"
#define GG22D7_8_10 "1d1cbc3ddf8b7f5b6be9b4ae1dd30881fad16cca5ff244374d45975a59701cfc172631fdde29bc91b6d768cac30bc7e26220d20613b07c8042b"
#define GG22D7_8_11 "371e99445027bad8d32509dc6fe097beea043f51ba6616dbcf9c66f367431cd40d7b5bb4139518b49fad3fb1424776938c1bc5f5bdb0ee11ae"
#define GG22D7_9_1 "2197b1ec5ec0a8b83f113f43ae356fb74d3242356edf810d52c6a6a6c323722d3ac5f4b92e796c846d6aef2eb7c9609922f41b53a720abb624"
#define GG22D7_9_2 "10805882caa33c7f0ca258a8ba26e3c5817252974073a81ceead1cfbd15a7e00366942c2f3ce2df42758ab9401765be3dc6b43a8c83bddaad89"
#define GG22D7_9_3 "15d0091c3dc363cf83cfe8ff58c97ba89f68667d1d607ebd9427c927d2380694fdcc181780c9b2c6119d5cbcaf31cb16a0da52a54e63fc388bf"
#define GG22D7_9_4 "3ccdb466c4396925b49a01f74f118d7b21418771db264b7821a86058e2e78f4573d14702db03bbdfb02590091bb07c21ed5deadf2f2e0cf3e9"
#define GG22D7_9_5 "bff51b1368d90199d04d5f11c413bc00901b1551e3fc2cd70f72b4b573bdd4b045a95419caa66bed6bebadcb7c79184ce0d83026bf838a1698"
#define GG22D7_9_6 "10659a061233670007316868666f938536793ce1547e708ce5a3982fcaa63edc2a929c5062bbec031de381c5f78718c90a7dc9fac7cbd6e6c53"
#define GG22D7_9_7 "72bbafab60cb5a8b2defbbf0356b4738c9b35dbf0ab4b3ad0788327db5f809c79055b0a5d5c640e406142c8881798e0a3e475058b797b669de"
#define GG22D7_9_8 "ee2048abfc6218dafe773d7d513aef9079e8bf36b6051aa97033245eb28910ceeb076dc84bfcb72abc687d956fbb978c6a8c6007255fab1304"
#define GG22D7_9_9 "35f9c0272a8254ed1ee54a949384d5a67055472276b96b9d627799d7be7b66c8f4f2e8ac47e01cc5ade836d5f8f0fe3efd87ca7dd504f3bd09"
#define GG22D7_9_10 "d29c774a2be9761eeab87417a4883f8336517e3523162faf860181989bb4c267bdf84547131b85b895c01db230b72826b0e03799eb9297ef16"
#define GG22D7_9_11 "6f186a50382083e0193d5cec63ff88b886ba7198f470c4cf9206ee165bcceffa6aca964a9ab680acb7b94e69fa331bdabb77c40508792b9797"
#define GG22D7_10_1 "104c3f3cb4f1e8e08c6f58c4524ffee611a7e256a79a5bbbb16d4d7341159855b4642cb4cfa14263c9944d3eb7dc684140e13e4c246fd96b42e"
#define GG22D7_10_2 "9c469623cc75a558bbbb7bc2238f4e022fcee80530e8c844733d2d4272206f330d6a2ca0566c75cc2f6c6d164b204bdfd45031d97c113c2161"
#define GG22D7_10_3 "fa9f25298a19e61a354573c2e8b04e8d3a1cad738c9914935b5df6636e999f164a2eb23bea98a720293f1d525d4cc53c5daa084840dd13af1"
#define GG22D7_10_4 "23d4a32d5c1cfcaa34c12b21a49ab8912627b382a06a46b6d578a55be02aadaed709b378ac6acda32caaa471d14b5823dbd766839d0c681bc2"
#define GG22D7_10_5 "6b513ad5118fda677e4e578cbc8468569c48a9ded5f56330865731e5774a0071462b4ea8fba89309a084d16571290894ff06b38bb1879b9b62"
#define GG22D7_10_6 "1c83756e94784f1aea971e1d3a843a3cfa872d6210c5aea8d8ea04a11f5774585ee3e77e4834097db0bcc784171932669c9bc685030c4c134d8"
#define GG22D7_10_7 "91b551a6edb1d23bcd55302c13f765047fc6b11afdef6b0f69cf6af08bf1234cb61b4a4dbe342f26c560ddba7904effb9bb82a8294dbad2228"
#define GG22D7_10_8 "174aacf77df003c95ce9194da75f1668344c16cc904508c6428f4815a7526d7353798a8240c54a49448060374fc3858b9e929398cf9dd33e668"
#define GG22D7_10_9 "7dbde9b15cd53430e6193e958dd280849fa273ae697ef93096829ed351ae1ae28937d8a3f4351856fadc67da4d6a75a884e066c686b1345a3f"
#define GG22D7_10_10 "175759318230636c7ec6ba71a9ad7bb88db4afff54fa461f7247da17b66ac34abc60a983e3403e75c65591780be6bd09f9a951568db4be26251"
#define GG22D7_10_11 "8354d6e11fd0e18e11d67bec05b3039422363c23bf42008f22842cddb5ca5373fdba3db8b2f591e99d1214cbcb74f5b9df75895f37566170bd"

//For the Im of Fp22
#define GG22D7_1b_1 "1283379a5e5e4fec60aefd7e71c6431194029cd53aeb14210388d1c77e8d51e15ab60e86cf5afe989d3cb22233205647dc19a5b7722d632046"
#define GG22D7_1b_2 "67a81087ac759589b9523bb2bf78feace65caf0d71d62cefdd7dc2d38d7affff33bc7d3474a072f77347458e22c0a7e2dc71a8a48ebfcb373"
#define GG22D7_1b_3 "17cd6ddc3064148d8306a79a8aec0fa105bbe13939dec26b767e812e018de6b56a3bce6b078514f7424cf6a149283bdf0bf9f20cde13378a547"
#define GG22D7_1b_4 "12ee7d568f2862bd2586b89255e91aef667e8ba69f6b50b47676c8d8a6ce594e1eb87e4c0ddb8aff3519ffdcfab22e989f6e43965457e60a311"
#define GG22D7_1b_5 "2802f6ef298c1d3f39c3c99c48913f6b98db98b29cd6c334c7de70b82c480c3367b350a70b22deb006579672255cfd703bf3f7421c798db1aa"
#define GG22D7_1b_6 "a188f790cb511701d5294b89e0f0bc3e3124f68781a3527b847c4f313359e99e34c659310f6d20edf5f9ea801b2b87697dd437cec2a80ff5ad"
#define GG22D7_1b_7 "1732959444cec9e07e3472384d143db317bc1ce7da4d5034ca794df6744a901802fd455b1f407c052a0d6662babf9e9de529e6be0ea1de27051"
#define GG22D7_1b_8 "ac1c9e8f045af09535f5ab5a875deb416644e75b5327365d479cbf90e9078c8d51f9869d7e618c0b4e92baf6ffc7589ecb05f5693d616e99c7"
#define GG22D7_1b_9 "143a77a71d9388279f9e07cb0e259bdc28b2745e257a3df7540de442116a8bc666d0424efae1a436916f2ff75cbead5a8aff542bbf0c8692c10"
#define GG22D7_1b_10 "1bec5c6a972fb3c736f5164785d5601bbd86ca1d6637d26c3eabe0639ddae06e56ad258df26b9e2d9851af48d173790e8b3866ba4700553a803"
#define GG22D7_1b_11 "10353f5cdc46b268377acf206476779a5cddf746c417f0557952d2db8d12930ede5b2a054d3fdbdd1c7892c27244951cf976fa95e7af903f023"
#define GG22D7_2b_1 "1640547581b7473ff1945fe961bf0b4e9ba7b39ffca290570345ef769ae2a5b6ceadbbc3ec4e067fbd4e443a5c186c63865952c4dec28c68cb2"
#define GG22D7_2b_2 "117411026bb0724590fdc7c67e8cf428da703239dcf8db49e3e10ab40bd5453c3a3c609db2ddbb9b5c45bea5bdc7746db06b9e7dd63269a3ef5"
#define GG22D7_2b_3 "8f25ec32c59c22852a9027dc696b1438f67647cee37730abf819db2d88364afcb0af62b069c56e41bd5a77a00e051e8c7dbed555885285aaca"
#define GG22D7_2b_4 "167b925bfb8953a2f42086305a4041314583f8837ec2df17dd7b66db3ac6914ad42290ffdceb2d6ab0ba10a44ea7976b5d1dbf38d154fbac4a"
#define GG22D7_2b_5 "e7b9c1e1b75d9016b41c24b54207fa86321aec9004d0fdd57b006050e046460e6798974bca2f4dc2f78274b34d0862120ed6a43aad7f363ff5"
#define GG22D7_2b_6 "811ae08b14224921e10bed22e016505ba84c5b8a86b5d5116245d56be0f9cdfef06f022428434f1159ad6d587efbe88a226c083ae55b8315a5"
#define GG22D7_2b_7 "4f806ef6d854f95a8e1de3e9087e3e7a64b00847e6bca942082f184905c0b85253d4f9f9f9d56a182ca3e3adac4f296f744abf4834e3ba9e93"
#define GG22D7_2b_8 "12973e744b82d6b8b3bc69b67693b0488c42334c2f224e57ab8f615b24a3279a7a709cd87e27509da9e2f5e30f60a0d8cb2c09fde343c3e5e4d"
#define GG22D7_2b_9 "6bd719034872354521078392f9ac5e63a7fcecb01c0ae6d75a4aa0be5e521c44438ef308eaf0a54f7285559560804a70cc67dad9859ff9b2ac"
#define GG22D7_2b_10 "15698b25afb053a9d9f5a8445e1cc7c64ddd9fd74de5c1f2b12cdd2246186af701c677e6d307bd1f2fa5dd563654245e1143bdcfe365fc2c1f"
#define GG22D7_2b_11 "178c3149a2e5b45af903455fbf97e170d2c8a4c2527d66b2268dd63139fdde3ea1756dc5b2d03dc4c7cc97324e367d6910b1f14d9481e26b6c6"
#define GG22D7_3b_1 "48be339ef30657478f2eb111b944996dfb7a600a142f5ffa275a53f3fc5aaa3c6e18e784a84ec088790598296f42a9b0582a3df8673b8d769a"
#define GG22D7_3b_2 "1226d2eda988624cc9b3c4ae9def224c2af3307922e47b6b5e9f739b8f99a0db4ab13a54bef633fd730b06513639e779dc2eb1859311e421e8c"
#define GG22D7_3b_3 "1af9fffb16c34da52f4357197c1cd1925136453a54c11d2371e528884f8ccba37cef10dbd08fd0d397acf6cc5e34cb82003cc31eb928480e144"
#define GG22D7_3b_4 "1a797c81e227b4fca6345a13d01ed6705c027cb14eff019b78cf0bc08ab2930a872ff8affc21260d13d6d48f74eac24980c82162edf6c67521c"
#define GG22D7_3b_5 "9d21e007934ac02da0d3ebe44059c50bfdb5dae05374abb81861dd7e43763f1f9655a78b4c558e68f230c03642d57cbaf78b5cfd2da860c2dc"
#define GG22D7_3b_6 "88eea84d03d53fe6dd236c18fc13ee31d326f47127061a3fece0c5224549bc783ec10fb9d588a4dda98b74ceb3eb8a0e3e035006421e2b70f5"
#define GG22D7_3b_7 "db3e5f6932b746221af3d47d547f8e384045ee665830c76bbb0207410c96a32ae189eaf14499daf720c14e51296f1a2e9434bff03a0539a73d"
#define GG22D7_3b_8 "1b58c50b87476742c935ed5d8e163098baa4a3d1eb590b4847eae0a56b44b187c62683232201bab0f7aea664cf9280d4e76e4fa6aa2fdf2e292"
#define GG22D7_3b_9 "994c55ed0d69798960f678ad2f6e73967b15adb28388f24944422d88793c18e3119a751c02c3213476413e9e75860a27d7f85ec307dd5e582b"
#define GG22D7_3b_10 "fb574dc7d9d66f09a03b6f44444b3c337d5d8cade2652620924f8df5b8b6f4f68b238e260130fe25def1c7c8a48ef94da7e36f39e353063517"
#define GG22D7_3b_11 "1c32c8e5a6a2fdda20546762fa2590d93b85cbe67c646efdaaee95118d40a27c9a06d073d815532301aecd1c67c687c76ce473e4c80fd8b3621"
#define GG22D7_4b_1 "c1ce63de22c4b0d30684ac09937bb8a7923117fcb539e40600ee4949f1e5365de6da041e2f27ab1a737bc65f5afd48e83448a3e0d4d20c790a"
#define GG22D7_4b_2 "d04e7687e52f32bd4edd16c95283c0918e20c25e5033635316049975284cea6f01fd52882ebb8e37c72b4cbe6f94b06b91fbf74776a84e3963"
#define GG22D7_4b_3 "a9b91e014b91f65d4858060290c744f0773aa8da4976c142b07c1d9d59561b286d57a8cb1f0154bc5a00e7d3b5090239e501f3fccab6253b3f"
#define GG22D7_4b_4 "16205baeea43f19da0424f7e136746e1049be08232b60c994bd6ff6a928b537cc1ad11cb85fddfd9919a5df4d1baea6de6f7b19eed3d21f9332"
#define GG22D7_4b_5 "ce47cf8a844c90b50b9713f7c251d363b1e4cc69748c1e669ae5ad699a976946798792825964b849bcd4e435db3b0d13ff9f6ada5eb89befc7"
#define GG22D7_4b_6 "e8f488db8c770f06e9721fb2ea1bdb90340e8322220f8a0a1affc8cb587cf1213d55fdd0e4079ad82203f711bebfd30feb078d93dcbc07b489"
#define GG22D7_4b_7 "147ff976013fc9f081632726f0b1a7d6bb875394e28dad8ecdda6c220acf7e7af61343bc9dbc159ccb2587f5160672504db949279ab27386018"
#define GG22D7_4b_8 "114a8fccf738bb15ec7f91d3b1bfa6fce7e50e52c8c93362528f3ee46baac7f59f09beeae059f97e2e8ba1eee563ac0b5011b4a3024dba8515a"
#define GG22D7_4b_9 "5fce2478b50fbde152a1f9c6c8f495c3b8c22c6e6d9cf26d3d6aa8c16500409a4413cedef5850303ff97c330bc34165a579bd9e219be9cb6d2"
#define GG22D7_4b_10 "10dde2bf558911c09d3c6841b4489121ffd2e2193f42e08965d47e7d435f8121f71a6f566059c9eed9325ffb0797ae581832c8427435923ffbd"
#define GG22D7_4b_11 "143c6e5e9c270e9ccf67f3ad3fbd650c3f9d31fdfb57be1a826ddd42ca17652b8d32f59e90e5b0e96c8d15a5fcec0b069b2af0bd17f42653f94"
#define GG22D7_5b_1 "161883450c35e1c7b717770c3cdaa8799dcb1f6e7123a544ee6997ebfa2047209ed2fb8cadec00e7f6062b31ba19171db26609def1e4a44dc67"
#define GG22D7_5b_2 "173a65e9e980dc0ba2d16922871802284c1d0c70e16b2631f3b5362f3d1a76de156c0dedbebddacf3201356dbab7c37404f7534b08dbaaf0fa4"
#define GG22D7_5b_3 "a3f28c0840240748237a52234e718844f8380322d82a5b89c9435fffe4288957cfe2a00d8b065660f3d2d80de3bc7ecf1d168bf0e2b2219ccb"
#define GG22D7_5b_4 "5d8892e2fac65c5cb429dea6fcfae539c629f2c57c5df51eccb286b21cec9b0e5a3b3cf3831ed36b954b2522b5eba5b7d4af82fa52b718515b"
#define GG22D7_5b_5 "1cf5bb56b3219b7170161fb548b5ee38e6a594b30d67fc65f987d619b0e742f9744c7c3071300d359082b6a301fc68f4cf08c210a2b06cb5637"
#define GG22D7_5b_6 "f9e60c2c0e07f8fa7e810d48b3899e5733f72c60738d5583b556325c8a2ebc8470bf91ea08151c1a7ebd2c92589b59937c46e2bd21565b8c64"
#define GG22D7_5b_7 "58eedd3ad07abbe1e32b68aa28fea80d23854f2efd3ed1a88154e41c20aa219cb0ecbe1823e46542f68e67033f9056070928263a93c7593003"
#define GG22D7_5b_8 "19e79f6a1f39a7f076713aef5381fddab82f1a6994ef664b095d1763cbcbb715e87f16640b5cd75361b61845628f8ee99db880757861d36dd42"
#define GG22D7_5b_9 "1afa3dea263110c7bd794a3849ab7a095c091e46a7b9cf22df667ceed273448e24b05ad2e91d58b2e849bab1200925804b8b576a86c70264848"
#define GG22D7_5b_10 "11a447b02221fa83e21792d240bba781afdc297c6c52d001a05dc818e494ee603e37b35daf24504eb6bea4f2dd67aa1ca7ac9c04947d1b4d20c"
#define GG22D7_5b_11 "1a0724e07556e4fbf1abde997ba80ee7f3bc369d2b1d14625b06fa11ccd27f8d5873fc4fd3bc09446cbe229f39da9b89c601ef8ab077ce63df7"
#define GG22D7_6b_1 "1d2bdecf3a1bba455f754d2272a3621629dbcc5e04edc35dda7d77b571644e6a85020559a11ae9cf8280b4d44c2284d14a81cbae6857b9bc906"
#define GG22D7_6b_2 "1bbe03587b8fdb658a396ddc7cdf12c911ccb585db18b18a063fcd8fb09f1ce18447e379a171a7bb619085fdb42e14b12ea35c3abddcc2b398c"
#define GG22D7_6b_3 "1817f21a9f82b5f79e67534d30bfeca980addf74ea99fefbe49d7741c3b64451c2b00c11d983a92f1b95b6300abf1da771b9caebf10d919582d"
#define GG22D7_6b_4 "46749ad8a0187a3e6f1623f2ae541bfb887fbc14476111c0f5042f16f3cc1eba90de0762a2b24bfaeac40c50e187e9e031f511a6781cd76257"
#define GG22D7_6b_5 "1c446d941f312592478d74799faebcc52fa3af639803a13bbaaf7ff196ba968f77fdd7649e5a328f45bc928add662e0bb75e1729e850d537fcf"
#define GG22D7_6b_6 "e3b7c0dab328e97b991fdc328a732ae1b65ee76b4feb0cea52e11a68d2315f0afe3aa3d50d65df102ce3e2a36e157267769ede14cbcb7d0140"
#define GG22D7_6b_7 "544c7fb68ade4dddccd1a9bdf88a1bc7169b62a1ba754e7b725e559ad146ce84080907edda43c73be49296d3a73c973fecaa437ef16e2f479"
#define GG22D7_6b_8 "10e15c9698a236288f5d19ad89436376c7976dd87cca8a44fae6ad6c1bc845b598a858ef38b681faa91d21c2005f5beda791883a7c379ecbae4"
#define GG22D7_6b_9 "1b12ba5a3c2e4f26bbe5c7e801bcd532f840164455159b7e5009011a50d053b8a798b6facfd038bb9c87e896276ac5ce4950bd0d9bce3b8ae28"
#define GG22D7_6b_10 "650ffd46968682fd1a213acac326a87343e69eed8b01e8136fe1cc6e05ee04ba4f245f4b2dce6f031ba51e629e2ecaeaa1eaecd9406a560b7e"
#define GG22D7_6b_11 "186c4fd9a4835a13073ed852b7c91029831c06c54e98b5f512d94c204f3ae170ac1b4f560ebce90b4f022ba32c75bc94cac772a5de55df834e8"
#define GG22D7_7b_1 "fdb99e583ccd33c510ebf8fcc42e7f49b4490fc354a86fd4bc54b0e0b9baa39ca28799a256abddec6c34d5e4bee184f30c4f0a0c7018fffdf7"
#define GG22D7_7b_2 "3731c36caa4a2b1c889d352d0cc488bae314e4f957dec18bcb693b95dc4c470ee1c9e598d753176ca5b9e3267edc1ec391635da2a83e75b53d"
#define GG22D7_7b_3 "9560519c34e86cd8410756315927a1d495fc5c392c7f837f77317647c39743f636bf43da2fe5e5f230a92ad3e9266db860b9a4062241f53816"
#define GG22D7_7b_4 "29bf1af64bf0068241fd26e0511e5f9103012420fe49764e7f2e7c5536663ce9a0819cf243e8507e3ff6c6e9be7729ca76b10c7509825eadaa"
#define GG22D7_7b_5 "3b51eb440dad0b18194a7de4cd0cf8a0fc8f89b88b7259fec12fc5f10f4959aca336291fb1b53c1fac84952801cc12019932fbda45debef57d"
#define GG22D7_7b_6 "1216b5c02cee80e50dcf81c4b6ba1389f7932b3be16396cd1ab8f9955f04122adaf3ebb74a9687ac3b7e3ac6210cd1b0fafbc17a1f0ded9acd4"
#define GG22D7_7b_7 "e38e4459113daab350df238070de1ada95524d9c5b3dd56c2a58f6e48b1818a5a199c5d4cbcfdd62fc656c90acc6154b1b02be7455c985871c"
#define GG22D7_7b_8 "131b00e5ce5058bad404bed33a3f37fd31627d45c8bd2158ad22a7a160de65fe3c6d92e8d336a95999f32fcca5c152cf052b264c50d88baeeb4"
#define GG22D7_7b_9 "3f891a487b7cab2d5d7047f87ad71bbf79151d7c97242a164d50b6d653eff1fc3bc895eef3fb86b4f1910ae34ef9fac4a7b6473caaa7363678"
#define GG22D7_7b_10 "157e6ecd4ee32d7448a587019a9ce8e86cb4a3e29f7dcb5cdf01f28face8506aac3a512dba6d708c59c0ebd7adf717bc8a6184a65eb49d04402"
#define GG22D7_7b_11 "13cc076ab0d46f5eaef46945d5a7da35e651d2212fbcb2d2e6d9d30db54242ade89d916b392a71cdff8db279ee08764c7a8b213ca2cf24c3926"
#define GG22D7_8b_1 "8eada7e1f5f747c62c49292da779e6ed7b392b06343b5645bc7e3b04499b5d6a47851c0731c4a14c5ab396608c5a39fbc68200652b4dca6b9f"
#define GG22D7_8b_2 "d7f8d763d71ee2c5e823fc4ddeda2e2db79852bf43239152d4fbc5f379e5b1fbbb1cbd8a39669771affa18b6262cf7450011f4f7e7d3ae704"
#define GG22D7_8b_3 "11527a0d819591b7fcf2400b81004819f35b3a1c4c238daca4a224217fa6911998911dc472fe5e5c11ab515c5b8898fef210634dd5c93ef06fd"
#define GG22D7_8b_4 "e0ca0ef07804e039d7e2ddb78977dee1ad4aa36217f2a99f1043fc8fc03793fc0f49fc5921c07ff3f9c1d50c80d55d7797f64c3f445f612971"
#define GG22D7_8b_5 "1ab9aa985cb10f42f2ee6c2a1bba2223b3d6bcbc71f62705e00c03c2dd58e04eead8a71084ae1ffa94a442450cd863bba64106fab4ae0196c52"
#define GG22D7_8b_6 "17562b08c03adc6886157b88046859438ad2e18b8a843e34ddfe4bd731e40ed42d1c9a868f34852eb5e0377aef18188ee9fff76f74cc5ea9140"
#define GG22D7_8b_7 "15c564f573c8b56bd535560f1c4564f6fcf9251aef4cb729a3a191cde51996fda71cada31739a946971734dc0150e8737a9259c7a95dc4924f7"
#define GG22D7_8b_8 "97a2cdc5fa9c17081429ef82eb686ee726b6be474d0ba7c1c70aab35a5895bb2ffd2da72a44013b60fb755cd6fe2dad17165912302c8560de8"
#define GG22D7_8b_9 "8f0d5c2a51e5730fee41ce67d2bf708f8b685e3b903fbb90f591823d80cb12937e6b00bac459c4385dd9b399083661ff6b80380d1aa8070899"
#define GG22D7_8b_10 "1a55518ff47fb02e6079ba5416d05d9e438eee49815d0d41a6aeb1ff2ef43f172956720f3dca52c0fee16451fc2a33c87819fba014b8bcbc22f"
#define GG22D7_8b_11 "3232e48294ef2e5b6ab235e82deb2fb6f669f45052e8a7a1396b67d21fac0ea1ff8d1f8e75044617383df86509a13a10efb98783c391f70865"
#define GG22D7_9b_1 "daf9e06ca2e73c92b0494b1de7a59a54c6b0661140ba2be3d1b2575846acbd0d5df01bfbb5b4a6448d930e5975fcab44a45c0b359d9ee4b072"
#define GG22D7_9b_2 "f17b032282a8fbd1659dbbb5c0133339d213d59a450edbb5776dd5635de2ffef05ae04b12fa2ef601ae502d197d18d3318e2037a70e59fbdb5"
#define GG22D7_9b_3 "12b321c2880e8681d7d64c5c4ec89fe028d4c600366111af0b73e4b4aff6aaa92f651b2105fd48f423585b2d7d73f6be88243b9b689ed3e9293"
#define GG22D7_9b_4 "1ca5f051a29a612857c24b83852a27bc47d4420d9560068152a38a6aa7475a2402b34c1ba9085f71fa0e12c05b2b6140e7136dd7aa7b0f2afaf"
#define GG22D7_9b_5 "1a26d23811fbf00158584cc83e3c21a80cef582d483455d374e02dc023b34133aed32cc2900240640e4c0fba785f16aba8002d0297687a3253f"
#define GG22D7_9b_6 "ad81cd9b1d22fcf985e02c3e520932d717807620d03a9ef6e8a2eff37556d673373bf4d25b57ccec74a4da9a4a2908322eaccec6c22d6517cc"
#define GG22D7_9b_7 "118a18597ecaf6d84efff5b086f1ac8a53c799788dd4fb8d53863c882fa922eeb433440370c2428a84dbc224272a074f77d538dbdb272802465"
#define GG22D7_9b_8 "57d5bdeb51578d07052911ab0ff5b0e10bec75a0201e365d4887ad2224423176e905033f4a5b8808e0c925a884c5052b9a107253c5087cd6be"
#define GG22D7_9b_9 "180319d281413edad298df99c1f23269470010499b4a8d936220bb0ea5b434bdc05af3cbfdda1f1d8ede95a8292fe0bc0c20dac275f48194c3f"
#define GG22D7_9b_10 "19720e2ff6e77169cad996b04a60713da10e4bd27350ce08defabf6e65bfe87ee817cf43e1de38491f4ae95233b109fbee20a67aad94e814304"
#define GG22D7_9b_11 "1b02828eddd8dd0c9a442d151a61c331b06dce1072ecc5ff782aa9362ddbc5f08d0450ae735c02043660c0cefe42cfc0f3b7f5ea49d1ddf833c"
#define GG22D7_10b_1 "9ccce83b3e17f9c8e217d550bdcc9bafb6fa5425e3399ea3f9d228fcb54cae609025e03e3d7c43b5c2aa32499666793d94b10d96e661a57c16"
#define GG22D7_10b_2 "1719eeccbf00d325c532d2a4e201fbae1da1f159f68fe2512fc81564c00d095bcfd5c1aa446ca07a161cb59ba88e5fe86a248e21a6514d63439"
#define GG22D7_10b_3 "1162da4f928f1290558e8ee3ac59bede836c222843266864188387c714c3ee2f69c722f4c6fb54140acd1568ba8caa91c1cf07481c0394f0bd1"
#define GG22D7_10b_4 "1958c13e3ad7ec0d712cc16219aef1866e2248add535019ad6079f92ebd9bf6e8948eec76e7b0897f54af00cd1a695b1bddfb40da58c66748a"
#define GG22D7_10b_5 "1af75418cc0b586f1b163ac3441c21b9ac6fc4522354c3a80d10e79d5b91148d7f5891331b1dd28e06d94bf8f5ed796d2989a1f8fe533087309"
#define GG22D7_10b_6 "1212d194d427a5dc66c74726af594f89ede659166a7e6a4d7f4e64fdd0231807798e23939f8d82b7f7217ccf78c8d1a93b709c50f22963943f2"
#define GG22D7_10b_7 "12a3f982f799d2ab6ff9fffd02e405f3f0a826ff517087c4be940fdcb554a5d3f20e4bc1a13d59fbf6383df0ae369ac56d82691fb8c1be0b4d4"
#define GG22D7_10b_8 "18cd20d948323908fbd989e7e0a92557054a013c9e7ec00ea285399043f511a23cb8e7c1a7c81e5e4a47a1a94200cc32d4375743c206c38904"
#define GG22D7_10b_9 "2eeeafb791870a83e07ff67e5bc958e849a4778d67ab876aea475284f17ca698c838d4193f50eaea7560ee19df3208482f647f82f5c252cd09"
#define GG22D7_10b_10 "1c0d006e0dec2fa146b44171308cbdb98b357c85eecec4fc8cec216c33d579daf12aedd655c09cde2224beb67aefc21fc39a60288a9fd1f4846"
#define GG22D7_10b_11 "523393f54e0ba8e09dbd25174ee3fa604500045d59a86482dcb2a7c54db093184dbacd62111fbb86c0437e8cc80ff04b14b99ba0d61636149a"

//(p+1)/2
#define GG22D7_half	"EB6B7FC15A0289AD8BA7CFBC493B9136D1E01E0D03D472A0A1D8488FAC7941AAFBCAA9E159D9A67BE554722C0DBCF0F6F02D143C5C70979F3E"

#define ASSIGNFP11(FP11)												\
	RLC_GET(str,FP11##_1, sizeof(FP11##_1));								\
	fp_read_str(c[0], str, strlen(FP11##_1),16);\
	RLC_GET(str,FP11##_2, sizeof(FP11##_2));								\
    fp_read_str(c[1], str,strlen(FP11##_2),16);\
	RLC_GET(str,FP11##_3, sizeof(FP11##_3));								\
	fp_read_str(c[2], str, strlen(FP11##_3),16);\
	RLC_GET(str,FP11##_4, sizeof(FP11##_4));	\
	fp_read_str(c[3], str, strlen(FP11##_4),16);\
	RLC_GET(str,FP11##_5, sizeof(FP11##_5));								\
	fp_read_str(c[4], str, strlen(FP11##_5),16);\
	RLC_GET(str,FP11##_6, sizeof(FP11##_6));								\
	fp_read_str(c[5], str, strlen(FP11##_6),16);\
	RLC_GET(str,FP11##_7, sizeof(FP11##_7));						\
	fp_read_str(c[6], str, strlen(FP11##_7),16);\
	RLC_GET(str,FP11##_8, sizeof(FP11##_8));								\
	fp_read_str(c[7], str, strlen(FP11##_8),16);\
	RLC_GET(str,FP11##_9, sizeof(FP11##_9));								\
	fp_read_str(c[8], str, strlen(FP11##_9),16);\
	RLC_GET(str,FP11##_10, sizeof(FP11##_10));								\
	fp_read_str(c[9], str, strlen(FP11##_10),16);\
	RLC_GET(str,FP11##_11, sizeof(FP11##_11));								\
	fp_read_str(c[10], str, strlen(FP11##_11),16);\

void fp11_field_init() {
	char str[2 * RLC_FP_BYTES + 1];
	fp11_t c;	//temp for d[]
	fp11_t u[8]; //(2u+2)^(index+1)
	fp11_t d[10];
	ctx_t *ctx = core_get();
	int jprime[10] = {3, 9, 5, 4, 1, 3, 9, 5, 4, 1}
	fp11_null(c);
	for (int i=0;i<8;i++) fp11_null(u[i]);
	for (int i=0;i<10;i++) fp11_null(d[i]);
	RLC_TRY{
		fp11_new(c);
		for (int i=0;i<8;i++) fp11_new(u[i]);
		for (int i=0;i<10;i++) fp11_new(d[i]);
		RLC_GET(str,GG22D7_half, sizeof(GG22D7_half));
		fp_read_str(ctx->fp11_half, str, strlen(GG22D7_half),16);

		fp11_set_dig(u[0],2);
		fp_copy(u[0][1],u[0][0]);
		for(int i=1;i<8;i++) fp11_mul(u[i],u[i-1],u[0]);
		ASSIGNFP11(GG22D7_1);
		fp11_copy(d[0], c);
		ASSIGNFP11(GG22D7_2);
		fp11_copy(d[1], c);
		ASSIGNFP11(GG22D7_3);
		fp11_copy(d[2], c);
		ASSIGNFP11(GG22D7_4);
		fp11_copy(d[3], c);
		ASSIGNFP11(GG22D7_5);
		fp11_copy(d[4], c);
		ASSIGNFP11(GG22D7_6);
		fp11_copy(d[5], c);
		ASSIGNFP11(GG22D7_7);
		fp11_copy(d[6], c);
		ASSIGNFP11(GG22D7_8);
		fp11_copy(d[7], c);
		ASSIGNFP11(GG22D7_9);
		fp11_copy(d[8], c);
		ASSIGNFP11(GG22D7_10);
		fp11_copy(d[9], c);

		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[0]*(i+1)/11;
			int e2=jprime[0]*(i+1)-e1*11;
			fp11_mul(c,d[0],c);		//d[0]^i
			fp11_copy(ctx->frb11_1[i],c);
			if (e1) fp11_mul(ctx->frb11_1[i],ctx->frb11_1[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_1[i],ctx->frb11_1[i]);
			}
		}
		
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[1]*(i+1)/11;
			int e2=jprime[1]*(i+1)-e1*11;
			fp11_mul(c,d[1],c);		//d[1]^i
			fp11_copy(ctx->frb11_2[i],c);
			if (e1) fp11_mul(ctx->frb11_2[i],ctx->frb11_2[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_2[i],ctx->frb11_2[i]);
			}
		}

		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[2]*(i+1)/11;
			int e2=jprime[2]*(i+1)-e1*11;
			fp11_mul(c,d[2],c);    //d[2]^i
			fp11_copy(ctx->frb11_3[i],c);
			if (e1) fp11_mul(ctx->frb11_3[i],ctx->frb11_3[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_3[i],ctx->frb11_3[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[3]*(i+1)/11;
			int e2=jprime[3]*(i+1)-e1*11;
			fp11_mul(c,d[3],c);    //d[3]^i
			fp11_copy(ctx->frb11_4[i],c);
			if (e1) fp11_mul(ctx->frb11_4[i],ctx->frb11_4[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_4[i],ctx->frb11_4[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[4]*(i+1)/11;
			int e2=jprime[4]*(i+1)-e1*11;
			fp11_mul(c,d[4],c);    //d[4]^i
			fp11_copy(ctx->frb11_5[i],c);
			if (e1) fp11_mul(ctx->frb11_5[i],ctx->frb11_5[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_5[i],ctx->frb11_5[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[5]*(i+1)/11;
			int e2=jprime[5]*(i+1)-e1*11;
			fp11_mul(c,d[5],c);    //d[5]^i
			fp11_copy(ctx->frb11_6[i],c);
			if (e1) fp11_mul(ctx->frb11_6[i],ctx->frb11_6[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_6[i],ctx->frb11_6[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[6]*(i+1)/11;
			int e2=jprime[6]*(i+1)-e1*11;
			fp11_mul(c,d[6],c);    //d[6]^i
			fp11_copy(ctx->frb11_7[i],c);
			if (e1) fp11_mul(ctx->frb11_7[i],ctx->frb11_7[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_7[i],ctx->frb11_7[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[7]*(i+1)/11;
			int e2=jprime[7]*(i+1)-e1*11;
			fp11_mul(c,d[7],c);    //d[7]^i
			fp11_copy(ctx->frb11_8[i],c);
			if (e1) fp11_mul(ctx->frb11_8[i],ctx->frb11_8[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_8[i],ctx->frb11_8[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[8]*(i+1)/11;
			int e2=jprime[8]*(i+1)-e1*11;
			fp11_mul(c,d[8],c);    //d[8]^i
			fp11_copy(ctx->frb11_9[i],c);
			if (e1) fp11_mul(ctx->frb11_9[i],ctx->frb11_9[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_9[i],ctx->frb11_9[i]);
			}
		}
		fp11_set_dig(c,1);
		for(int i=0;i<10;i++){
			int e1=jprime[9]*(i+1)/11;
			int e2=jprime[9]*(i+1)-e1*11;
			fp11_mul(c,d[9],c);    //d[9]^i
			fp11_copy(ctx->frb11_10[i],c);
			if (e1) fp11_mul(ctx->frb11_10[i],ctx->frb11_10[i],u[e1-1]);
			for(int ii=0;ii<e2;ii++){
				fp11_mul_art(ctx->frb11_10[i],ctx->frb11_10[i]);
			}
		}

		ASSIGNFP11(GG22D7_1b);
		fp11_copy(ctx->frb11_1b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_1b[i],ctx->frb11_1[i-1],c);
		ASSIGNFP11(GG22D7_2b);
		fp11_copy(ctx->frb11_2b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_2b[i],ctx->frb11_2[i-1],c);
		ASSIGNFP11(GG22D7_3b);
		fp11_copy(ctx->frb11_3b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_3b[i],ctx->frb11_3[i-1],c);
		ASSIGNFP11(GG22D7_4b);
		fp11_copy(ctx->frb11_4b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_4b[i],ctx->frb11_4[i-1],c);
		ASSIGNFP11(GG22D7_5b);
		fp11_copy(ctx->frb11_5b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_5b[i],ctx->frb11_5[i-1],c);
		ASSIGNFP11(GG22D7_6b);
		fp11_copy(ctx->frb11_6b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_6b[i],ctx->frb11_6[i-1],c);
		ASSIGNFP11(GG22D7_7b);
		fp11_copy(ctx->frb11_7b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_7b[i],ctx->frb11_7[i-1],c);
		ASSIGNFP11(GG22D7_8b);
		fp11_copy(ctx->frb11_8b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_8b[i],ctx->frb11_8[i-1],c);
		ASSIGNFP11(GG22D7_9b);
		fp11_copy(ctx->frb11_9b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_9b[i],ctx->frb11_9[i-1],c);
		ASSIGNFP11(GG22D7_10b);
		fp11_copy(ctx->frb11_10b[0],c);
		for(int i=1;i<11;i++) fp11_mul(ctx->frb11_10b[i],ctx->frb11_10[i-1],c);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(c);
		for (int i=0;i<8;i++) fp11_free(u[i]);
		for (int i=0;i<10;i++) fp11_free(d[i]);
	}
}