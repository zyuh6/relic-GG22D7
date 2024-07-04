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
 * Implementation of frobenius action in extensions defined over prime fields.
 *
 * @ingroup fpx
 */

#include "relic_core.h"
#include "relic_fpx_low.h"
#include "relic_fp_low.h"

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void fp2_frb(fp2_t c, const fp2_t a, int i) {
	switch (i % 2) {
		case 0:
			fp2_copy(c, a);
			break;
		case 1:
			/* (a_0 + a_1 * u)^p = a_0 - a_1 * u. */
			fp_copy(c[0], a[0]);
			fp_neg(c[1], a[1]);
			break;
	}
}

void fp3_frb(fp3_t c, const fp3_t a, int i) {
	fp3_copy(c, a);
	switch (i % 3) {
		case 1:
			fp3_mul_frb(c, c, 0, 1);
			break;
		case 2:
			fp3_mul_frb(c, c, 0, 2);
			break;
	}
}

void fp4_frb(fp4_t c, const fp4_t a, int i) {
	/* Cost of a single multiplication in Fp^2 per Frobenius. */
	fp4_copy(c, a);
	for (; i % 4 > 0; i--) {
		fp2_frb(c[0], c[0], 1);
		fp2_frb(c[1], c[1], 1);
		if (fp_prime_get_mod18() % 3 == 1) {
			fp2_mul_frb(c[1], c[1], 1, 3);
		} else {
			fp2_mul_frb(c[1], c[1], 2, 1);
			fp2_mul_frb(c[1], c[1], 2, 1);
		}
	}
}

void fp6_frb(fp6_t c, const fp6_t a, int i) {
	/* Cost of two multiplication in Fp^2 per Frobenius. */
	fp6_copy(c, a);
	for (; i % 6 > 0; i--) {
		fp2_frb(c[0], c[0], 1);
		fp2_frb(c[1], c[1], 1);
		fp2_frb(c[2], c[2], 1);
		fp2_mul_frb(c[1], c[1], 1, 2);
		fp2_mul_frb(c[2], c[2], 1, 4);
	}
}

void fp8_frb(fp8_t c, const fp8_t a, int i) {
	/* Cost of four multiplication in Fp^2 per Frobenius. */
	fp8_copy(c, a);
	for (; i % 8 > 0; i--) {
		fp4_frb(c[0], c[0], 1);
		fp4_frb(c[1], c[1], 1);
		fp2_mul_frb(c[1][0], c[1][0], 2, 1);
		fp2_mul_frb(c[1][1], c[1][1], 2, 1);
		if (fp_prime_get_mod8() % 4 != 1) {
			fp4_mul_art(c[1], c[1]);
		}
	}
}

void fp9_frb(fp9_t c, const fp9_t a, int i) {
	/* Cost of two multiplication in Fp^3 per Frobenius. */
	fp9_copy(c, a);
	for (; i % 9 > 0; i--) {
		fp3_frb(c[0], c[0], 1);
		fp3_frb(c[1], c[1], 1);
		fp3_frb(c[2], c[2], 1);
		fp3_mul_frb(c[1], c[1], 1, 2);
		fp3_mul_frb(c[2], c[2], 1, 4);
	}
}

void fp11_frb(fp11_t c, const fp11_t a, int j) {
	ctx_t *ctx = core_get();
	fp_t b;
	dv11_t t0;
	dv11_t t1;
	fp_null(b);
	dv11_null(t0);
	dv11_null(t1);
	RLC_TRY{
		fp_new(b);
		dv11_new(t0);
		dv11_new(t1);
		for(int i=0; i<11;i++)dv_zero(t0[i],2*RLC_FP_DIGS);
		switch(j%11){
			case 0:
				fp11_copy(c,a);
				break;
			case 1:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_1[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 2:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_2[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 3:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_3[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 4:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_4[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 5:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_5[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 6:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_6[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 7:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_7[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 8:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_8[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 9:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_9[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
			case 10:
				for(int i=0;i<10;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i+1],ctx->frb11_10[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp_copy(b,a[0]);
				fp11_rdc(c,t0);
				fp_add(c[0],b,c[0]);
				break;
		}
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		dv11_free(t0);
		dv11_free(t1);
	}
}

void fp11_frb2(fp11_t c, const fp11_t a, int j) {
	ctx_t *ctx = core_get();
	dv11_t t0;
	dv11_t t1;
	dv11_null(t0);
	dv11_null(t1);
	RLC_TRY{
		dv11_new(t0);
		dv11_new(t1);
		for(int i=0; i<11;i++)dv_zero(t0[i],2*RLC_FP_DIGS);
		switch(j%11){
			case 0:
				fp11_copy(c,a);
				break;
			case 1:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_1b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 2:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_2b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 3:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_3b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 4:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_4b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 5:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_5b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 6:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_6b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 7:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_7b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 8:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_8b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 9:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_9b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
			case 10:
				for(int i=0;i<11;i++){
					for(int k=0;k<11;k++)fp_muln_low(t1[k],a[i],ctx->frb11_10b[i][k]);
					fp11_addd_low(t0,t0,t1);
				}
				fp11_rdc(c,t0);
				break;
		}
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		dv11_free(t0);
		dv11_free(t1);
	}
}


void fp12_frb(fp12_t c, const fp12_t a, int i) {
	/* Cost of five multiplication in Fp^2 per Frobenius. */
	fp12_copy(c, a);
	for (; i % 12 > 0; i--) {
		fp6_frb(c[0], c[0], 1);
		fp2_frb(c[1][0], c[1][0], 1);
		fp2_frb(c[1][1], c[1][1], 1);
		fp2_frb(c[1][2], c[1][2], 1);
		fp2_mul_frb(c[1][0], c[1][0], 1, 1);
		fp2_mul_frb(c[1][1], c[1][1], 1, 3);
		fp2_mul_frb(c[1][2], c[1][2], 1, 5);
	}
}

void fp16_frb(fp16_t c, const fp16_t a, int i) {
	/* Cost of four multiplication in Fp^2 per Frobenius. */
	fp16_copy(c, a);
	for (; i % 8 > 0; i--) {
		fp8_frb(c[0], c[0], 1);
		fp8_frb(c[1], c[1], 1);
		fp2_mul_frb(c[1][0][0], c[1][0][0], 2, 2);
		fp2_mul_frb(c[1][0][1], c[1][0][1], 2, 2);
		fp2_mul_frb(c[1][1][0], c[1][1][0], 2, 2);
		fp2_mul_frb(c[1][1][1], c[1][1][1], 2, 2);
		if (fp_prime_get_mod8() % 4 != 1) {
			fp8_mul_art(c[1], c[1]);
		}
		if (fp_prime_get_mod8() == 5) {
			fp4_mul_art(c[1][0], c[1][0]);
			fp4_mul_art(c[1][1], c[1][1]);
		}
	}
}

void fp18_frb(fp18_t c, const fp18_t a, int i) {
	/* Cost of five multiplication in Fp^3 per Frobenius. */
	fp18_copy(c, a);
	for (; i % 18 > 0; i--) {
		fp9_frb(c[0], c[0], 1);
		fp3_frb(c[1][0], c[1][0], 1);
		fp3_frb(c[1][1], c[1][1], 1);
		fp3_frb(c[1][2], c[1][2], 1);
		fp3_mul_frb(c[1][0], c[1][0], 1, 1);
		fp3_mul_frb(c[1][1], c[1][1], 1, 3);
		fp3_mul_frb(c[1][2], c[1][2], 1, 5);
	}
}

void fp22_frb(fp22_t c, const fp22_t a, int j){
	fp11_t c1;
	fp11_null(c1);
	RLC_TRY{
		fp11_new(c1);
		fp11_frb(c[0], a[0],j);
		fp11_frb2(c[1], a[1],j);
		fp11_neg(c1,c[1]);
		if (j%22-11>=0) fp11_copy(c[1],c1);
		else fp11_copy(c[1],c[1]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(c1);
	}

}

void fp24_frb(fp24_t c, const fp24_t a, int i) {
	/* Cost of 20 multiplication in Fp^2 per Frobenius. */
	fp24_copy(c, a);
	for (; i % 24 > 0; i--) {
		fp8_frb(c[0], c[0], 1);
		fp8_frb(c[1], c[1], 1);
		fp8_frb(c[2], c[2], 1);
		for (int j = 0; j < 2; j++) {
			for (int l = 0; l < 2; l++) {
				fp2_mul_frb(c[1][j][l], c[1][j][l], 2, 3);
				fp2_mul_frb(c[2][j][l], c[2][j][l], 1, 1);
			}
			if ((fp_prime_get_mod8() % 4) == 3) {
				fp4_mul_art(c[1][j], c[1][j]);
			}
		}
	}
}

void fp48_frb(fp48_t c, const fp48_t a, int i) {
	/* Cost of 52 multiplication in Fp^2 per Frobenius. */
	fp48_copy(c, a);
	for (; i % 48 > 0; i--) {
		fp24_frb(c[0], c[0], 1);
		fp24_frb(c[1], c[1], 1);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					fp2_mul_frb(c[1][j][k][l], c[1][j][k][l], 2, 4);
				}
				if (fp_prime_get_mod8() == 3) {
					fp4_mul_art(c[1][j][k], c[1][j][k]);
				}
			}
			if ((fp_prime_get_mod8() % 4) == 3) {
				fp8_mul_art(c[1][j], c[1][j]);
			}
		}
	}
}

void fp54_frb(fp54_t c, const fp54_t a, int i) {
	/* Cost of 20 multiplication in Fp^2 per Frobenius. */
	fp54_copy(c, a);
	for (; i % 54 > 0; i--) {
		fp18_frb(c[0], c[0], 1);
		fp18_frb(c[1], c[1], 1);
		fp18_frb(c[2], c[2], 1);
		for (int j = 0; j < 2; j++) {
			for (int l = 0; l < 3; l++) {
				fp3_mul_frb(c[1][j][l], c[1][j][l], 2, 3);
				fp3_mul_frb(c[2][j][l], c[2][j][l], 2, 1);
			}
			/* This is not general enough, so hard code parameters needing the
			tweak. */
#if FP_PRIME == 256
			fp9_mul_art(c[1][j], c[1][j]);
			fp9_mul_art(c[1][j], c[1][j]);
			fp9_mul_art(c[2][j], c[2][j]);
#endif
#if FP_PRIME == 446
			fp9_mul_art(c[1][j], c[1][j]);
			fp9_mul_art(c[2][j], c[2][j]);
			fp9_mul_art(c[2][j], c[2][j]);
#endif
		}
	}
}
