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
 * Implementation of the pairings over prime curves.
 *
 * @ingroup pp
 */

#include "relic_core.h"
#include "relic_pp.h"
#include "relic_util.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/* Precompute the constants in the dual map of 2-isogeny */
#define phihat_alpha1		"EF3D2552498A7063D6BFFCCE4B1F0DCCBD8527D0AF776CC823E38B4BEE4665D039EB548F1C2DC001D02B1769DB4AE1CDB86859E04605971BD4"
#define phihat_alpha323		"25AC06450C1B586E69988FD90FCC583BB54E6A57437DD3EC75FD603234CE353D18E38E2D882324164895B279968B161A1482E19F36F864EB1E"
#define phihat_U1			"11A81ADE0C01202FE2B2BA882AF7AD25E150823944631F5196727BEF46B322DA70926B66904CC77B27C781AE6467489CB8809059B853FFCFAAA"
#define phihat_U13			"BC5551A1F3F3105CEC23F6F5E2FC500F8EB81885C176F027DC88D22AEDC055AEEE6E9D59AEE6D5454E30C971D5055822585122DD33A13243D3"
/**
 * Compute the values of phihat(P) and [2]P
 * where phihat is the dual of 2-isogeny
 * @param[out] R1			- the resulting point phihat(P).
 * @param[out] R2			- the resulting point [2]P.
 * @param[in] P1			- P \in G1
 */
static void ep11_phihat_2(ep_t R1, ep_t R2, const ep_t P) {
	char str[2 * RLC_FP_BYTES + 1];
	ep_t RR1, RR2;
	fp_t alpha1, alpha323, U1, U13;
	fp_t t0, t1, t2, t3, t4, t5, t6;
	ep_null(RR1);
	ep_null(RR2);
	fp_null(alpha1);
	fp_null(alpha323);
	fp_null(U1);
	fp_null(U13);
	fp_null(t0);
	fp_null(t1);
	fp_null(t2);
	fp_null(t3);
	fp_null(t4);
	fp_null(t5);
	fp_null(t6);
	RLC_TRY {
		ep_new(RR1);
		ep_new(RR2);
		fp_new(alpha1);
		fp_new(alpha323);
		fp_new(U1);
		fp_new(U13);
		fp_new(t0);
		fp_new(t1);
		fp_new(t2);
		fp_new(t3);
		fp_new(t4);
		fp_new(t5);
		fp_new(t6);

		RLC_GET( str, phihat_alpha1, sizeof(phihat_alpha1));
		fp_read_str(alpha1, str, strlen(str), 16);
		RLC_GET( str, phihat_alpha323, sizeof(phihat_alpha323));
		fp_read_str(alpha323, str, strlen(str), 16);
		RLC_GET( str, phihat_U1, sizeof(phihat_U1));
		fp_read_str(U1, str, strlen(str), 16);
		RLC_GET( str, phihat_U13, sizeof(phihat_U13));
		fp_read_str(U13, str, strlen(str), 16);

		fp_sub(t0,P->x,alpha1);
		fp_sqr(t1,t0);
		fp_mul(RR1->x,t0,P->x);
		fp_add(RR1->x,RR1->x,alpha323);
		fp_sub(RR1->y,t1,alpha323);
		fp_mul(RR1->y,RR1->y,P->y);

		fp_sqr(t3,P->x);
		fp_sqr(t4,P->y);
		fp_sqr(t2,t4);
		fp_add(t5,P->x,t4);
		fp_sqr(t5,t5);
		fp_sub(t5,t5,t3);
		fp_sub(t5,t5,t2);
		fp_add(t5,t5,t5);

		fp_sub_dig(t4,t3,1);
		fp_add(t3,t4,t4);
		fp_add(t4,t4,t3);

		fp_sqr(RR2->x,t4);
		fp_sub(RR2->x,RR2->x,t5);
		fp_sub(RR2->x,RR2->x,t5);

		fp_sub(RR2->y,t5,RR2->x);
		fp_mul(RR2->y,RR2->y,t4);
		fp_add(t2,t2,t2);
		fp_add(t2,t2,t2);
		fp_add(t2,t2,t2);
		fp_sub(RR2->y,RR2->y,t2);

		//batch inverse
		fp_mul(t3,U13,t1);
		fp_add(t5,P->y,P->y);
		fp_mul(t4,t5,t3);
		fp_inv(t4,t4);
		fp_mul(t2,t4,t5);
		fp_mul(RR1->y,RR1->y,t2);
		fp_mul(t2,t2,U1);
		fp_mul(t2,t2,t0);
		fp_mul(RR1->x,RR1->x,t2);

		fp_mul(t2,t4,t3);
		fp_sqr(t4,t2);
		fp_mul(RR2->x,RR2->x,t4);
		fp_mul(t4,t4,t2);
		fp_mul(RR2->y,RR2->y,t4);

		fp_set_dig(RR1->z,1);
		fp_set_dig(RR2->z,1);
		RR1->coord=BASIC;
		RR2->coord=BASIC;

		ep_copy(R1,RR1);
		ep_copy(R2,RR2);

	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep_free(RR1);
		ep_free(RR2);
		fp_free(alpha1);
		fp_free(alpha323);
		fp_free(U1);
		fp_free(U13);
		fp_free(t0);
		fp_free(t1);
		fp_free(t2);
		fp_free(t3);
		fp_free(t4);
		fp_free(t5);
		fp_free(t6);
	}
}

/**
 * Compute the point Q + T2
 * where T2 is a 2-order point on the curve
 * @param[out] RR			- the resulting point Q+T2 \in G2.
 * @param[in] Q				- Q \in the twist of G2
 */
static void ep11_QpT2(ep11_t RR, const ep11_t Q) {
	fp_t T2;
	fp11_t t0,t1,t2,t3,t4;
	ep11_t R;
	fp_null(T2);
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	fp11_null(t3);
	fp11_null(t4);
	ep11_null(R);

	RLC_TRY{
		fp_new(T2);
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		fp11_new(t3);
		fp11_new(t4);
		ep11_new(R);

		ep11_curve_get_T2x(T2);
		fp11_mul_art(t0,Q->x);
		fp11_mul_art(t1,Q->y);
		fp11_neg(t2,t0);
		fp_add(t2[0],T2,t2[0]);
		fp11_add(R->z,t2,t2);
		fp11_sqr(t3,R->z);
		fp11_mul(t2,t2,t3);
		fp11_mul(t4,t3,t0);
		fp11_add(t3,t1,t1);
		fp11_sqr(R->x,t3);
		fp11_mul_art(R->x,R->x);
		fp11_add(t0,t2,t4);
		fp11_sub(R->x,R->x,t0);
		fp11_sub(R->x,R->x,t4);
		fp11_sub(t0,R->x,t0);
		fp11_mul(R->y,t0,t3);

		R->coord = PROJC;
		ep11_copy(RR,R);	
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp_free(T2);
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		fp11_free(t3);
		fp11_free(t4);
		ep11_free(R);
	}
}

/**
 * Compute the point [-x]Q + T2
 * where T2 is a 2-order point on the curve
 * @param[out] RR			- the resulting point [-x]Q+T2.
 * @param[in] xQ			- [-x]Q \in the twist of G2
 */
static void ep11_xQpT2(ep11_t RR, const ep11_t xQ){
	fp11_t t0, t1, ZZ, U, H, I, J, V, K;
	ep11_t R;
	fp_t T2;
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(ZZ);
    fp11_null(U);
    fp11_null(H);
    fp11_null(I);
    fp11_null(J);
    fp11_null(V);
	fp11_null(K);
	ep11_null(R);
	fp_null(T2);
	RLC_TRY {
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(ZZ);
		fp11_new(U);
		fp11_new(H);
		fp11_new(I);
		fp11_new(J);
		fp11_new(V);
		fp11_new(K);
		ep11_new(R);
		fp_new(T2);

		fp11_mul_art(t0,xQ->x);
		fp11_mul_art(t1,xQ->y);
		ep11_curve_get_T2x(T2);
		fp11_sqr(ZZ, xQ->z);
		for(int i=0;i<11;i++) fp_mul(U[i],T2,ZZ[i]);
		fp11_sub(H,U,t0);
		fp11_add(K,H,H);
		fp11_mul(R->z, K, xQ->z);

		fp11_sqr(I,K);
		fp11_mul(J,H,I);
		fp11_mul(V,t0,I);
		fp11_add(J,J,V);
		fp11_add(K,t1,t1);
		fp11_sqr(R->x,K);
		fp11_mul_art(R->x,R->x);
		fp11_sub(R->x,R->x,J);
		fp11_sub(R->x,R->x,V);

		fp11_sub(R->y,R->x,J);
		fp11_mul(R->y,R->y,K);
		
		R->coord = PROJC;
		ep11_copy(RR,R);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(ZZ);
		fp11_free(U);
		fp11_free(H);
		fp11_free(I);
		fp11_free(J);
		fp11_free(V);
		fp11_free(K);
		ep11_free(R);
		fp_free(T2);
	}
}

/**
 * Compute the  Miller loop  for optimal ate pairings of type G_2 x G_1 with k=22.
 * where seed = -779523
 * @param[out] value		- the value of f^{p-x}_{-x,Q}(P) f_{-x,-xQ}(P)
 * @param[out] R1			- the resulting point [-x]Q.
 * @param[out] R2			- the resulting point [2]Q.
 * @param[out] l1			- the value of l^{p^2}_{Q,Q}(P)
 * @param[in] Q				- Q \in the twist of G2
 * @param[in] P				- P \in G1
 * @param[in] x				- x = 779523 is positive
 */
static void pp_mil_k22(fp22_t value, ep11_t R1, ep11_t R2, fp22_t l1, const ep11_t Q, const ep_t P, bn_t x) {
	ep11_t T, _Q, _xQ;
	fp22_t f, g, _g, L;
	size_t len = bn_bits(x) + 1;
	int8_t s[RLC_FP_BITS + 1];
	ep11_null(T);
	ep11_null(_Q);
	ep11_null(_xQ);
	fp22_null(f);
	fp22_null(g);
	fp22_null(_g);
	fp22_null(L);

	RLC_TRY{
		ep11_new(T);
		ep11_new(_Q);
		ep11_new(_xQ);
		fp22_new(f);
		fp22_new(g);
		fp22_new(_g);
		fp22_new(L);

		ep11_copy(T,Q);
		ep11_neg(_Q,Q);
		bn_rec_naf(s, &len, x, 2);

		pp_dbl1_k22(L, T, T, P);
		fp22_copy(f,L);
		fp22_frb(l1,L,2);
		ep11_copy(R2,T);

		if (s[len-2] > 0) {
			pp_add_k22(L, T, T, Q, P);
			fp22_mul(f,f,L);
		} else if (s[len-2] < 0) {
			pp_add_k22(L, T, T, _Q, P);
			fp22_mul(f,f,L);
		}
		for(int i = len-3;i>=0;i--) {
			pp_dbl_k22(L, T, T, P);
			fp22_sqr(f,f);
			fp22_mul(f,f,L);

			if(s[i] > 0){
				pp_add_k22(L, T, T, Q, P);
				fp22_mul(f,f,L);
			} else if (s[i] < 0) {
				pp_add_k22(L, T, T, _Q, P);
				fp22_mul(f,f,L);
			}
		}

		fp22_copy(g,f);
		fp22_inv_cyc(_g,g);
		ep11_norm(T,T);//it will be faster
		ep11_copy(R1,T);
		ep11_neg(_xQ,R1);

		pp_dbl1_k22(L, T, T, P);
		fp22_sqr(f,f);
		fp22_mul(f,f,L);
		if (s[len-2] > 0) {
			pp_add_k22(L, T, T, R1, P);
			fp22_mul(f,f,L);
		} else if (s[len-2] < 0) {
			pp_add_k22(L, T, T, _xQ, P);
			fp22_mul(f,f,L);
		}

		for(int i = len-3; i>=0; i--) {
			pp_dbl_k22(L, T, T, P);
			fp22_sqr(f,f);
			fp22_mul(f,f,L);
			if(s[i] > 0){
				pp_add_k22(L, T, T, R1, P);
				fp22_mul(f,f,L);
				fp22_mul(f,f,g);
			} else if (s[i] < 0) {
				pp_add_k22(L, T, T, _xQ, P);
				fp22_mul(f,f,L);
				fp22_mul(f,f,_g);
			}
		}
		
		fp22_frb(g,g,1);
		fp22_mul(f,f,g);
		fp22_copy(value,f);

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		ep11_free(T);
		ep11_free(_Q);
		ep11_free(_xQ);
		fp22_free(f);
		fp22_free(g);
		fp22_free(_g);
		fp22_free(L);
	}
}

/**
 * Compute the value of line function l_{Q+T2,Q+T2}(phihat(P))
 * where phihat is the dual of 2-isogeny and T2 is a 2-order point on the curve 
 * @param[out] M1			- the value of l_{Q+T2,Q+T2}(phihat(P))
 * @param[in] Q				- Q+T2 \in G2
 * @param[in] P				- P \in G1
 */
static void pp_M1_k22(fp22_t M1, const ep11_t Q, const ep_t P){
	fp11_t t0,t2,t4;
	dv11_t u0, u1;
	fp11_null(t0);
	fp11_null(t2);
	fp11_null(t4);
	dv11_null(u0);
	dv11_null(u1);

	RLC_TRY{
		fp11_new(t0);
		fp11_new(t2);
		fp11_new(t4);
		dv11_new(u0);
		dv11_new(u1);

		fp11_sqr(t2,Q->z);
		fp11_mul(t0,t2,Q->z);
		fp11_mul(M1[0],Q->y,t0);
		fp11_mul_art(M1[0],M1[0]);
		for(int i=0;i<11;i++) fp_mul(M1[0][i],M1[0][i],P->y);
		fp11_add(M1[0],M1[0],M1[0]);

		fp11_sqrn_low(u0,Q->y);
		fp11_mul_nor_low(u0,u0);
		fp11_addd_low(u0,u0,u0);
		fp11_add(t0,t2,Q->x);
		fp11_sub(t4,Q->x,t2);
		fp11_mul(t0,t0,t4);
		fp11_add(t4,t0,t0);
		fp11_add(t4,t4,t0);
		for(int i=0;i<11;i++) fp_mul(t0[i],P->x,t2[i]);
		fp11_sub(t0,t0,Q->x);
		fp11_muln_low(u1,t0,t4);
		fp11_addd_low(u0,u0,u1);
		fp11_rdc(t0,u0);
		fp11_neg(M1[1],t0);
		

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t2);
		fp11_free(t4);
		dv11_free(u0);
		dv11_free(u1);
	}
}

/**
 * Compute the value of line function l_{[-x]Q,Q+T2}(phihat(P))
 * where phihat is the dual of 2-isogeny and T2 is a 2-order point on the curve 
 * @param[out] M34			- the value of l_{[-x]Q,Q+T2}(phihat(P))
 * @param[in] QQ1			- [-x]Q \in the twist of G2
 * @param[in] Q2			- Q+T2 \in G2
 * @param[in] P				- P \in G1
 */
static void pp_M34_k22(fp22_t M34, const ep11_t QQ1, const ep11_t Q2, const ep_t P){
	ep11_t Q1;
	fp11_t t0,t1,t2,t3;
	dv11_t u0,u1,u2;
	ep11_null(Q1);
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	fp11_null(t3);
	dv11_null(u0);
	dv11_null(u1);
	dv11_null(u2);


	RLC_TRY{
		ep11_new(Q1);
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		fp11_new(t3);
		dv11_new(u0);
		dv11_new(u1);
		dv11_new(u2);

		fp11_copy(Q1->z,QQ1->z);
		fp11_mul_art(Q1->x,QQ1->x);
		fp11_mul_art(Q1->y,QQ1->y);

		fp11_sqr(t0,Q1->z);
		fp11_mul(t1,t0,Q1->z);
		fp11_sqr(t2,Q2->z);
		fp11_mul(t2,t2,Q2->z);

		fp11_mul(t3,Q2->x,Q2->z);
		fp11_muln_low(u0,t0,t3);
		fp11_muln_low(u1,Q1->x,t2);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(t3,u0);
		fp11_mul(M34[0],t3,t1);
		for(int i=0;i<11;i++)fp_mul(M34[0][i],M34[0][i],P->y);

		fp11_muln_low(u0,t1,Q2->y);
		fp11_muln_low(u1,t2,Q1->y);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(t1,u0);
		for(int i=0;i<11;i++)fp_mul(t0[i],t0[i],P->x);
		fp11_sub(t0,t0,Q1->x);
		fp11_muln_low(u1,t0,t1);
		fp11_muln_low(u0,Q1->y,t3);
		fp11_addd_low(u0,u0,u1);
		fp11_rdc(t0,u0);
		fp11_neg(M34[1],t0);

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		ep11_free(Q1);
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		fp11_free(t3);
		dv11_free(u0);
		dv11_free(u1);
		dv11_free(u2);
	}
}

/**
 * Compute the value of line function l_{[-x]Q+T2,Q}(phihat(P))
 * where phihat is the dual of 2-isogeny and T2 is a 2-order point on the curve 
 * @param[out] M42			- the value of l_{[-x]Q+T2,Q}(phihat(P))
 * @param[in] QQ1			- Q \in the twist of G2
 * @param[in] Q2			- [-x]Q+T2 \in G2
 * @param[in] P				- P \in G1
 */
static void pp_M42_k22(fp22_t M42, const ep11_t QQ1, const ep11_t Q2, const ep_t P){
	fp11_t t0,t1,t2;
	dv11_t u0,u1;
	ep11_t Q1;
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	dv11_null(u0);
	dv11_null(u1);
	ep11_null(Q1);

	RLC_TRY{
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		dv11_new(u0);
		dv11_new(u1);
		ep11_new(Q1);

		fp11_mul_art(Q1->x,QQ1->x);
		fp11_mul_art(Q1->y,QQ1->y);

		fp11_sqr(t1,Q2->z);
		fp11_mul(t1,t1,Q2->z);
		fp11_muln_low(u0,Q2->x,Q2->z);
		fp11_muln_low(u1,Q1->x,t1);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(t0,u0);
		for(int i=0;i<11;i++)fp_mul(M42[0][i],t0[i],P->y);

		fp11_mul(t1,Q1->y,t1);
		fp11_sub(t1,Q2->y,t1);
		fp11_neg(t2,Q1->x);
		fp_add(t2[0],t2[0],P->x);
		fp11_muln_low(u1,t1,t2);
		fp11_muln_low(u0,t0,Q1->y);
		fp11_addd_low(u0,u1,u0);
		fp11_rdc(t0,u0);
		fp11_neg(M42[1],t0);

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		dv11_free(u0);
		dv11_free(u1);
	}
}


/**
 * Compute the  Miller loop  for super-optimal ate pairings of type G_2 x G_1 with k=22.
 * where seed = -779523
 * @param[out] value		- the value of (f^{2p^21}_{-x,Q}([2]P)*M1/M2)^{(-x+1)/4} f_{-x,Q}(phihat(P))
 * @param[out] R1			- the resulting point [-x]Q.
 * @param[out] R2			- the resulting point [2]Q.
 * @param[out] M2			- the value of l_{Q,Q}(phihat(P)) for M2
 * @param[out] M32			- the value of l^{p}_{Q,Q}([2]P) for M3
 * @param[out] M31			- the value of f_{-x,Q}(2P) for M3
 * @param[out] M41			- the value of f^{p^21}_{-x,Q}(2P) for M4
 * @param[in] Q				- Q \in G2
 * @param[in] P1			- phihat(P) \in G1
 * @param[in] P2			- [2]P \in G1
 * @param[in] x				- x = 779523 is positive
 * @param[in] M1			- M1 = l_{Q+T2,Q+T2}(phihat(P))
 */
static void pp_mil_k22_2iso(fp22_t value, ep11_t R1, ep11_t R2, fp22_t M2, fp22_t M32, fp22_t M31, fp22_t M41, const ep11_t Q, const ep_t P1, const ep_t P2, bn_t x, fp22_t M1){
	ep11_t T, _Q, _xQ;
	fp22_t f, g, h1, _h1, L1, L2, M2c;
	size_t len = bn_bits(x) + 1;
	int8_t s[RLC_FP_BITS + 1];
	fp22_t tab[2*len];
	ep11_null(T);
	ep11_null(_Q);
	ep11_null(_xQ);
	fp22_null(f);
	fp22_null(g);
	fp22_null(h1);
	fp22_null(_h1);
	fp22_null(L1);
	fp22_null(L2);
	fp22_null(M2c);
	for(int i=0;i<2*len;i++) fp22_null(tab[i]);

	RLC_TRY{
		ep11_new(T);
		ep11_new(_Q);
		ep11_new(_xQ);
		fp22_new(f);
		fp22_new(g);
		fp22_new(h1);
		fp22_new(_h1);
		fp22_new(L1);
		fp22_new(L2);
		fp22_new(M2c);
		for(int i=0;i<2*len;i++) fp22_new(tab[i]);

		ep11_copy(T,Q);
		ep11_neg(_Q,Q);
		bn_rec_naf(s, &len, x, 2);

		pp_dbl1_k22_2iso(M2, L2, T, T, P1, P2);
		fp22_copy(f,L2);
		fp22_copy(g,M2);
		fp22_frb(M32,L2,1);
		ep11_copy(R2,T);

		if (s[len-2] > 0) {
			pp_add_k22_2iso(L1, L2, T, T, Q, P1, P2);
			fp22_mul(f,f,L2);
			fp22_mul(g,g,L1);
		} else if (s[len-2] < 0) {
			pp_add_k22_2iso(L1, L2, T, T, _Q, P1, P2);
			fp22_mul(f,f,L2);
			fp22_mul(g,g,L1);
		}

		pp_dbl_k22_2iso(L1, L2, T, T, P1, P2);
		fp22_sqr(f,f);
		fp22_mul(f,f,L2);
		fp22_sqr(g,g);
		fp22_mul(g,g,L1);

		if (s[len-3] > 0) {
			pp_add_k22_2iso(L1, L2, T, T, Q, P1, P2);
			fp22_mul(f,f,L2);
			fp22_mul(g,g,L1);
		} else if (s[len-3] < 0) {
			pp_add_k22_2iso(L1, L2, T, T, _Q, P1, P2);
			fp22_mul(f,f,L2);
			fp22_mul(g,g,L1);
		}

		int j=0;
		for(int i = len-4;i>=0;i--) {
			pp_dbl_k22_2iso(tab[j], L2, T, T, P1, P2);
			j++;
			fp22_sqr(f,f);
			fp22_mul(f,f,L2);
			if(s[i] > 0){
				pp_add_k22_2iso(tab[j], L2, T, T, Q, P1, P2);
				j++;
				fp22_mul(f,f,L2);
			} else if (s[i] < 0) {
				pp_add_k22_2iso(tab[j], L2, T, T, _Q, P1, P2);
				j++;
				fp22_mul(f,f,L2);
			}
		}
		ep11_copy(R1,T);
		fp22_copy(M31,f);
		fp22_frb(f,f,21);

		fp22_copy(M41,f);
		fp22_sqr(f,f);
		fp22_mul(f,f,M1);
		fp22_inv_cyc(M2c,M2);
		fp22_mul(f,f,M2c);
		fp22_copy(h1,f);
		fp22_inv_cyc(_h1,h1);
		fp22_mul(f,f,g);

		j=0;
		for(int i=len-4;i>=0;i--){
			fp22_sqr(f,f);
			fp22_mul(f,f,tab[j]);
			j++;
			if(s[i]){
				fp22_mul(f,f,tab[j]);
				j++;
			}
			if(s[i+2]>0) fp22_mul(f,f,h1);
			else if(s[i+2]<0) fp22_mul(f,f,_h1);
		}

		fp22_copy(value,f);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		ep11_free(T);
		ep11_free(_Q);
		ep11_free(_xQ);
		fp22_free(f);
		fp22_free(g);
		fp22_free(h1);
		fp22_free(_h1);
		fp22_free(L1);
		fp22_free(L2);
		fp22_free(M2c);
		for(int i=0;i<2*len;i++) fp22_free(tab[i]);
	}
}

/**
 * Compute the value of line function l_{pi[2]Q,[-x]Q}(P])
 * @param[out] l2			- the value of l_{pi[2]Q,[-x]Q}(P])
 * @param[in] QQ1			- [-x]Q \in the twist of G2, normed
 * @param[in] Q2			- pi[2]Q \in G2
 * @param[in] P				- P \in G1, normed
 */
static void pp_l2_k22(fp22_t l2, const ep11_t QQ1, const ep11_t Q2, const ep_t P){
	ep11_t piQ2, Q1;
	fp11_t t0,t1,t2,t3,t4;
	dv11_t u0,u1;
	ep11_null(piQ2);
	ep11_null(Q1);
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	fp11_null(t3);
	fp11_null(t4);
	dv11_null(u0);
	dv11_null(u1);
	RLC_TRY{
		ep11_new(piQ2);
		ep11_new(Q1);
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		fp11_new(t3);
		fp11_new(t4);
		dv11_new(u0);
		dv11_new(u1);

		ep11_frb(piQ2,Q2,1);
		fp11_mul_art(Q1->x,QQ1->x);
		fp11_mul_art(Q1->y,QQ1->y);

		fp11_sqr(t2,piQ2->z);
		fp11_mul(t3,t2,piQ2->z);

		fp11_muln_low(u0,piQ2->x,piQ2->z);
		fp11_muln_low(u1,Q1->x,t3);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(t4,u0);

		for(int i=0;i<11;i++) fp_mul(l2[0][i],P->y,t4[i]);
		
		fp11_mul(t1,Q1->y,t3);
		fp11_sub(t1,piQ2->y,t1);
		fp11_copy(t2,Q1->x);
		fp_sub(t2[0],t2[0],P->x);
		fp11_muln_low(u0,t1,t2);
		fp11_muln_low(u1,t4,Q1->y);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(l2[1],u0);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		ep11_free(piQ2);
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		fp11_free(t3);
		fp11_free(t4);
		dv11_free(u0);
		dv11_free(u1);
	}
}


/**
 * Compute the value of line function l_{pi[2]Q,[-x]Q}([2]P])
 * @param[out] l2			- the value of l_{pi[2]Q,[-x]Q}([2]P])
 * @param[in] QQ1			- [-x]Q \in the twist of G2
 * @param[in] Q2			- pi[2]Q \in G2
 * @param[in] P				- P \in G1, normed
 */
static void pp_l2_k22_2iso(fp22_t l2, const ep11_t QQ1, const ep11_t Q2, const ep_t P){
	ep11_t piQ2, Q1;
	fp11_t t0,t1,t2,t3,t4;
	dv11_t u0,u1;
	ep11_null(piQ2);
	ep11_null(Q1);
	fp11_null(t0);
	fp11_null(t1);
	fp11_null(t2);
	fp11_null(t3);
	fp11_null(t4);
	dv11_null(u0);
	dv11_null(u1);
	RLC_TRY{
		ep11_new(piQ2);
		ep11_new(Q1);
		fp11_new(t0);
		fp11_new(t1);
		fp11_new(t2);
		fp11_new(t3);
		fp11_new(t4);
		dv11_new(u0);
		dv11_new(u1);

		ep11_frb(piQ2,Q2,1);
		fp11_copy(Q1->z,QQ1->z);
		fp11_mul_art(Q1->x,QQ1->x);
		fp11_mul_art(Q1->y,QQ1->y);

		fp11_sqr(t0,Q1->z);
		fp11_mul(t1,t0,Q1->z);
		fp11_sqr(t2,piQ2->z);
		fp11_mul(t3,t2,piQ2->z);

		fp11_mul(t4,piQ2->x,piQ2->z);
		fp11_muln_low(u0,t4,t0);
		fp11_muln_low(u1,Q1->x,t3);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(t4,u0);

		for(int i=0;i<11;i++) fp_mul(l2[0][i],P->y,t1[i]);
		fp11_mul(l2[0],l2[0],t4);
		
		fp11_muln_low(u0,piQ2->y,t1);
		fp11_muln_low(u1,Q1->y,t3);
		fp11_subc_low(u0,u0,u1);
		fp11_rdc(t1,u0);
		for(int i=0;i<11;i++) fp_mul(t2[i],P->x,t0[i]);
		fp11_sub(t2,t2,Q1->x);
		fp11_muln_low(u0,t1,t2);
		fp11_muln_low(u1,t4,Q1->y);
		fp11_addd_low(u0,u0,u1);
		fp11_rdc(t1,u0);
		fp11_neg(l2[1],t1);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		ep11_free(piQ2);
		fp11_free(t0);
		fp11_free(t1);
		fp11_free(t2);
		fp11_free(t3);
		fp11_free(t4);
		dv11_free(u0);
		dv11_free(u1);
	}
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void pp_map_oatep_k22(fp22_t r, const ep11_t Q, const ep_t P){
	bn_t x;
	ep11_t R1, R2;
	fp22_t l1, l2;
	bn_null(x);
	ep11_null(R1);
	ep11_null(R2);
	fp22_null(l1);
	fp22_null(l2);

	RLC_TRY{
		if(ep_is_infty(P) || ep11_is_infty(Q)) {
			fp22_set_dig(r,1);
			return;
		}
		bn_new(x);
		ep11_new(R1);
		ep11_new(R2);
		fp22_new(l1);
		fp22_new(l2);

		fp_prime_get_par(x);
		bn_neg(x,x);
		pp_mil_k22(r, R1, R2, l1, Q, P, x);
		pp_l2_k22(l2, R1, R2, P);
		fp22_frb(l2,l2,1);

		fp22_mul(r,r,l1);
		fp22_mul(r,r,l2);
		pp_exp_k22(r, r);

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		bn_free(x);
		ep11_free(R1);
		ep11_free(R2);
		fp22_free(l1);
		fp22_free(l2);
	}
	
}

void pp_map_soatep_k22_2iso(fp22_t r, const ep11_t Q, const ep_t P){
	bn_t x;
	ep11_t R1, R2, QpT2, xQpT2;
	fp22_t M1, M2, M31, M32, M33, M34, M41, M42;
	ep_t P1, P2;
	bn_null(x);
	ep11_null(R1);
	ep11_null(R2);
	ep11_null(QpT2);
	ep11_null(xQpT2);
	fp22_null(M1);
	fp22_null(M2);
	fp22_null(M31);
	fp22_null(M32);
	fp22_null(M33);
	fp22_null(M34);
	fp22_null(M41);
	fp22_null(M42);
	ep_null(P1);
	ep_null(P2);

	RLC_TRY{
		if(ep_is_infty(P) || ep11_is_infty(Q)) {
			fp22_set_dig(r,1);
			return;
		}
		bn_new(x);
		ep11_new(R1);
		ep11_new(R2);
		ep11_new(QpT2);
		ep11_new(xQpT2);
		fp22_new(M1);
		fp22_new(M2);
		fp22_new(M31);
		fp22_new(M32);
		fp22_new(M33);
		fp22_new(M34);
		fp22_new(M41);
		fp22_new(M42);
		ep_new(P1);
		ep_new(P2);

		fp_prime_get_par(x);
		bn_neg(x,x);
		ep11_phihat_2(P1,P2,P);
		ep11_QpT2(QpT2,Q);
		pp_M1_k22(M1,QpT2,P1);
		pp_mil_k22_2iso(r, R1, R2, M2, M32, M31, M41, Q, P1, P2, x, M1);
		fp22_sqr(r,r);
		ep11_xQpT2(xQpT2,R1);
		pp_l2_k22_2iso(M33, R1, R2, P2);
		pp_M34_k22(M34, R1, QpT2, P1);
		pp_M42_k22(M42, Q, xQpT2, P1);

		fp22_mul(r,r,M2);
		fp22_mul(r,r,M31);
		fp22_mul(r,r,M32);
		fp22_mul(r,r,M33);
		fp22_mul(r,r,M34);
		fp22_mul(M1,M1,M41);
		fp22_mul(M1,M1,M42);
		fp22_inv_cyc(M1,M1);
		fp22_mul(r,r,M1);
		pp_exp_k22(r, r);

	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
		bn_free(x);
		ep11_free(R1);
		ep11_free(R2);
		ep11_free(QpT2);
		ep11_free(xQpT2);
		fp22_free(M1);
		fp22_free(M2);
		fp22_free(M31);
		fp22_free(M32);
		fp22_free(M33);
		fp22_free(M34);
		fp22_free(M41);
		fp22_free(M42);
		ep_free(P1);
		ep_free(P2);
	}
	
}
