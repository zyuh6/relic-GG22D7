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
 * Implementation of point multiplication on prime elliptic curves over a
 * cubic extensions.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if defined(EP_ENDOM)

static void ep3_psi(ep3_t r, const ep3_t p) {
	ep3_t q;

	ep3_null(q);

	if (ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	RLC_TRY {
		ep3_new(q);

		switch (ep_curve_is_pairf()) {
			case EP_SG18:
				/* -3*u = (2*p^2 - p^5) mod r */
				ep3_frb(q, p, 5);
				ep3_frb(r, p, 2);
				ep3_dbl(r, r);
				ep3_sub(r, r, q);
				break;
			case EP_K18:
				/* For KSS18, we have that u = (p^4 - 3*p) mod r. */
				ep3_dbl(q, p);
				ep3_add(q, q, p);
				ep3_frb(r, p, 3);
				ep3_sub(r, r, q);
				ep3_frb(r, r, 1);
				break;
			case EP_FM18:
				/* For FM18, we have that u = (p^4-p) mod r. */
				ep3_frb(q, p, 3);
				ep3_sub(r, q, p);
				ep3_frb(r, r, 1);
				break;
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep3_free(q);
	}
}

#if EP_MUL == LWNAF || !defined(STRIP)

static void ep3_mul_glv_imp(ep3_t r, const ep3_t p, const bn_t k) {
	int i, j;
	size_t l, _l[6];
	bn_t n, _k[6], u;
	int8_t naf[6][RLC_FP_BITS + 1];
	ep3_t q[6];

	bn_null(n);
	bn_null(u);

	RLC_TRY {
		bn_new(n);
		bn_new(u);
		for (i = 0; i < 6; i++) {
			bn_null(_k[i]);
			ep3_null(q[i]);
			bn_new(_k[i]);
			ep3_new(q[i]);
		}

		fp_prime_get_par(u);
		if (ep_curve_is_pairf() == EP_SG18) {
			/* Compute base -3*u for the recoding below. */
			bn_dbl(n, u);
			bn_add(u, u, n);
			bn_neg(u, u);
		}
		ep3_curve_get_ord(n);
		bn_mod(_k[0], k, n);
		bn_rec_frb(_k, 6, _k[0], u, n, ep_curve_is_pairf() == EP_BN);

		ep3_norm(q[0], p);
		for (int i = 1; i < 6; i++) {
			ep3_psi(q[i], q[i - 1]);
		}
#if defined(EP_MIXED)
		ep3_norm_sim(q + 1, q + 1, 5);
#endif

		l = 0;
		for (i = 0; i < 6; i++) {
			if (bn_sign(_k[i]) == RLC_NEG) {
				ep3_neg(q[i], q[i]);
			}
			_l[i] = RLC_FP_BITS + 1;
			bn_rec_naf(naf[i], &_l[i], _k[i], 2);
			l = RLC_MAX(l, _l[i]);
		}

		ep3_set_infty(r);
		for (j = l - 1; j >= 0; j--) {
			ep3_dbl(r, r);

			for (i = 0; i < 6; i++) {
				if (naf[i][j] > 0) {
					ep3_add(r, r, q[i]);
				}
				if (naf[i][j] < 0) {
					ep3_sub(r, r, q[i]);
				}
			}
		}

		/* Convert r to affine coordinates. */
		ep3_norm(r, r);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(u);
		for (i = 0; i < 3; i++) {
			bn_free(_k[i]);
			ep3_free(q[i]);
		}

	}
}

#endif /* EP_MUL == LWNAF */

#if EP_MUL == LWREG || !defined(STRIP)

static void ep3_mul_reg_gls(ep3_t r, const ep3_t p, const bn_t k) {
	int8_t reg[6][RLC_FP_BITS + 1], b[6], s[6], c0, n0;
	ep3_t q, w, t[6][1 << (RLC_WIDTH - 2)];
	bn_t n, _k[6], u;
	size_t l, len, _l[6];

	bn_null(n);
	bn_null(u);
	ep3_null(q);
	ep3_null(w);

	RLC_TRY {
		bn_new(n);
		bn_new(u);
		ep3_new(q);
		ep3_new(w);
		for (size_t i = 0; i < 6; i++) {
			bn_null(_k[i]);
			bn_new(_k[i]);
			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				ep3_null(t[i][j]);
				ep3_new(t[i][j]);
			}
		}

		ep3_curve_get_ord(n);
		fp_prime_get_par(u);
		bn_mod(_k[0], k, n);
		bn_rec_frb(_k, 6, _k[0], u, n, ep_curve_is_pairf() == EP_BN);

		l = 0;
		/* Make some extra room for BN curves that grow subscalars by 1. */
		len = bn_bits(u) + (ep_curve_is_pairf() == EP_BN);
		ep3_norm(t[0][0], p);
		for (size_t i = 0; i < 6; i++) {
			s[i] = bn_sign(_k[i]);
			bn_abs(_k[i], _k[i]);
			b[i] = bn_is_even(_k[i]);
			_k[i]->dp[0] |= b[i];

			_l[i] = RLC_FP_BITS + 1;
			bn_rec_reg(reg[i], &_l[i], _k[i], len, RLC_WIDTH);
			l = RLC_MAX(l, _l[i]);
			
			/* Apply Frobenius before flipping sign to build table. */
			if (i > 0) {
				ep3_psi(t[i][0], t[i - 1][0]);
			}
		}

		for (size_t i = 0; i < 6; i++) {
			ep3_neg(q, t[i][0]);
			fp3_copy_sec(q->y, t[i][0]->y, s[i] == RLC_POS);
			ep3_tab(t[i], q, RLC_WIDTH);
		}

#if defined(EP_MIXED)
		fp3_set_dig(w->z, 1);
		w->coord = BASIC;
#else
		w->coord = = EP_ADD;
#endif

		ep3_set_infty(r);
		for (int j = l - 1; j >= 0; j--) {
			for (size_t i = 0; i < RLC_WIDTH - 1; i++) {
				ep3_dbl(r, r);
			}

			for (size_t i = 0; i < 6; i++) {
				n0 = reg[i][j];
				c0 = (n0 >> 7);
				n0 = ((n0 ^ c0) - c0) >> 1;

				for (size_t m = 0; m < (1 << (RLC_WIDTH - 2)); m++) {
					fp3_copy_sec(w->x, t[i][m]->x, m == n0);
					fp3_copy_sec(w->y, t[i][m]->y, m == n0);
	#if !defined(EP_MIXED)
					fp3_copy_sec(w->z, t[i][m]->z, m == n0);
	#endif
				}

				ep3_neg(q, w);
				fp3_copy_sec(q->y, w->y, c0 == 0);
				ep3_add(r, r, q);
			}
		}

		for (size_t i = 0; i < 6; i++) {
			/* Tables are built with points already negated, so no need here. */
			ep3_sub(q, r, t[i][0]);
			fp3_copy_sec(r->x, q->x, b[i]);
			fp3_copy_sec(r->y, q->y, b[i]);
			fp3_copy_sec(r->z, q->z, b[i]);
		}

		/* Convert r to affine coordinates. */
		ep3_norm(r, r);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(u);
		ep3_free(q);
		ep3_free(w);
		for (int i = 0; i < 6; i++) {
			bn_free(_k[i]);
			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				ep3_free(t[i][j]);
			}
		}
	}
}

#endif /* EP_MUL == LWREG */
#endif /* EP_ENDOM */

#if defined(EP_PLAIN) || defined(EP_SUPER)

#if EP_MUL == LWNAF || !defined(STRIP)

static void ep3_mul_naf_imp(ep3_t r, const ep3_t p, const bn_t k) {
	int i, n;
	int8_t naf[RLC_FP_BITS + 1];
	ep3_t t[1 << (RLC_WIDTH - 2)];
	size_t l;

	RLC_TRY {
		/* Prepare the precomputation table. */
		for (i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep3_null(t[i]);
			ep3_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep3_tab(t, p, RLC_WIDTH);

		/* Compute the w-NAF representation of k. */
		l = sizeof(naf);
		bn_rec_naf(naf, &l, k, RLC_WIDTH);

		ep3_set_infty(r);
		for (i = l - 1; i >= 0; i--) {
			ep3_dbl(r, r);

			n = naf[i];
			if (n > 0) {
				ep3_add(r, r, t[n / 2]);
			}
			if (n < 0) {
				ep3_sub(r, r, t[-n / 2]);
			}
		}
		/* Convert r to affine coordinates. */
		ep3_norm(r, r);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		/* Free the precomputation table. */
		for (i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep3_free(t[i]);
		}
	}
}

#endif /* EP_MUL == LWNAF */

#if EP_MUL == LWREG || !defined(STRIP)

static void ep3_mul_reg_imp(ep3_t r, const ep3_t p, const bn_t k) {
	bn_t _k;
	int8_t s, reg[1 + RLC_CEIL(RLC_FP_BITS + 1, RLC_WIDTH - 1)];
	ep3_t t[1 << (RLC_WIDTH - 2)], u, v;
	size_t l, n;

	bn_null(_k);

	RLC_TRY {
		bn_new(_k);
		ep3_new(u);
		ep3_new(v);
		/* Prepare the precomputation table. */
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep3_null(t[i]);
			ep3_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep3_tab(t, p, RLC_WIDTH);

		ep3_curve_get_ord(_k);
		n = bn_bits(_k);

		/* Make a copy of the scalar for processing. */
		bn_abs(_k, k);
		_k->dp[0] |= 1;

		/* Compute the regular w-NAF representation of k. */
		l = RLC_CEIL(n, RLC_WIDTH - 1) + 1;
		bn_rec_reg(reg, &l, _k, n, RLC_WIDTH);

#if defined(EP_MIXED)
		fp3_set_dig(u->z, 1);
		u->coord = BASIC;
#else
		u->coord = EP_ADD;
#endif
		ep3_set_infty(r);
		for (int i = l - 1; i >= 0; i--) {
			for (size_t j = 0; j < RLC_WIDTH - 1; j++) {
				ep3_dbl(r, r);
			}

			n = reg[i];
			s = (n >> 7);
			n = ((n ^ s) - s) >> 1;

			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				fp3_copy_sec(u->x, t[j]->x, j == n);
				fp3_copy_sec(u->y, t[j]->y, j == n);
#if !defined(EP_MIXED)
				fp_copy_sec(u->z, t[j]->z, j == n);
#endif
			}
			ep3_neg(v, u);
			fp3_copy_sec(u->y, v->y, s != 0);
			ep3_add(r, r, u);
		}
		/* t[0] has an unmodified copy of p. */
		ep3_sub(u, r, t[0]);
		fp3_copy_sec(r->x, u->x, bn_is_even(k));
		fp3_copy_sec(r->y, u->y, bn_is_even(k));
		fp3_copy_sec(r->z, u->z, bn_is_even(k));
		/* Convert r to affine coordinates. */
		ep3_norm(r, r);
		ep3_neg(u, r);
		fp3_copy_sec(r->y, u->y, bn_sign(k) == RLC_NEG);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		/* Free the precomputation table. */
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep3_free(t[i]);
		}
		bn_free(_k);
		ep3_free(u);
		ep3_free(v);
	}
}

#endif /* EP_MUL == LWREG */
#endif /* EP_PLAIN || EP_SUPER */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep3_mul_basic(ep3_t r, const ep3_t p, const bn_t k) {
	ep3_t t;
	int8_t u, *naf = RLC_ALLOCA(int8_t, bn_bits(k) + 1);
	size_t l;

	ep3_null(t);

	if (bn_is_zero(k) || ep3_is_infty(p)) {
		RLC_FREE(naf);
		ep3_set_infty(r);
		return;
	}

	if (bn_bits(k) <= RLC_DIG) {
		ep3_mul_dig(r, p, k->dp[0]);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
		RLC_FREE(naf);
		return;
	}

	RLC_TRY {
		ep3_new(t);
		if (naf == NULL) {
			RLC_THROW(ERR_NO_BUFFER);
		}

		l = bn_bits(k) + 1;
		bn_rec_naf(naf, &l, k, 2);
		ep3_set_infty(t);
		for (int i = l - 1; i >= 0; i--) {
			ep3_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep3_add(t, t, p);
			} else if (u < 0) {
				ep3_sub(t, t, p);
			}
		}

		ep3_norm(r, t);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep3_free(t);
		RLC_FREE(naf);
	}
}

#if EP_MUL == SLIDE || !defined(STRIP)

void ep3_mul_slide(ep3_t r, const ep3_t p, const bn_t k) {
	ep3_t t[1 << (RLC_WIDTH - 1)], q;
	int i, j;
	size_t l;
	uint8_t win[RLC_FP_BITS + 1];

	ep3_null(q);

	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	RLC_TRY {
		for (i = 0; i < (1 << (RLC_WIDTH - 1)); i ++) {
			ep3_null(t[i]);
			ep3_new(t[i]);
		}

		ep3_new(q);

		ep3_copy(t[0], p);
		ep3_dbl(q, p);

#if defined(EP_MIXED)
		ep3_norm(q, q);
#endif

		/* Create table. */
		for (i = 1; i < (1 << (RLC_WIDTH - 1)); i++) {
			ep3_add(t[i], t[i - 1], q);
		}

#if defined(EP_MIXED)
		ep3_norm_sim(t + 1, t + 1, (1 << (RLC_WIDTH - 1)) - 1);
#endif

		ep3_set_infty(q);
		l = RLC_FP_BITS + 1;
		bn_rec_slw(win, &l, k, RLC_WIDTH);
		for (i = 0; i < l; i++) {
			if (win[i] == 0) {
				ep3_dbl(q, q);
			} else {
				for (j = 0; j < util_bits_dig(win[i]); j++) {
					ep3_dbl(q, q);
				}
				ep3_add(q, q, t[win[i] >> 1]);
			}
		}

		ep3_norm(r, q);
		if (bn_sign(k) == RLC_NEG) {
			ep3_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		for (i = 0; i < (1 << (RLC_WIDTH - 1)); i++) {
			ep3_free(t[i]);
		}
		ep3_free(q);
	}
}

#endif

#if EP_MUL == MONTY || !defined(STRIP)

void ep3_mul_monty(ep3_t r, const ep3_t p, const bn_t k) {
	ep3_t t[2];
	bn_t n, l, _k;
	size_t bits;

	bn_null(n);
	bn_null(l);
	bn_null(_k);
	ep3_null(t[0]);
	ep3_null(t[1]);

	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	RLC_TRY {
		bn_new(n);
		bn_new(l);
		bn_new(_k);
		ep3_new(t[0]);
		ep3_new(t[1]);

		ep3_curve_get_ord(n);
		bits = bn_bits(n);

		bn_mod(_k, k, n);
		bn_abs(l, _k);
		bn_add(l, l, n);
		bn_add(n, l, n);
		dv_swap_sec(l->dp, n->dp, RLC_MAX(l->used, n->used),
			bn_get_bit(l, bits) == 0);
		l->used = RLC_SEL(l->used, n->used, bn_get_bit(l, bits) == 0);

		ep3_norm(t[0], p);
		ep3_dbl(t[1], t[0]);

		/* Blind both points independently. */
		ep3_blind(t[0], t[0]);
		ep3_blind(t[1], t[1]);

		for (int i = bits - 1; i >= 0; i--) {
			int j = bn_get_bit(l, i);
			dv_swap_sec(t[0]->x[0], t[1]->x[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[1], t[1]->x[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[2], t[1]->x[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[0], t[1]->y[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[1], t[1]->y[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[2], t[1]->y[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[0], t[1]->z[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[1], t[1]->z[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[2], t[1]->z[2], RLC_FP_DIGS, j ^ 1);
			ep3_add(t[0], t[0], t[1]);
			ep3_dbl(t[1], t[1]);
			dv_swap_sec(t[0]->x[0], t[1]->x[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[1], t[1]->x[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[2], t[1]->x[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[0], t[1]->y[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[1], t[1]->y[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[2], t[1]->y[2], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[0], t[1]->z[0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[1], t[1]->z[1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[2], t[1]->z[2], RLC_FP_DIGS, j ^ 1);
		}

		ep3_norm(r, t[0]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(l);
		bn_free(_k);
		ep3_free(t[1]);
		ep3_free(t[0]);
	}
}

#endif

#if EP_MUL == LWNAF || !defined(STRIP)

void ep3_mul_lwnaf(ep3_t r, const ep3_t p, const bn_t k) {
	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep3_mul_glv_imp(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep3_mul_naf_imp(r, p, k);
#endif
}

#endif

#if EP_MUL == LWREG || !defined(STRIP)

void ep3_mul_lwreg(ep3_t r, const ep3_t p, const bn_t k) {
	if (bn_is_zero(k) || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep3_mul_reg_gls(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep3_mul_reg_imp(r, p, k);
#endif
}

#endif

void ep3_mul_gen(ep3_t r, const bn_t k) {
	if (bn_is_zero(k)) {
		ep3_set_infty(r);
		return;
	}

#ifdef EP_PRECO
	ep3_mul_fix(r, ep3_curve_get_tab(), k);
#else
	ep3_t g;

	ep3_null(g);

	RLC_TRY {
		ep3_new(g);
		ep3_curve_get_gen(g);
		ep3_mul(r, g, k);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep3_free(g);
	}
#endif
}

void ep3_mul_dig(ep3_t r, const ep3_t p, const dig_t k) {
	ep3_t t;
	bn_t _k;
	int8_t u, naf[RLC_DIG + 1];
	size_t l;

	ep3_null(t);
	bn_null(_k);

	if (k == 0 || ep3_is_infty(p)) {
		ep3_set_infty(r);
		return;
	}

	RLC_TRY {
		ep3_new(t);
		bn_new(_k);

		bn_set_dig(_k, k);

		l = RLC_DIG + 1;
		bn_rec_naf(naf, &l, _k, 2);

		ep3_copy(t, p);
		for (int i = l - 2; i >= 0; i--) {
			ep3_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep3_add(t, t, p);
			} else if (u < 0) {
				ep3_sub(t, t, p);
			}
		}

		ep3_norm(r, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep3_free(t);
		bn_free(_k);
	}
}
