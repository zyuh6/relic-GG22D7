/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2021 RELIC Authors
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
 * quartic extension field.
 *
 * @ingroup epx
 */
#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if defined(EP_ENDOM)

static void ep4_psi(ep4_t r, const ep4_t p) {
	ep4_t q;

	ep4_null(q);

	if (ep4_is_infty(p)) {
		ep4_set_infty(r);
		return;
	}

	RLC_TRY {
		ep4_new(q);

		ep4_copy(r, p);

		switch (ep_curve_is_pairf()) {
			case EP_K16:
				/* u = (2*p^5 - p) mod r */
				ep4_frb(q, p, 1);
				ep4_frb(r, q, 4);
				ep4_dbl(r, r);
				ep4_sub(r, r, q);
				break;
			case EP_N16:
				/* u = -p^5 mod r */
				ep4_frb(r, p, 5);
				ep4_neg(r, r);
				break;
			case EP_FM16:
				/* u = p mod r */
			default:
				ep4_frb(r, p, 1);
				break;
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep4_free(q);
	}
}

#if EP_MUL == LWNAF || !defined(STRIP)

static void ep4_mul_glv_imp(ep4_t r, const ep4_t p, const bn_t k) {
	size_t l, _l[8];
	bn_t n, _k[8], u;
	int8_t naf[8][RLC_FP_BITS + 1];
	ep4_t q[8];

	bn_null(n);
	bn_null(u);

	RLC_TRY {
		bn_new(n);
		bn_new(u);
		for (int i = 0; i < 8; i++) {
			bn_null(_k[i]);
			ep4_null(q[i]);
			bn_new(_k[i]);
			ep4_new(q[i]);
		}

		ep4_curve_get_ord(n);
		fp_prime_get_par(u);
		bn_mod(_k[0], k, n);
		bn_rec_frb(_k, 8, _k[0], u, n, ep_curve_is_pairf() == EP_BN);

		ep4_norm(q[0], p);
		for (size_t i = 1; i < 8; i++) {
            ep4_psi(q[i], q[i - 1]);
		}
#if defined(EP_MIXED)
		ep4_norm_sim(q + 1, q + 1, 7);
#endif

		l = 0;
		for (size_t i = 0; i < 8; i++) {
			if (bn_sign(_k[i]) == RLC_NEG) {
				ep4_neg(q[i], q[i]);
			}
			_l[i] = RLC_FP_BITS + 1;
			bn_rec_naf(naf[i], &_l[i], _k[i], 2);
			l = RLC_MAX(l, _l[i]);
		}

		ep4_set_infty(r);
		for (int j = l - 1; j >= 0; j--) {
			ep4_dbl(r, r);

			for (int i = 0; i < 8; i++) {
				if (naf[i][j] > 0) {
					ep4_add(r, r, q[i]);
				}
				if (naf[i][j] < 0) {
					ep4_sub(r, r, q[i]);
				}
			}
		}

		/* Convert r to affine coordinates. */
		ep4_norm(r, r);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(u);
		for (int i = 0; i < 8; i++) {
			bn_free(_k[i]);
			ep4_free(q[i]);
		}

	}
}

#endif /* EP_MUL == LWNAF */

#if EP_MUL == LWREG || !defined(STRIP)

static void ep4_mul_reg_gls(ep4_t r, const ep4_t p, const bn_t k) {
	int8_t reg[8][RLC_FP_BITS + 1], b[8], s[8], c0, n0;
	ep4_t q, w, t[8][1 << (RLC_WIDTH - 2)];
	bn_t n, _k[8], u;
	size_t l, len, _l[8];

	bn_null(n);
	bn_null(u);
	ep4_null(q);
	ep4_null(w);

	RLC_TRY {
		bn_new(n);
		bn_new(u);
		ep4_new(q);
		ep4_new(w);
		for (size_t i = 0; i < 8; i++) {
			bn_null(_k[i]);
			bn_new(_k[i]);
			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				ep4_null(t[i][j]);
				ep4_new(t[i][j]);
			}
		}

		ep4_curve_get_ord(n);
		fp_prime_get_par(u);
		bn_mod(_k[0], k, n);
		bn_rec_frb(_k, 8, _k[0], u, n, ep_curve_is_pairf() == EP_BN);

		l = 0;
		/* Make some extra room for BN curves that grow subscalars by 1. */
		len = bn_bits(u) + (ep_curve_is_pairf() == EP_BN);
		ep4_norm(t[0][0], p);
		for (size_t i = 0; i < 8; i++) {
			s[i] = bn_sign(_k[i]);
			bn_abs(_k[i], _k[i]);
			b[i] = bn_is_even(_k[i]);
			_k[i]->dp[0] |= b[i];

			_l[i] = RLC_FP_BITS + 1;
			bn_rec_reg(reg[i], &_l[i], _k[i], len, RLC_WIDTH);
			l = RLC_MAX(l, _l[i]);
			
			/* Apply Frobenius before flipping sign to build table. */
			if (i > 0) {
				ep4_psi(t[i][0], t[i - 1][0]);
			}
		}

		for (size_t i = 0; i < 8; i++) {
			ep4_neg(q, t[i][0]);
			fp4_copy_sec(q->y, t[i][0]->y, s[i] == RLC_POS);
			ep4_tab(t[i], q, RLC_WIDTH);
		}

#if defined(EP_MIXED)
		fp4_set_dig(w->z, 1);
		w->coord = BASIC;
#else
		w->coord = = EP_ADD;
#endif

		ep4_set_infty(r);
		for (int j = l - 1; j >= 0; j--) {
			for (size_t i = 0; i < RLC_WIDTH - 1; i++) {
				ep4_dbl(r, r);
			}

			for (size_t i = 0; i < 8; i++) {
				n0 = reg[i][j];
				c0 = (n0 >> 7);
				n0 = ((n0 ^ c0) - c0) >> 1;

				for (size_t m = 0; m < (1 << (RLC_WIDTH - 2)); m++) {
					fp4_copy_sec(w->x, t[i][m]->x, m == n0);
					fp4_copy_sec(w->y, t[i][m]->y, m == n0);
	#if !defined(EP_MIXED)
					fp4_copy_sec(w->z, t[i][m]->z, m == n0);
	#endif
				}

				ep4_neg(q, w);
				fp4_copy_sec(q->y, w->y, c0 == 0);
				ep4_add(r, r, q);
			}
		}

		for (size_t i = 0; i < 8; i++) {
			/* Tables are built with points already negated, so no need here. */
			ep4_sub(q, r, t[i][0]);
			fp4_copy_sec(r->x, q->x, b[i]);
			fp4_copy_sec(r->y, q->y, b[i]);
			fp4_copy_sec(r->z, q->z, b[i]);
		}

		/* Convert r to affine coordinates. */
		ep4_norm(r, r);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(u);
		ep4_free(q);
		ep4_free(w);
		for (int i = 0; i < 4; i++) {
			bn_free(_k[i]);
			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				ep4_free(t[i][j]);
			}
		}
	}
}

#endif /* EP_MUL == LWREG */
#endif /* EP_ENDOM */

#if defined(EP_PLAIN) || defined(EP_SUPER)

#if EP_MUL == LWNAF || !defined(STRIP)

static void ep4_mul_naf_imp(ep4_t r, const ep4_t p, const bn_t k) {
	int i, n;
	int8_t naf[RLC_FP_BITS + 1];
	ep4_t t[1 << (RLC_WIDTH - 2)];
	size_t l;

	RLC_TRY {
		/* Prepare the precomputation table. */
		for (i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep4_null(t[i]);
			ep4_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep4_tab(t, p, RLC_WIDTH);

		/* Compute the w-NAF representation of k. */
		l = sizeof(naf);
		bn_rec_naf(naf, &l, k, RLC_WIDTH);

		ep4_set_infty(r);
		for (i = l - 1; i >= 0; i--) {
			ep4_dbl(r, r);

			n = naf[i];
			if (n > 0) {
				ep4_add(r, r, t[n / 2]);
			}
			if (n < 0) {
				ep4_sub(r, r, t[-n / 2]);
			}
		}
		/* Convert r to affine coordinates. */
		ep4_norm(r, r);
		if (bn_sign(k) == RLC_NEG) {
			ep4_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		/* Free the precomputation table. */
		for (i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep4_free(t[i]);
		}
	}
}

#endif /* EP_MUL == LWNAF */

#if EP_MUL == LWREG || !defined(STRIP)

static void ep4_mul_reg_imp(ep4_t r, const ep4_t p, const bn_t k) {
	bn_t _k;
	int8_t s, reg[1 + RLC_CEIL(RLC_FP_BITS + 1, RLC_WIDTH - 1)];
	ep4_t t[1 << (RLC_WIDTH - 2)], u, v;
	size_t l, n;

	bn_null(_k);

	RLC_TRY {
		bn_new(_k);
		ep4_new(u);
		ep4_new(v);
		/* Prepare the precomputation table. */
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep4_null(t[i]);
			ep4_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep4_tab(t, p, RLC_WIDTH);

		ep4_curve_get_ord(_k);
		n = bn_bits(_k);

		/* Make a copy of the scalar for processing. */
		bn_abs(_k, k);
		_k->dp[0] |= 1;

		/* Compute the regular w-NAF representation of k. */
		l = RLC_CEIL(n, RLC_WIDTH - 1) + 1;
		bn_rec_reg(reg, &l, _k, n, RLC_WIDTH);

#if defined(EP_MIXED)
		fp4_set_dig(u->z, 1);
		u->coord = BASIC;
#else
		u->coord = EP_ADD;
#endif
		ep4_set_infty(r);
		for (int i = l - 1; i >= 0; i--) {
			for (size_t j = 0; j < RLC_WIDTH - 1; j++) {
				ep4_dbl(r, r);
			}

			n = reg[i];
			s = (n >> 7);
			n = ((n ^ s) - s) >> 1;

			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				fp4_copy_sec(u->x, t[j]->x, j == n);
				fp4_copy_sec(u->y, t[j]->y, j == n);
#if !defined(EP_MIXED)
				fp_copy_sec(u->z, t[j]->z, j == n);
#endif
			}
			ep4_neg(v, u);
			fp4_copy_sec(u->y, v->y, s != 0);
			ep4_add(r, r, u);
		}
		/* t[0] has an unmodified copy of p. */
		ep4_sub(u, r, t[0]);
		fp4_copy_sec(r->x, u->x, bn_is_even(k));
		fp4_copy_sec(r->y, u->y, bn_is_even(k));
		fp4_copy_sec(r->z, u->z, bn_is_even(k));
		/* Convert r to affine coordinates. */
		ep4_norm(r, r);
		ep4_neg(u, r);
		fp4_copy_sec(r->y, u->y, bn_sign(k) == RLC_NEG);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		/* Free the precomputation table. */
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep4_free(t[i]);
		}
		bn_free(_k);
		ep4_free(u);
		ep4_free(v);
	}
}

#endif /* EP_MUL == LWREG */
#endif /* EP_PLAIN || EP_SUPER */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep4_mul_basic(ep4_t r, const ep4_t p, const bn_t k) {
	ep4_t t;
	int8_t u, *naf = RLC_ALLOCA(int8_t, bn_bits(k) + 1);
	size_t l;

	ep4_null(t);

	if (bn_is_zero(k) || ep4_is_infty(p)) {
		RLC_FREE(naf);
		ep4_set_infty(r);
		return;
	}

	if (bn_bits(k) <= RLC_DIG) {
		ep4_mul_dig(r, p, k->dp[0]);
		if (bn_sign(k) == RLC_NEG) {
			ep4_neg(r, r);
		}
		RLC_FREE(naf);
		return;
	}

	RLC_TRY {
		ep4_new(t);
		if (naf == NULL) {
			RLC_THROW(ERR_NO_BUFFER);
		}

		l = bn_bits(k) + 1;
		bn_rec_naf(naf, &l, k, 2);
		ep4_set_infty(t);
		for (int i = l - 1; i >= 0; i--) {
			ep4_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep4_add(t, t, p);
			} else if (u < 0) {
				ep4_sub(t, t, p);
			}
		}

		ep4_norm(r, t);
		if (bn_sign(k) == RLC_NEG) {
			ep4_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep4_free(t);
		RLC_FREE(naf);
	}
}

#if EP_MUL == SLIDE || !defined(STRIP)

void ep4_mul_slide(ep4_t r, const ep4_t p, const bn_t k) {
	ep4_t t[1 << (RLC_WIDTH - 1)], q;
	uint8_t win[RLC_FP_BITS + 1];
	size_t l;

	ep4_null(q);

	if (bn_is_zero(k) || ep4_is_infty(p)) {
		ep4_set_infty(r);
		return;
	}

	RLC_TRY {
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 1)); i ++) {
			ep4_null(t[i]);
			ep4_new(t[i]);
		}

		ep4_new(q);

		ep4_copy(t[0], p);
		ep4_dbl(q, p);

#if defined(EP_MIXED)
		ep4_norm(q, q);
#endif

		/* Create table. */
		for (size_t i = 1; i < (1 << (RLC_WIDTH - 1)); i++) {
			ep4_add(t[i], t[i - 1], q);
		}

#if defined(EP_MIXED)
		ep4_norm_sim(t + 1, t + 1, (1 << (RLC_WIDTH - 1)) - 1);
#endif

		ep4_set_infty(q);
		l = RLC_FP_BITS + 1;
		bn_rec_slw(win, &l, k, RLC_WIDTH);
		for (size_t i = 0; i < l; i++) {
			if (win[i] == 0) {
				ep4_dbl(q, q);
			} else {
				for (size_t j = 0; j < util_bits_dig(win[i]); j++) {
					ep4_dbl(q, q);
				}
				ep4_add(q, q, t[win[i] >> 1]);
			}
		}

		ep4_norm(r, q);
		if (bn_sign(k) == RLC_NEG) {
			ep4_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 1)); i++) {
			ep4_free(t[i]);
		}
		ep4_free(q);
	}
}

#endif

#if EP_MUL == MONTY || !defined(STRIP)

void ep4_mul_monty(ep4_t r, const ep4_t p, const bn_t k) {
	ep4_t t[2];
	bn_t n, l, _k;
	size_t bits;

	bn_null(n);
	bn_null(l);
	bn_null(_k);
	ep4_null(t[0]);
	ep4_null(t[1]);

	if (bn_is_zero(k) || ep4_is_infty(p)) {
		ep4_set_infty(r);
		return;
	}

	RLC_TRY {
		bn_new(n);
		bn_new(l);
		bn_new(_k);
		ep4_new(t[0]);
		ep4_new(t[1]);

		ep4_curve_get_ord(n);
		bits = bn_bits(n);

		bn_mod(_k, k, n);
		bn_abs(l, _k);
		bn_add(l, l, n);
		bn_add(n, l, n);
		dv_swap_sec(l->dp, n->dp, RLC_MAX(l->used, n->used),
			bn_get_bit(l, bits) == 0);
		l->used = RLC_SEL(l->used, n->used, bn_get_bit(l, bits) == 0);

		ep4_norm(t[0], p);
		ep4_dbl(t[1], t[0]);

		/* Blind both points independently. */
		ep4_blind(t[0], t[0]);
		ep4_blind(t[1], t[1]);

		for (int i = bits - 1; i >= 0; i--) {
			int j = bn_get_bit(l, i);
			dv_swap_sec(t[0]->x[0][0], t[1]->x[0][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[0][1], t[1]->x[0][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[1][0], t[1]->x[1][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[1][1], t[1]->x[1][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[0][0], t[1]->y[0][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[0][1], t[1]->y[0][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[1][0], t[1]->y[1][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[1][1], t[1]->y[1][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[0][0], t[1]->z[0][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[0][1], t[1]->z[0][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[1][0], t[1]->z[1][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[1][1], t[1]->z[1][1], RLC_FP_DIGS, j ^ 1);
			ep4_add(t[0], t[0], t[1]);
			ep4_dbl(t[1], t[1]);
			dv_swap_sec(t[0]->x[0][0], t[1]->x[0][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[0][1], t[1]->x[0][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[1][0], t[1]->x[1][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->x[1][1], t[1]->x[1][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[0][0], t[1]->y[0][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[0][1], t[1]->y[0][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[1][0], t[1]->y[1][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->y[1][1], t[1]->y[1][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[0][0], t[1]->z[0][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[0][1], t[1]->z[0][1], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[1][0], t[1]->z[1][0], RLC_FP_DIGS, j ^ 1);
			dv_swap_sec(t[0]->z[1][1], t[1]->z[1][1], RLC_FP_DIGS, j ^ 1);
		}

		ep4_norm(r, t[0]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(l);
		bn_free(_k);
		ep4_free(t[1]);
		ep4_free(t[0]);
	}
}

#endif

#if EP_MUL == LWNAF || !defined(STRIP)

void ep4_mul_lwnaf(ep4_t r, const ep4_t p, const bn_t k) {
	if (bn_is_zero(k) || ep4_is_infty(p)) {
		ep4_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep4_mul_glv_imp(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep4_mul_naf_imp(r, p, k);
#endif
}

#endif

#if EP_MUL == LWREG || !defined(STRIP)

void ep4_mul_lwreg(ep4_t r, const ep4_t p, const bn_t k) {
	if (bn_is_zero(k) || ep4_is_infty(p)) {
		ep4_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep4_mul_reg_gls(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep4_mul_reg_imp(r, p, k);
#endif
}

#endif

void ep4_mul_gen(ep4_t r, const bn_t k) {
	if (bn_is_zero(k)) {
		ep4_set_infty(r);
		return;
	}

#ifdef EP_PRECO
	ep4_mul_fix(r, ep4_curve_get_tab(), k);
#else
	ep4_t g;

	ep4_null(g);

	RLC_TRY {
		ep4_new(g);
		ep4_curve_get_gen(g);
		ep4_mul(r, g, k);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep4_free(g);
	}
#endif
}

void ep4_mul_dig(ep4_t r, const ep4_t p, const dig_t k) {
	ep4_t t;
	bn_t _k;
	int8_t u, naf[RLC_DIG + 1];
	size_t l;

	ep4_null(t);
	bn_null(_k);

	if (k == 0 || ep4_is_infty(p)) {
		ep4_set_infty(r);
		return;
	}

	RLC_TRY {
		ep4_new(t);
		bn_new(_k);

		bn_set_dig(_k, k);

		l = RLC_DIG + 1;
		bn_rec_naf(naf, &l, _k, 2);

		ep4_copy(t, p);
		for (int i = l - 2; i >= 0; i--) {
			ep4_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep4_add(t, t, p);
			} else if (u < 0) {
				ep4_sub(t, t, p);
			}
		}

		ep4_norm(r, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep4_free(t);
		bn_free(_k);
	}
}