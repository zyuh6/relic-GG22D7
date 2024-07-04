/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (c) 2023 RELIC Authors
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
 * Implementation of point multiplication on prime elliptic curves over an
 * octic extension field.
 *
 * @ingroup epx
 */

#include "relic_core.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

#if defined(EP_ENDOM)

#if EP_MUL == LWNAF || !defined(STRIP)

static void ep8_mul_glv_imp(ep8_t r, const ep8_t p, const bn_t k) {
	size_t l, _l[16];
	bn_t n, _k[16], u;
	int8_t naf[16][RLC_FP_BITS + 1];
	ep8_t q[16];

	bn_null(n);
	bn_null(u);

	RLC_TRY {
		bn_new(n);
		bn_new(u);
		for (int i = 0; i < 16; i++) {
			bn_null(_k[i]);
			ep8_null(q[i]);
			bn_new(_k[i]);
			ep8_new(q[i]);
		}

		ep8_curve_get_ord(n);
		fp_prime_get_par(u);
		bn_mod(_k[0], k, n);
		bn_rec_frb(_k, 16, _k[0], u, n, ep_curve_is_pairf() == EP_BN);

		ep8_norm(q[0], p);
		for (size_t i = 1; i < 16; i++) {
            ep8_frb(q[i], q[i - 1], 1);
		}

		l = 0;
		for (size_t i = 0; i < 16; i++) {
			if (bn_sign(_k[i]) == RLC_NEG) {
				ep8_neg(q[i], q[i]);
			}
			_l[i] = RLC_FP_BITS + 1;
			bn_rec_naf(naf[i], &_l[i], _k[i], 2);
			l = RLC_MAX(l, _l[i]);
		}

		ep8_set_infty(r);
		for (int j = l - 1; j >= 0; j--) {
			ep8_dbl(r, r);

			for (int i = 0; i < 16; i++) {
				if (naf[i][j] > 0) {
					ep8_add(r, r, q[i]);
				}
				if (naf[i][j] < 0) {
					ep8_sub(r, r, q[i]);
				}
			}
		}

		/* Convert r to affine coordinates. */
		ep8_norm(r, r);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(u);
		for (int i = 0; i < 16; i++) {
			bn_free(_k[i]);
			ep8_free(q[i]);
		}
	}
}

#endif /* EP_MUL == LWNAF */

#if EP_MUL == LWREG || !defined(STRIP)

static void ep8_mul_reg_gls(ep8_t r, const ep8_t p, const bn_t k) {
	int8_t reg[16][RLC_FP_BITS + 1], b[16], s[16], c0, n0;
	ep8_t q, w, t[16][1 << (RLC_WIDTH - 2)];
	bn_t n, _k[16], u;
	size_t l, len, _l[16];

	bn_null(n);
	bn_null(u);
	ep8_null(q);
	ep8_null(w);

	RLC_TRY {
		bn_new(n);
		bn_new(u);
		ep8_new(q);
		ep8_new(w);
		for (size_t i = 0; i < 16; i++) {
			bn_null(_k[i]);
			bn_new(_k[i]);
			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				ep8_null(t[i][j]);
				ep8_new(t[i][j]);
			}
		}

		ep8_curve_get_ord(n);
		fp_prime_get_par(u);
		bn_mod(_k[0], k, n);
		bn_rec_frb(_k, 16, _k[0], u, n, ep_curve_is_pairf() == EP_BN);

		l = 0;
		/* Make some extra room for BN curves that grow subscalars by 1. */
		len = bn_bits(u) + (ep_curve_is_pairf() == EP_BN);
		ep8_norm(t[0][0], p);
		for (size_t i = 0; i < 16; i++) {
			s[i] = bn_sign(_k[i]);
			bn_abs(_k[i], _k[i]);
			b[i] = bn_is_even(_k[i]);
			_k[i]->dp[0] |= b[i];

			_l[i] = RLC_FP_BITS + 1;
			bn_rec_reg(reg[i], &_l[i], _k[i], len, RLC_WIDTH);
			l = RLC_MAX(l, _l[i]);
			
			/* Apply Frobenius before flipping sign to build table. */
			if (i > 0) {
				ep8_frb(t[i][0], t[i - 1][0], 1);
			}
		}

		for (size_t i = 0; i < 16; i++) {
			ep8_neg(q, t[i][0]);
			fp8_copy_sec(q->y, t[i][0]->y, s[i] == RLC_POS);
			ep8_tab(t[i], q, RLC_WIDTH);
		}

#if defined(EP_MIXED)
		fp8_set_dig(w->z, 1);
		w->coord = BASIC;
#else
		w->coord = = EP_ADD;
#endif

		ep8_set_infty(r);
		for (int j = l - 1; j >= 0; j--) {
			for (size_t i = 0; i < RLC_WIDTH - 1; i++) {
				ep8_dbl(r, r);
			}

			for (size_t i = 0; i < 16; i++) {
				n0 = reg[i][j];
				c0 = (n0 >> 7);
				n0 = ((n0 ^ c0) - c0) >> 1;

				for (size_t m = 0; m < (1 << (RLC_WIDTH - 2)); m++) {
					fp8_copy_sec(w->x, t[i][m]->x, m == n0);
					fp8_copy_sec(w->y, t[i][m]->y, m == n0);
	#if !defined(EP_MIXED)
					fp8_copy_sec(w->z, t[i][m]->z, m == n0);
	#endif
				}

				ep8_neg(q, w);
				fp8_copy_sec(q->y, w->y, c0 == 0);
				ep8_add(r, r, q);
			}
		}

		for (size_t i = 0; i < 16; i++) {
			/* Tables are built with points already negated, so no need here. */
			ep8_sub(q, r, t[i][0]);
			fp8_copy_sec(r->x, q->x, b[i]);
			fp8_copy_sec(r->y, q->y, b[i]);
			fp8_copy_sec(r->z, q->z, b[i]);
		}

		/* Convert r to affine coordinates. */
		ep8_norm(r, r);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(u);
		ep8_free(q);
		ep8_free(w);
		for (int i = 0; i < 16; i++) {
			bn_free(_k[i]);
			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				ep8_free(t[i][j]);
			}
		}
	}
}

#endif /* EP_MUL == LWREG */
#endif /* EP_ENDOM */

#if defined(EP_PLAIN) || defined(EP_SUPER)

#if EP_MUL == LWNAF || !defined(STRIP)

static void ep8_mul_naf_imp(ep8_t r, const ep8_t p, const bn_t k) {
	int i, n;
	int8_t naf[RLC_FP_BITS + 1];
	ep8_t t[1 << (RLC_WIDTH - 2)];
	size_t l;

	RLC_TRY {
		/* Prepare the precomputation table. */
		for (i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep8_null(t[i]);
			ep8_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep8_tab(t, p, RLC_WIDTH);

		/* Compute the w-NAF representation of k. */
		l = sizeof(naf);
		bn_rec_naf(naf, &l, k, RLC_WIDTH);

		ep8_set_infty(r);
		for (i = l - 1; i >= 0; i--) {
			ep8_dbl(r, r);

			n = naf[i];
			if (n > 0) {
				ep8_add(r, r, t[n / 2]);
			}
			if (n < 0) {
				ep8_sub(r, r, t[-n / 2]);
			}
		}
		/* Convert r to affine coordinates. */
		ep8_norm(r, r);
		if (bn_sign(k) == RLC_NEG) {
			ep8_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		/* Free the precomputation table. */
		for (i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep8_free(t[i]);
		}
	}
}

#endif /* EP_MUL == LWNAF */

#if EP_MUL == LWREG || !defined(STRIP)

static void ep8_mul_reg_imp(ep8_t r, const ep8_t p, const bn_t k) {
	bn_t _k;
	int8_t s, reg[1 + RLC_CEIL(RLC_FP_BITS + 1, RLC_WIDTH - 1)];
	ep8_t t[1 << (RLC_WIDTH - 2)], u, v;
	size_t l, n;

	bn_null(_k);

	RLC_TRY {
		bn_new(_k);
		ep8_new(u);
		ep8_new(v);
		/* Prepare the precomputation table. */
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep8_null(t[i]);
			ep8_new(t[i]);
		}
		/* Compute the precomputation table. */
		ep8_tab(t, p, RLC_WIDTH);

		ep8_curve_get_ord(_k);
		n = bn_bits(_k);

		/* Make a copy of the scalar for processing. */
		bn_abs(_k, k);
		_k->dp[0] |= 1;

		/* Compute the regular w-NAF representation of k. */
		l = RLC_CEIL(n, RLC_WIDTH - 1) + 1;
		bn_rec_reg(reg, &l, _k, n, RLC_WIDTH);

#if defined(EP_MIXED)
		fp8_set_dig(u->z, 1);
		u->coord = BASIC;
#else
		u->coord = EP_ADD;
#endif
		ep8_set_infty(r);
		for (int i = l - 1; i >= 0; i--) {
			for (size_t j = 0; j < RLC_WIDTH - 1; j++) {
				ep8_dbl(r, r);
			}

			n = reg[i];
			s = (n >> 7);
			n = ((n ^ s) - s) >> 1;

			for (size_t j = 0; j < (1 << (RLC_WIDTH - 2)); j++) {
				fp8_copy_sec(u->x, t[j]->x, j == n);
				fp8_copy_sec(u->y, t[j]->y, j == n);
#if !defined(EP_MIXED)
				fp_copy_sec(u->z, t[j]->z, j == n);
#endif
			}
			ep8_neg(v, u);
			fp8_copy_sec(u->y, v->y, s != 0);
			ep8_add(r, r, u);
		}
		/* t[0] has an unmodified copy of p. */
		ep8_sub(u, r, t[0]);
		fp8_copy_sec(r->x, u->x, bn_is_even(k));
		fp8_copy_sec(r->y, u->y, bn_is_even(k));
		fp8_copy_sec(r->z, u->z, bn_is_even(k));
		/* Convert r to affine coordinates. */
		ep8_norm(r, r);
		ep8_neg(u, r);
		fp8_copy_sec(r->y, u->y, bn_sign(k) == RLC_NEG);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		/* Free the precomputation table. */
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 2)); i++) {
			ep8_free(t[i]);
		}
		bn_free(_k);
		ep8_free(u);
		ep8_free(v);
	}
}

#endif /* EP_MUL == LWREG */
#endif /* EP_PLAIN || EP_SUPER */

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep8_mul_basic(ep8_t r, const ep8_t p, const bn_t k) {
	ep8_t t;
	int8_t u, *naf = RLC_ALLOCA(int8_t, bn_bits(k) + 1);
	size_t l;

	ep8_null(t);

	if (bn_is_zero(k) || ep8_is_infty(p)) {
		RLC_FREE(naf);
		ep8_set_infty(r);
		return;
	}

	if (bn_bits(k) <= RLC_DIG) {
		ep8_mul_dig(r, p, k->dp[0]);
		if (bn_sign(k) == RLC_NEG) {
			ep8_neg(r, r);
		}
		RLC_FREE(naf);
		return;
	}

	RLC_TRY {
		ep8_new(t);
		if (naf == NULL) {
			RLC_THROW(ERR_NO_BUFFER);
		}

		l = bn_bits(k) + 1;
		bn_rec_naf(naf, &l, k, 2);
		ep8_set_infty(t);
		for (int i = l - 1; i >= 0; i--) {
			ep8_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep8_add(t, t, p);
			} else if (u < 0) {
				ep8_sub(t, t, p);
			}
		}

		ep8_norm(r, t);
		if (bn_sign(k) == RLC_NEG) {
			ep8_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep8_free(t);
		RLC_FREE(naf);
	}
}

#if EP_MUL == SLIDE || !defined(STRIP)

void ep8_mul_slide(ep8_t r, const ep8_t p, const bn_t k) {
	ep8_t t[1 << (RLC_WIDTH - 1)], q;
	uint8_t win[RLC_FP_BITS + 1];
	size_t l;

	ep8_null(q);

	if (bn_is_zero(k) || ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

	RLC_TRY {
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 1)); i ++) {
			ep8_null(t[i]);
			ep8_new(t[i]);
		}

		ep8_new(q);

		ep8_copy(t[0], p);
		ep8_dbl(q, p);

#if defined(EP_MIXED)
		ep8_norm(q, q);
#endif

		/* Create table. */
		for (size_t i = 1; i < (1 << (RLC_WIDTH - 1)); i++) {
			ep8_add(t[i], t[i - 1], q);
		}

#if defined(EP_MIXED)
		ep8_norm_sim(t + 1, t + 1, (1 << (RLC_WIDTH - 1)) - 1);
#endif

		ep8_set_infty(q);
		l = RLC_FP_BITS + 1;
		bn_rec_slw(win, &l, k, RLC_WIDTH);
		for (size_t i = 0; i < l; i++) {
			if (win[i] == 0) {
				ep8_dbl(q, q);
			} else {
				for (size_t j = 0; j < util_bits_dig(win[i]); j++) {
					ep8_dbl(q, q);
				}
				ep8_add(q, q, t[win[i] >> 1]);
			}
		}

		ep8_norm(r, q);
		if (bn_sign(k) == RLC_NEG) {
			ep8_neg(r, r);
		}
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		for (size_t i = 0; i < (1 << (RLC_WIDTH - 1)); i++) {
			ep8_free(t[i]);
		}
		ep8_free(q);
	}
}

#endif

#if EP_MUL == MONTY || !defined(STRIP)

void ep8_mul_monty(ep8_t r, const ep8_t p, const bn_t k) {
	ep8_t t[2];
	bn_t n, l, _k;
	size_t bits;

	bn_null(n);
	bn_null(l);
	bn_null(_k);
	ep8_null(t[0]);
	ep8_null(t[1]);

	if (bn_is_zero(k) || ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

	RLC_TRY {
		bn_new(n);
		bn_new(l);
		bn_new(_k);
		ep8_new(t[0]);
		ep8_new(t[1]);

		ep8_curve_get_ord(n);
		bits = bn_bits(n);

		bn_mod(_k, k, n);
		bn_abs(l, _k);
		bn_add(l, l, n);
		bn_add(n, l, n);
		dv_swap_sec(l->dp, n->dp, RLC_MAX(l->used, n->used),
			bn_get_bit(l, bits) == 0);
		l->used = RLC_SEL(l->used, n->used, bn_get_bit(l, bits) == 0);

		ep8_norm(t[0], p);
		ep8_dbl(t[1], t[0]);

		/* Blind both points independently. */
		ep8_blind(t[0], t[0]);
		ep8_blind(t[1], t[1]);

		for (int i = bits - 1; i >= 0; i--) {
			int j = bn_get_bit(l, i);
			for (int l = 0; l < 2; l++) {
				for (int m = 0; m < 2; m++) {
					for (int n = 0; n < 2; n++) {
						dv_swap_sec(t[0]->x[l][m][n], t[1]->x[l][m][n], RLC_FP_DIGS, j ^ 1);
						dv_swap_sec(t[0]->y[l][m][n], t[1]->y[l][m][n], RLC_FP_DIGS, j ^ 1);
						dv_swap_sec(t[0]->z[l][m][n], t[1]->z[l][m][n], RLC_FP_DIGS, j ^ 1);
					}
				}
			}
			ep8_add(t[0], t[0], t[1]);
			ep8_dbl(t[1], t[1]);
			for (int l = 0; l < 2; l++) {
				for (int m = 0; m < 2; m++) {
					for (int n = 0; n < 2; n++) {
						dv_swap_sec(t[0]->x[l][m][n], t[1]->x[l][m][n], RLC_FP_DIGS, j ^ 1);
						dv_swap_sec(t[0]->y[l][m][n], t[1]->y[l][m][n], RLC_FP_DIGS, j ^ 1);
						dv_swap_sec(t[0]->z[l][m][n], t[1]->z[l][m][n], RLC_FP_DIGS, j ^ 1);
					}
				}
			}
		}

		ep8_norm(r, t[0]);
	} RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		bn_free(n);
		bn_free(l);
		bn_free(_k);
		ep8_free(t[1]);
		ep8_free(t[0]);
	}
}

#endif

#if EP_MUL == LWNAF || !defined(STRIP)

void ep8_mul_lwnaf(ep8_t r, const ep8_t p, const bn_t k) {
	if (bn_is_zero(k) || ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep8_mul_glv_imp(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep8_mul_naf_imp(r, p, k);
#endif
}

#endif

#if EP_MUL == LWREG || !defined(STRIP)

void ep8_mul_lwreg(ep8_t r, const ep8_t p, const bn_t k) {
	if (bn_is_zero(k) || ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

#if defined(EP_ENDOM)
	if (ep_curve_is_endom()) {
		ep8_mul_reg_gls(r, p, k);
		return;
	}
#endif

#if defined(EP_PLAIN) || defined(EP_SUPER)
	ep8_mul_reg_imp(r, p, k);
#endif
}

#endif

void ep8_mul_gen(ep8_t r, const bn_t k) {
	if (bn_is_zero(k)) {
		ep8_set_infty(r);
		return;
	}

#ifdef EP_PRECO
	ep8_mul_fix(r, ep8_curve_get_tab(), k);
#else
	ep8_t g;

	ep8_null(g);

	RLC_TRY {
		ep8_new(g);
		ep8_curve_get_gen(g);
		ep8_mul(r, g, k);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep8_free(g);
	}
#endif
}

void ep8_mul_dig(ep8_t r, const ep8_t p, const dig_t k) {
	ep8_t t;
	bn_t _k;
	int8_t u, naf[RLC_DIG + 1];
	size_t l;

	ep8_null(t);
	bn_null(_k);

	if (k == 0 || ep8_is_infty(p)) {
		ep8_set_infty(r);
		return;
	}

	RLC_TRY {
		ep8_new(t);
		bn_new(_k);

		bn_set_dig(_k, k);

		l = RLC_DIG + 1;
		bn_rec_naf(naf, &l, _k, 2);

		ep8_copy(t, p);
		for (int i = l - 2; i >= 0; i--) {
			ep8_dbl(t, t);

			u = naf[i];
			if (u > 0) {
				ep8_add(t, t, p);
			} else if (u < 0) {
				ep8_sub(t, t, p);
			}
		}

		ep8_norm(r, t);
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep8_free(t);
		bn_free(_k);
	}
}