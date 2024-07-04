/*
 * RELIC is an Efficient LIbrary for Cryptography
 * Copyright (C) 2007-2019 RELIC Authors
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
#include "relic.h"
#include "relic_epx.h"

#if defined(EP_PLAIN) && FP_PRIME == 457

#define GG22D7_P457_A0		"1D6D6FF82B405135B174F9F789277226DA3C03C1A07A8E54143B0911F58F28355F79553C2B3B34CF7CAA8E4581B79E1EDE05A2878B8E12F3E78"
#define GG22D7_P457_A1		"0"
#define GG22D7_P457_A2		"0"
#define GG22D7_P457_A3		"0"
#define GG22D7_P457_A4		"0"
#define GG22D7_P457_A5		"0"
#define GG22D7_P457_A6		"0"
#define GG22D7_P457_A7		"0"
#define GG22D7_P457_A8		"0"
#define GG22D7_P457_A9		"EB6B7FC15A0289AD8BA7CFBC493B9136D1E01E0D03D472A0A1D8488FAC7941AAFBCAA9E159D9A67BE554722C0DBCF0F6F02D143C5C70979F3C"
#define GG22D7_P457_A10		"EB6B7FC15A0289AD8BA7CFBC493B9136D1E01E0D03D472A0A1D8488FAC7941AAFBCAA9E159D9A67BE554722C0DBCF0F6F02D143C5C70979F3F"
#define GG22D7_P457_B0		"1BFB7A6370BC5AE0A257A012BD7C36B0FC179E919956509031CF5C6F1632E369A24333D211C27E514DAAD1F8CB59BD30A42852BD57011ACFF7B"
#define GG22D7_P457_B1		"0"
#define GG22D7_P457_B2		"0"
#define GG22D7_P457_B3		"0"
#define GG22D7_P457_B4		"0"
#define GG22D7_P457_B5		"0"
#define GG22D7_P457_B6		"0"
#define GG22D7_P457_B7		"0"
#define GG22D7_P457_B8		"B8FACA5D41FB2A878EACF265D59DBAEF12329803921EE1F135D6516FAE2265DE9B10B50CBC5B3F177FDE265B2EF0771CEEA7E51A467C11F80"
#define GG22D7_P457_B9		"1CB4752DCDFE560B29E64D052351D46BEB29D1299CE86F72230532C085E105CF80DE44871E7ED990652AB01F2688ADA7C116FAA2714796E1EFB"
#define GG22D7_P457_B10		"B8FACA5D41FB2A878EACF265D59DBAEF12329803921EE1F135D6516FAE2265DE9B10B50CBC5B3F177FDE265B2EF0771CEEA7E51A467C11F80"

#define GG22D7_P457_X0		"FF67F96C3DEC4D0E01D8FC6E3C44726AE5B757FF5BC6BBD90DB50659B1DF7A1DC23CB06DA719702B9A1C59D700B58910A2E664A54BA0E23C2A"
#define GG22D7_P457_X1      "19C90AA58695BC208083657ADDA9D91E09290EE03F29F4ED3DBD582C990497C53DAF944E01CE713BDA69C32602237E652D75C5A19CAE9D5824E"
#define GG22D7_P457_X2		"16B16355DEF41094C94B29D843BE4EAECB42C8101E870D6A9F5C9B350909BE2EC4D09AF644F2A216149599D57C1A36C1900DFF3BD48BB78FF45"
#define GG22D7_P457_X3		"66AA407E4F6226EE61191FC496B72BFD63C3AC2AD9948D2558C6EA45E1EDF7A30A2E9C0ED33A73E7D4E49C7BC5E4DAE537A7B849446EA3C8E7"
#define GG22D7_P457_X4		"E604241BEE3EA898103CA25D03E0D82FE31724C99ECF57B503329AE1D169572A9E9C0810AC8ADAAC15CD271B470532FDCD63CA526F9B52CD94"
#define GG22D7_P457_X5		"185EA69DEAAAC792EF22E8664F85D27F8681E39AD146E8A9FFC7B1322B306EF5DCA6085632BF14F4ED2A2F165417EC3B3E79E536CCF084B0AE1"
#define GG22D7_P457_X6		"1BAE0A085016C81F6D730DA492200299271E2AD2B575A683C36C894C6D47C8A0D1EED43FE45B6B20C73DACF4E5BC0DFEAD83E341046F46E92E8"
#define GG22D7_P457_X7		"64008108F1BE73A284C01AC8301CB5E804A2460739D0611AA2FA8DF3B228B06D3454E7AF5ECF99ABD310459AFDE62FE3DDBF007FAECF39B544"
#define GG22D7_P457_X8		"74F3CCA01E35B5A8358F438A7AA89D8B665396F23BCFA28186B65F7003D3E23CD991FE5CC448AF033E306C39E568D72B04B37493B3D80376F"
#define GG22D7_P457_X9		"901E0EC029591D1D6A59DC044F3140D43A94FF7BADC28E416C196F5982D399B5F0A5ACCD5756F95E5DA2DFE46BA059D4CC5B07E82B771CBB49"
#define GG22D7_P457_X10		"13BEF21F5224B0F4F7E2C0215EDEC75216B61BBA358050D60A9177EE3733266FDD7BC5270E6B8D2B3341C4AFCDCF38BBC9FA3886729DBF432C4"

#define GG22D7_P457_Y0		"8E1596BDE449FAB4855B3A95F0A63BE681C1E6FB17830F4EF25CAB6D58F489CFECEF12EF460600AC6E6D0FFFE4A671C0C06A4038A2FAD0EC9"
#define GG22D7_P457_Y1		"1C09A06C03629C9CC34BB2629BD5CD44F0638ECF2802DB95CDC9D90E271ECAD98045A43049BF677ECC453C08C4CADBDDF06AF2999C96E76834"
#define GG22D7_P457_Y2		"184B4B6717CFE15347B7A78ADE3C1E759D2D056DA92AE2571BFF651029CAB49C85BD8D032D461FBD5C90B9349B0EEE4BD475C0C9BE095918D9"
#define GG22D7_P457_Y3		"1D46119B751A3F6EDC8874F5BA2AAE0D6E8F63B53818F23F490989DD092D2FF18AE665A8D80F3E4BFDAEB8D5870B66A713DCE31AE79517CADDB"
#define GG22D7_P457_Y4		"1931890DE051993C2701F1C51E6CB60049005C673DA606DBB9BE20E8165ED55753BF6CF2AD9C50865F5F25421765733C47D3C589C5D74649063"
#define GG22D7_P457_Y5		"162D2D73DC2D491CE75513A28B4633C191101371754F46036F052851C38F41E909FCBE145E348896974BE7D2A0D76210D5D413B274CBE970260"
#define GG22D7_P457_Y6		"3F1004D34F5DF3F91BE41A694661CBA782BBF72B02016E7A48351854868FB16BD40297BF8C8F5DD07C8507263DA2363B74FC150185E37F531C"
#define GG22D7_P457_Y7		"757FE81CB0A9D54829DB9609755794834B7900A74A72047DA297D21ED9F85FDBA4D8FBC5B1A3D8E2537ABC105F3BF980CA6479894420668A88"
#define GG22D7_P457_Y8		"102F97374A252E987617D497109C1019A981717FF8DA66C9FC1A1E91D13817A80596DBBC97F2638B81445E276EB224942EE90FF0AA3BDF917EF"
#define GG22D7_P457_Y9		"56A31434210BCCAFFCCDB995D4C87AF7157B13EB1D4D4E16497CBDC5DC73B5CBC0BDA044845887BAA0EAA733B47DB100BAC824F7995028AB14"
#define GG22D7_P457_Y10		"113DD091744E9D81A90297981B12EDDE06841D5EF95BDE6ECD8CB869A0E0EC17339CBFCB2438832F654FE5ED33B5BCDF103826FB413CC9836A"

#define GG22D7_P457_R		"544DACE03351E476C00D591F5A370FC93D2C9666C22874C2514C2EE109931EF6C130BB95A1DF04DD7764B0A22F1F0FDD"
#define GG22D7_P457_H		"1"

#define GG22D7_P457_T2		"5DEC4F6E90031B9C57A71C175CC621F5418AF0165F870E41D42DD34499D431B045F4A0CEAF3DF48B29974F7A613F6BE3E0AFE8F8473D74ECFE"

#endif

#define ASSIGN(CURVE)   \
	RLC_GET(str, CURVE##_A0, sizeof(CURVE##_A0));							\
	fp_read_str(a[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A1, sizeof(CURVE##_A1));							\
	fp_read_str(a[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A2, sizeof(CURVE##_A2));							\
	fp_read_str(a[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A3, sizeof(CURVE##_A3));							\
	fp_read_str(a[3], str, strlen(str), 16);		                        \
	RLC_GET(str, CURVE##_A4, sizeof(CURVE##_A4));							\
	fp_read_str(a[4], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A5, sizeof(CURVE##_A5));							\
	fp_read_str(a[5], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A6, sizeof(CURVE##_A6));							\
	fp_read_str(a[6], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A7, sizeof(CURVE##_A7));							\
	fp_read_str(a[7], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A8, sizeof(CURVE##_A8));							\
	fp_read_str(a[8], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A9, sizeof(CURVE##_A9));							\
	fp_read_str(a[9], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_A10, sizeof(CURVE##_A10));							\
	fp_read_str(a[10], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B0, sizeof(CURVE##_B0));							\
	fp_read_str(b[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B1, sizeof(CURVE##_B1));							\
	fp_read_str(b[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B2, sizeof(CURVE##_B2));							\
	fp_read_str(b[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B3, sizeof(CURVE##_B3));							\
	fp_read_str(b[3], str, strlen(str), 16);	                         	\
	RLC_GET(str, CURVE##_B4, sizeof(CURVE##_B4));							\
	fp_read_str(b[4], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B5, sizeof(CURVE##_B5));							\
	fp_read_str(b[5], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_B6, sizeof(CURVE##_B6));                           \
	fp_read_str(b[6], str, strlen(str), 16);	                            \
	RLC_GET(str, CURVE##_B7, sizeof(CURVE##_B7));                           \
	fp_read_str(b[7], str, strlen(str), 16);	                            \
	RLC_GET(str, CURVE##_B8, sizeof(CURVE##_B8));                           \
	fp_read_str(b[8], str, strlen(str), 16);	                            \
	RLC_GET(str, CURVE##_B9, sizeof(CURVE##_B9));                           \
	fp_read_str(b[9], str, strlen(str), 16);	                            \
	RLC_GET(str, CURVE##_B10, sizeof(CURVE##_B10));                           \
	fp_read_str(b[10], str, strlen(str), 16);	                            \
    RLC_GET(str, CURVE##_X0, sizeof(CURVE##_X0));								\
    fp_read_str(g->x[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X1, sizeof(CURVE##_X1));								\
	fp_read_str(g->x[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X2, sizeof(CURVE##_X2));								\
	fp_read_str(g->x[2], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X3, sizeof(CURVE##_X3));								\
	fp_read_str(g->x[3], str, strlen(str), 16);								\
    RLC_GET(str, CURVE##_X4, sizeof(CURVE##_X4));								\
    fp_read_str(g->x[4], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X5, sizeof(CURVE##_X5));								\
    fp_read_str(g->x[5], str, strlen(str), 16);								\
    RLC_GET(str, CURVE##_X6, sizeof(CURVE##_X6));								\
	fp_read_str(g->x[6], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X7, sizeof(CURVE##_X7));								\
	fp_read_str(g->x[7], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X8, sizeof(CURVE##_X8));								\
	fp_read_str(g->x[8], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X9, sizeof(CURVE##_X9));								\
	fp_read_str(g->x[9], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_X10, sizeof(CURVE##_X10));								\
	fp_read_str(g->x[10], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y0, sizeof(CURVE##_Y0));								\
	fp_read_str(g->y[0], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y1, sizeof(CURVE##_Y1));								\
	fp_read_str(g->y[1], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y2, sizeof(CURVE##_Y2));								\
	fp_read_str(g->y[2], str, strlen(str), 16);								\
    RLC_GET(str, CURVE##_Y3, sizeof(CURVE##_Y3));								\
	fp_read_str(g->y[3], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y4, sizeof(CURVE##_Y4));								\
	fp_read_str(g->y[4], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y5, sizeof(CURVE##_Y5));								\
	fp_read_str(g->y[5], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y6, sizeof(CURVE##_Y6));					        \
	fp_read_str(g->y[6], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y7, sizeof(CURVE##_Y7));								\
	fp_read_str(g->y[7], str, strlen(str), 16);								\
 	RLC_GET(str, CURVE##_Y8, sizeof(CURVE##_Y8));							\
 	fp_read_str(g->y[8], str, strlen(str), 16);								\
 	RLC_GET(str, CURVE##_Y9, sizeof(CURVE##_Y9));							\
	fp_read_str(g->y[9], str, strlen(str), 16);								\
	RLC_GET(str, CURVE##_Y10, sizeof(CURVE##_Y10));							\
	fp_read_str(g->y[10], str, strlen(str), 16);							\
	RLC_GET(str, CURVE##_T2, sizeof(CURVE##_T2));							\
	fp_read_str(T2, str, strlen(str), 16);

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void ep11_curve_set_twist(int type){
    char str[2 * RLC_FP_BYTES + 1];
	ctx_t *ctx = core_get();
	fp11_t a,b;
	ep11_t g;
	fp_t T2;
	fp11_null(a);
	fp11_null(b);
	fp_null(T2);
	ep11_null(g);
	ctx->ep11_is_twist = 0;	
	if (type == RLC_EP_MTYPE || type == RLC_EP_DTYPE) {
		ctx->ep11_is_twist = type;
	} else {
		return;
	}		

	RLC_TRY {
		fp11_new(a);
		fp11_new(b);
		ep11_new(g);
		fp_new(T2);

		ASSIGN(GG22D7_P457);
		fp11_field_init();
		fp11_zero(g->z);
		fp_set_dig(g->z[0], 1);
		g->coord = BASIC;
		ep11_copy(&(ctx->ep11_g), g);
		fp11_copy(ctx->ep11_a, a);
		fp11_copy(ctx->ep11_b, b);
		fp_copy(ctx->ep11_T2, T2);
		
	}
	RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	}
	RLC_FINALLY {
		ep11_free(g);
		fp11_free(a);
		fp11_free(b);

	}
}				

void ep11_curve_get_a(fp11_t a) {
	fp11_copy(a, core_get()->ep11_a);
}

void ep11_curve_get_b(fp11_t b) {
	fp11_copy(b, core_get()->ep11_b);
}
				
void ep11_curve_get_gen(ep11_t g) {
	ep11_copy(g, &(core_get()->ep11_g));
}

void ep11_curve_get_T2x(fp_t T2) {
	fp_copy(T2, &(core_get()->ep11_T2));
}

void ep11_curve_get_half(fp_t h) {
	fp_copy(h, &(core_get()->fp11_half));
}

int ep11_curve_is_twist(void) {
	return core_get()->ep11_is_twist;
}

