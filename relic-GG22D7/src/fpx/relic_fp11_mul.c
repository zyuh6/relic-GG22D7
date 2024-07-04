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
 * Implementation of multiplication in an undecic extension of a prime field.
 *
 * @ingroup fpx
 */

#include"relic_core.h"
#include"relic_fp.h"
#include"relic_fp_low.h"
#include"relic_fpx_low.h"
/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

/**
 * Given a=a[0]+a[1]x and  b=b[0]+ b[1]x. compute a*b
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_mul_level2(dv_t* c, const fp_t* a, const fp_t* b) {
      dig_t k0[2 * RLC_FP_DIGS], k1[2 * RLC_FP_DIGS];

      fp_muln_low(c[0], a[0], b[0]);
      fp_muln_low(c[2], a[1], b[1]);

      fp_addn_low(k0, a[0], a[1]);
      fp_addn_low(k1, b[0], b[1]);
      fp_muln_low(c[1], k0, k1);
      fp_subd_low(c[1], c[1], c[0]);
      fp_subd_low(c[1], c[1], c[2]);
}

/**
 * Given a=a[0]+a[1]x+a[2]x^2 and  b=b[0]+ b[1]x+b[2]x^2. compute a*b
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_mul_level3(dv_t* c, const fp_t* a, const fp_t* b) {
      dig_t k0[2 * RLC_FP_DIGS], k1[2 * RLC_FP_DIGS], k2[2 * RLC_FP_DIGS];

      for(int i = 0; i<3;i++) fp_muln_low(c[2*i],a[i],b[i]);

      fp_addn_low(k0, a[0], a[1]);
      fp_addn_low(k1, b[0], b[1]);
      fp_muln_low(c[1], k0, k1);
      fp_subd_low(c[1], c[1], c[0]);
      fp_subd_low(c[1], c[1], c[2]);

      fp_addn_low(k0, a[1], a[2]);
      fp_addn_low(k1, b[1], b[2]);
      fp_muln_low(c[3], k0, k1);
      fp_subd_low(c[3], c[3], c[2]);
      fp_subd_low(c[3], c[3], c[4]);

      fp_addn_low(k0, a[0], a[2]);
      fp_addn_low(k1, b[0], b[2]);
      fp_muln_low(k2, k0, k1);
      fp_addd_low(c[2], c[2], k2);
      fp_subd_low(c[2], c[2], c[0]);
      fp_subd_low(c[2], c[2], c[4]);

}

/**
 * Given a=a[0]+a[1]x+...+a[4]x^4 and  b=b[0]+ b[1]x+...+b[4]x^4.
 *  compute a*b
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_mul_level5(dv_t* c, const fp_t* a, const fp_t* b) {

      fp_t u[3];
      fp_t uu[3];
      dv_t m[3],mm[5],mmm[5];  

      for(int i =0;i<3;i++){
            fp_null(u[i]);
            fp_null(uu[i]);
            dv_null(m[i]);
            dv_null(mm[i]);
            dv_null(mmm[i]);
      }
      dv_null(mm[3]);
      dv_null(mm[4]);
      dv_null(mmm[3]);
      dv_null(mmm[4]);
      

      RLC_TRY{
            for(int i =0;i<3;i++){
                  fp_new(u[i]);
                  fp_new(uu[i]);
                  dv_new(m[i]);
                  dv_new(mm[i]);
                  dv_new(mmm[i]);
            }
            dv_new(mm[3]);
            dv_new(mm[4]);
            dv_new(mmm[3]);
            dv_new(mmm[4]);

            for(int i = 0;i<2; i++){
                  fp_addn_low(u[i],a[i],a[i+2]);
                  fp_addn_low(uu[i],b[i],b[i+2]);
            }
            fp_copy(u[2],a[4]);
            fp_copy(uu[2],b[4]);
            fp_mul_level2(m,a,b);
            fp_mul_level3(mm,a+2,b+2);
            fp_mul_level3(mmm,u,uu);
            for(int i = 0;i<3; i++){
                  fp_subd_low(mmm[i],mmm[i],m[i]);
                  fp_subd_low(mmm[i],mmm[i],mm[i]);
            }
            fp_subd_low(mmm[3],mmm[3],mm[3]);
            fp_subd_low(mmm[4],mmm[4],mm[4]);
            for(int i = 0;i<2; i++){
                  dv_copy(c[i],m[i],16);
            }
            fp_addd_low(c[2],m[2],mmm[0]);
            dv_copy(c[3],mmm[1],16);
            for(int i = 0;i<3;i++){
                  fp_addd_low(c[i+4],mm[i],mmm[i+2]);
            }
            for(int i = 0;i<2; i++){
                  dv_copy(c[7+i],mm[i+3],16);
            }   
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            for(int i = 0; i <3; i++)dv_free(m[i]);
            for(int i = 0; i <5; i++){
                  dv_free(mm[i]);	
                  dv_free(mmm[i]);	
            }
            for(int i = 0; i < 3; i++){
                  fp_free(u[i]);
                  fp_free(uu[i]);	
            }
	}


}

/**
 * Given a=a[0]+a[1]x+...+a[5]x^5 and  b=b[0]+ b[1]x+...+b[5]x^5.
 *  compute a*b
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_mul_level6(dv_t* c, const fp_t* a, const fp_t* b) {
      fp_t u[3];
      fp_t uu[3];
      dv_t m[5],mm[5],mmm[5];  

      for(int i =0;i<3;i++){
            fp_null(u[i]);
            fp_null(uu[i]);
            dv_null(m[i]);
            dv_null(mm[i]);
            dv_null(mmm[i]);
      }
      dv_null(m[3]);
      dv_null(m[4]);
      dv_null(mm[3]);
      dv_null(mm[4]);
      dv_null(mmm[3]);
      dv_null(mmm[4]);

      RLC_TRY{
            for(int i =0;i<3;i++){
                  fp_new(u[i]);
                  fp_new(uu[i]);
                  dv_new(m[i]);
                  dv_new(mm[i]);
                  dv_new(mmm[i]);
            }
            dv_new(m[3]);
            dv_new(m[4]);
            dv_new(mm[3]);
            dv_new(mm[4]);
            dv_new(mmm[3]);
            dv_new(mmm[4]);

            for(int i = 0;i<3; i++){
                  fp_addn_low(u[i],a[i],a[i+3]);
                  fp_addn_low(uu[i],b[i],b[i+3]);
            }
            fp_mul_level3(m,a,b);
            fp_mul_level3(mm,a+3,b+3);
            fp_mul_level3(mmm,u,uu);
            for(int i = 0;i<5; i++){
                  fp_subd_low(mmm[i],mmm[i],m[i]);
                  fp_subd_low(mmm[i],mmm[i],mm[i]);
            }
            for(int i = 0;i<3; i++){ 
                  dv_copy(c[i], m[i],16);
            }
            dv_copy(c[5],mmm[2],16);
            for(int i = 0;i<3; i++){ 
                  dv_copy(c[i+8], mm[i+2],16);
            }

            fp_addd_low(c[3],m[3],mmm[0]);
            fp_addd_low(c[4],m[4],mmm[1]);

            fp_addd_low(c[7],mm[1],mmm[4]);
            fp_addd_low(c[6],mm[0],mmm[3]);
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            for(int i = 0; i < 5; i++){
            dv_free(m[i]);
            dv_free(mm[i]);	
            dv_free(mmm[i]);	
            }
            for(int i = 0; i < 3; i++){
            fp_free(u[i]);
            fp_free(uu[i]);	
            }
	}
}

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

void fp11_muln_low(dv11_t c, const fp11_t a, const fp11_t b){
      fp_t u[6],uu[6];
      dv_t m[9],mm[11],mmm[11];	
      for(int i = 0; i <6; i++){
            fp_null(u[i]);
            fp_null(uu[i]);
      }
      for(int i = 0; i <9; i++)dv_null(m[i]);
      for(int i = 0; i <11; i++){
            dv_null(mm[i]);	
            dv_null(mmm[i]);	
      }
      RLC_TRY{
            for(int i = 0; i <6; i++){
                  fp_new(u[i]);
                  fp_new(uu[i]);
            }
            for(int i = 0; i <9; i++)dv_new(m[i]);
            for(int i = 0; i <11; i++){
                  dv_new(mm[i]);	
                  dv_new(mmm[i]);	
            }
            fp_mul_level5(m,a,b);
            fp_mul_level6(mm,a+5,b+5);
            for(int i = 0;i<5; i++){
                  fp_addn_low(u[i],a[i],a[i+5]);
                  fp_addn_low(uu[i],b[i],b[i+5]);
            }       
            fp_copy(u[5],a[10]);
            fp_copy(uu[5],b[10]);
            fp_mul_level6(mmm,u,uu);
            for(int i = 0;i<9; i++){
                  fp_subd_low(mmm[i],mmm[i],m[i]);
                  fp_subd_low(mmm[i],mmm[i],mm[i]);
            }
            fp_subd_low(mmm[9],mmm[9],mm[9]);
            fp_subd_low(mmm[10],mmm[10],mm[10]);

            for(int i = 0;i<5; i++){   
                  fp_addd_low(c[i],mm[i+1],mmm[i+6]);
            }
            fp_addd_low(c[10],c[4],c[4]);//temp c5
            for(int i=4;i>0;i--){
                  fp_addd_low(c[i],c[i],c[i-1]);
                  fp_addd_low(c[i],c[i],c[i]);
                  fp_addd_low(c[i],m[i],c[i]);
            }
            fp_addd_low(c[0],c[0],c[0]);
            fp_addd_low(c[0],m[0],c[0]);

            for(int i = 0;i<4; i++){
                  fp_addd_low(c[i+5], mm[i+6], mm[i+6]);
            }
            dv_copy(m[0],c[8],16);//temp c9
            for(int i=3;i>0;i--){
                  fp_addd_low(c[i+5],c[i+5],c[i+4]);
            }
            fp_addd_low(c[5],c[5],c[10]);
            for(int i=0;i<4;i++){
                  fp_addd_low(m[1], m[i+5], mmm[i]);
                  fp_addd_low(c[i+5], m[1], c[i+5]);
            }

            fp_addd_low(m[1],mm[10],mm[10]);
            fp_addd_low(c[9], m[1], m[0]);
            fp_addd_low(c[9], mmm[4], c[9]);

            fp_addd_low(c[10],mmm[5],mm[0]);
            fp_addd_low(c[10],c[10],m[1]);
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            for(int i = 0; i <6; i++){
                  fp_free(u[i]);
                  fp_free(uu[i]);
            }
            for(int i = 0; i <9; i++)dv_free(m[i]);
            for(int i = 0; i <11; i++){
                  dv_free(mm[i]);	
                  dv_free(mmm[i]);	
            }
	}
}

void fp11_mul_lazyr(fp11_t d, const fp11_t a, const fp11_t b){
      dv11_t c;
      dv11_null(c);
      RLC_TRY{
            dv11_new(c);
            fp11_muln_low(c,a, b);
            for(int i = 0; i<11;i++)fp_rdc(d[i],c[i]);  
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            dv11_free(c);
	}
}

void fp11_mul_nor_low(dv11_t c, dv11_t a) {
      dv11_t b;
      dv11_null(b);
      RLC_TRY{
            dv11_new(b);
            for(int i = 1; i < 11; i++)dv_copy(b[i], a[i-1], 16);
            fp_addd_low(b[0],a[10],a[10]);
            fp_addd_low(b[1],b[1],b[0]);
            for(int i=0;i<11;i++)dv_copy(c[i],b[i],16);
      } RLC_CATCH_ANY {
		RLC_THROW (ERR_CAUGHT);
	} RLC_FINALLY {
            dv11_free(b);
	}
}

void fp11_mul_art(fp11_t c, const fp11_t a) {
      fp11_t b;
      fp11_null(b);
      RLC_TRY{
            fp11_new(b);
            for(int i = 1; i < 11; i++)fp_copy(b[i], a[i-1]);
            fp_dbl(b[0],a[10]);
            fp_add(b[1],b[1],b[0]);
            fp11_copy(c,b);
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            fp11_free(b);
	}     
}

void fp11_mul_invu(fp11_t c, const fp11_t a) {
      fp_t h;
      fp11_t b;
      fp_null(h);
      fp11_null(b);
      RLC_TRY{
            fp11_new(h);
            fp11_new(b);
            for(int i = 1; i < 11; i++)fp_copy(b[i-1], a[i]);
            ep11_curve_get_half(h);
            fp_mul(b[10],a[0],h);
            fp_sub(b[0],b[0],a[0]);
            fp11_copy(c,b);
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            fp_free(h);
            fp11_free(b);
	}     
}