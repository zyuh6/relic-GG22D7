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
 * Implementation of squaring in an undecic extension of a prime field.
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
 * Given a=a[0]+a[1]x and compute a^2
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */

static void fp_sqr_level2(dv_t* c, const fp_t* a){
      fp_addn_low(c[0],a[0],a[1]);
      fp_sqrn_low(c[1],c[0]);
      fp_sqrn_low(c[0],a[0]);
      fp_sqrn_low(c[2],a[1]);
      fp_subd_low(c[1],c[1],c[0]);
      fp_subd_low(c[1],c[1],c[2]);
}
/**
 * Given a=a[0]+a[1]x+a[2]x^2 and compute a^2
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_sqr_level3(dv_t* c, const fp_t* a) {
      dig_t k0[2 * RLC_FP_DIGS], k1[2 * RLC_FP_DIGS];
      for(int i = 0; i<3;i++)fp_sqrn_low(c[2*i],a[i]);

      fp_addn_low(k0,a[0],a[1]);
      fp_sqrn_low(c[1],k0);
      fp_subd_low(c[1],c[1],c[0]);
      fp_subd_low(c[1],c[1],c[2]);

      fp_addn_low(k0,a[1],a[2]);
      fp_sqrn_low(c[3],k0);
      fp_subd_low(c[3],c[3],c[2]);
      fp_subd_low(c[3],c[3],c[4]);

      fp_addn_low(k0,a[0],a[2]);
      fp_sqrn_low(k1,k0);       //
      fp_addd_low(c[2],c[2],k1);
      fp_subd_low(c[2],c[2],c[0]);
      fp_subd_low(c[2],c[2],c[4]);

}
/**
 * Compute a^2, where a=a[0]+a[1]x+...+a[4]x^4
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_sqr_level5(dv_t* c, const fp_t* a) {
      fp_t u[3];
      dv_t m[3],mm[5],mmm[5];  
      for(int i =0;i<3;i++){
            fp_null(u[i]);
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

            for(int i = 0;i<2; i++)fp_addn_low(u[i],a[i],a[i+2]);
            fp_copy(u[2],a[4]);
            fp_sqr_level2(m,a);
            fp_sqr_level3(mm,a+2);
            fp_sqr_level3(mmm,u);
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
            for(int i = 0; i < 3; i++)fp_free(u[i]);
            for(int i = 0; i < 3; i++)dv_free(m[i]);
            for(int i = 0; i < 5; i++){
                  dv_free(mm[i]);	
                  dv_free(mmm[i]);	
            }
      }
}

/**
 * Given a=a[0]+a[1]x+...+a[5]x^5 and compute a^2
 *
 * @param[out] c			- the result.
 * @param[in] a			- the first element.
 * @param[in] b			- the second element.
 */
static void fp_sqr_level6(dv_t* c, const fp_t* a) {

      fp_t u[3];
      dv_t m[5],mm[5],mmm[5];  
      for(int i =0;i<3;i++){
            fp_null(u[i]);
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

            for(int i = 0;i<3; i++)fp_addn_low(u[i],a[i],a[i+3]);
            fp_sqr_level3(m,a);
            fp_sqr_level3(mm,a+3);
            fp_sqr_level3(mmm,u);

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
            for(int i = 0; i < 3; i++)fp_free(u[i]);
            for(int i = 0; i < 5; i++){
                  dv_free(m[i]);
                  dv_free(mm[i]);	
                  dv_free(mmm[i]);	
            }
      }
}



/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/
void fp11_sqrn_low(dv11_t c, const fp11_t a){
      fp_t u[6];
      dv_t m[9],mm[11],mmm[11];	
      for(int i = 0; i <6; i++){
            fp_null(u[i]);
      }
      for(int i = 0; i <9; i++)dv_null(m[i]);
      for(int i = 0; i <11; i++){
            dv_null(mm[i]);	
            dv_null(mmm[i]);	
      }
      RLC_TRY{
            for(int i = 0; i <6; i++){
                  fp_new(u[i]);
            }
            for(int i = 0; i <9; i++)dv_new(m[i]);
            for(int i = 0; i <11; i++){
                  dv_new(mm[i]);	
                  dv_new(mmm[i]);	
            }
            fp_sqr_level5(m,a);
            fp_sqr_level6(mm,a+5);
            for(int i = 0;i<5; i++)fp_addn_low(u[i],a[i],a[i+5]);
            fp_copy(u[5],a[10]);
            fp_sqr_level6(mmm,u);
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
            for(int i = 0; i <6; i++)fp_free(u[i]);
            for(int i = 0; i <9; i++)dv_free(m[i]);
            for(int i = 0; i <11; i++){
                  dv_free(mm[i]);	
                  dv_free(mmm[i]);	
            }
      }
}

void fp11_sqr_lazyr(fp11_t d, const fp11_t a){
      dv11_t c;
      dv11_null(c);
      RLC_TRY{
            dv11_new(c);
            fp11_sqrn_low(c,a);
            for(int i = 0; i<11;i++)fp_rdc(d[i],c[i]);  
      } RLC_CATCH_ANY {
		RLC_THROW(ERR_CAUGHT);
	} RLC_FINALLY {
            dv11_free(c);
	}
}

