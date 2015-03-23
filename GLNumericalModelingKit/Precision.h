/*
 *  Precision.h
 *  FPSimulator
 *
 *  Created by Jeffrey Early on 3/16/10.
 *  Copyright 2010 Early Innovations, LLC. All rights reserved.
 *
 */

//#import <Accelerate/Accelerate.h>

#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#import <complex.h>
#include <Accelerate/Accelerate.h>
#include <CoreFoundation/CoreFoundation.h>
#import "fftw3.h"
//#include <libMatrixFFT/fftUtils.h>
//#include <libMatrixFFT/MatrixFFT.h>
//#include <libMatrixFFT/vdspUtils.h>
//#include <libMatrixFFT/complexBufUtils.h>

#ifndef _GLPRECISION_
#define _GLPRECISION_

// Must also change this in MatrixFFTConfig.h, line 71.
#define DOUBLE_PREC

#ifdef DOUBLE_PREC

typedef double GLFloat;
typedef double complex GLFloatComplex;
typedef	DSPDoubleSplitComplex	GLSplitComplex;
typedef	FFTSetupD	GLFFTSetup;

#define	vGL_vclr	vDSP_vclrD
#define vGL_vfill   vDSP_vfillD

#define vGL_vmaxmg		vDSP_vmaxmgD			// Vector maximum magnitude
#define	vGL_maxv		vDSP_maxvD
#define	vGL_minv		vDSP_minvD
#define vGL_vthr		vDSP_vthrD
#define	vGL_vclip		vDSP_vclipD
#define	vGL_vthres		vDSP_vthresD
#define vGL_meanv		vDSP_meanvD

#define vGL_vvsin		vvsin
#define vGL_vvcos		vvcos
#define vGL_vvsincos	vvsincos
#define vGL_vvatan      vvatan				// Inverse tanget
#define vGL_vvatan2     vvatan2				// Inverse tanget
#define vGL_vvtanh      vvtanh				// Hyperbolic tanget
#define vGL_vvsinh      vvsinh				// Hyperbolic sine
#define vGL_vvasinh     vvasinh				// Inverse hyperbolic sine
#define vGL_vvcosh      vvcosh				// Hyperbolic cosine
#define vGL_vvacosh     vvacosh				// Inverse hyperbolic cosine
#define vGL_vvexp		vvexp               // Exponential
#define vGL_vvlog       vvlog               // Natural logarithm
#define vGL_vvlog10		vvlog10
#define vGL_vvpow		vvpow
#define vGL_vvsqrt		vvsqrt

#define vGL_vvfloor		vvfloor				// Floor

#define vGL_dotpr		vDSP_dotprD
#define vGL_zdotpr		vDSP_zdotprD

#define vGL_vlint		vDSP_vlintD

#define vGL_vrvrs		vDSP_vrvrsD

#define vGL_zvzsml		vDSP_zvzsmlD		// Complex vector multiply by complex scalar
#define	vGL_vsmul		vDSP_vsmulD			// Vector scalar multipy
#define	vGL_vadd		vDSP_vaddD
#define	vGL_vsmsa		vDSP_vsmsaD
#define	vGL_vsq			vDSP_vsqD			// Square a real vector
#define vGL_zvmags		vDSP_zvmagsD		// Square a complex vector
#define	vGL_sve			vDSP_sveD
#define	vGL_vsadd		vDSP_vsaddD
#define	vGL_vmul		vDSP_vmulD
#define vGL_zvmul		vDSP_zvmulD			// Complex vector multiply by complex vector
#define	vGL_vneg		vDSP_vnegD
#define vGL_vabs		vDSP_vabsD
#define vGL_zvabs		vDSP_zvabsD
#define	vGL_svdiv		vDSP_svdivD
#define vGL_vdiv		vDSP_vdivD			// Vector divide
#define vGL_zvdiv		vDSP_zvdivD			// Complex vector divide
#define	vGL_vma			vDSP_vmaD
#define	vGL_vsma		vDSP_vsmaD
#define	vGL_vfrac		vDSP_vfracD
#define	vGL_vsub		vDSP_vsubD
#define	vGL_vindex		vDSP_vindexD		// Extract values at indices
#define	vGL_vramp		vDSP_vrampD
#define	vGL_vsbsm		vDSP_vsbsmD
#define	vGL_vmmsb		vDSP_vmmsbD
#define	vGL_vmma		vDSP_vmmaD			// Vector multiply, multiply, add  ( A*B + C*D)

#define vGL_mmov		vDSP_mmovD
#define vGL_mtrans		vDSP_mtransD		// Matrix transpose
#define vGL_vsorti		vDSP_vsortiD		// Vector index in-place sort

#define	vGL_fft2d_zip			vDSP_fft2d_zipD
#define	vGL_fft2d_zop			vDSP_fft2d_zopD
#define	vGL_create_fftsetup		vDSP_create_fftsetupD
#define	vGL_destroy_fftsetup	vDSP_destroy_fftsetupD


#define	vGL_ctoz        vDSP_ctozD          // Interleaved Complex->Split Complex
#define	vGL_ztoc        vDSP_ztocD          // Split Complex->Interleaved Complex

#define vGL_vtrapz		vDSP_vtrapzD		// Trapezoidal integration
#define	vGL_vsimps		vDSP_vsimpsD		// Simpson integration

#define vGL_mmul		vDSP_mmulD			// Real matrix multiplication
#define vGL_zmmul		vDSP_zmmulD			// Complex matrix multiplication

#define vGL_ggev		dggev_
#define vGL_geev		dgeev_
#define vGL_gesv		dgesv_
#define vGL_getrf		dgetrf_
#define vGL_getri		dgetri_

#define vGL_fftw_plan							fftw_plan
#define vGL_fftw_plan_with_nthreads				fftw_plan_with_nthreads
#define vGL_fftw_import_wisdom_from_filename	fftw_import_wisdom_from_filename
#define vGL_fftw_export_wisdom_to_filename		fftw_export_wisdom_to_filename
#define vGL_fftw_destroy_plan					fftw_destroy_plan

#define vGL_fftw_iodim							fftw_iodim
#define vGL_fftw_r2r_kind						fftw_r2r_kind

#define vGL_fftw_plan_guru_split_dft			fftw_plan_guru_split_dft
#define vGL_fftw_plan_guru_r2r					fftw_plan_guru_r2r
#define vGL_fftw_plan_guru_split_dft_r2c		fftw_plan_guru_split_dft_r2c
#define vGL_fftw_plan_guru_split_dft_c2r		fftw_plan_guru_split_dft_c2r

#define vGL_fftw_execute						fftw_execute
#define vGL_fftw_execute_split_dft				fftw_execute_split_dft
#define vGL_fftw_execute_r2r					fftw_execute_r2r
#define vGL_fftw_execute_split_dft_r2c			fftw_execute_split_dft_r2c
#define vGL_fftw_execute_split_dft_c2r			fftw_execute_split_dft_c2r

#define vGL_fftw_init_threads					fftw_init_threads

#else

typedef float complex GLFloatComplex;
typedef float GLFloat;
typedef	DSPSplitComplex	GLSplitComplex;
typedef	FFTSetup	GLFFTSetup;

#define	vGL_vclr	vDSP_vclr
#define vGL_vfill   vDSP_vfill

#define vGL_vmaxmg		vDSP_vmaxmg			// Vector maximum magnitude
#define	vGL_maxv		vDSP_maxv
#define	vGL_minv		vDSP_minv
#define vGL_vthr		vDSP_vthr
#define	vGL_vclip		vDSP_vclip
#define	vGL_vthres		vDSP_vthres
#define vGL_meanv		vDSP_meanv

#define vGL_vvsin		vvsinf
#define vGL_vvcos		vvcosf
#define vGL_vvsincos	vvsincosf
#define vGL_vvatan      vvatanf				// Inverse tanget
#define vGL_vvatan2     vvatan2f			// Inverse tanget 2
#define vGL_vvtanh      vvtanhf				// Hyperbolic tanget
#define vGL_vvsinh      vvsinhf				// Hyperbolic sine
#define vGL_vvasinh     vvasinhf			// Inverse hyperbolic sine
#define vGL_vvcosh      vvcoshf				// Hyperbolic cosine
#define vGL_vvacosh     vvacoshf			// Inverse hyperbolic cosine
#define vGL_vvexp		vvexpf              // Exponential
#define vGL_vvlog       vvlogf              // Natural logarithm
#define vGL_vvlog10		vvlog10f
#define vGL_vvpow		vvpowf
#define vGL_vvsqrt		vvsqrtf

#define vGL_vvfloor		vvfloorf			// Floor

#define vGL_dotpr		vDSP_dotpr
#define vGL_zdotpr		vDSP_zdotpr

#define vGL_vlint		vDSP_vlint

#define vGL_vrvrs		vDSP_vrvrs

#define vGL_zvzsml		vDSP_zvzsml			// Complex vector multiply by complex scalar
#define	vGL_vsmul		vDSP_vsmul
#define	vGL_vadd		vDSP_vadd
#define	vGL_vsmsa		vDSP_vsmsa
#define	vGL_vsq			vDSP_vsq			// Square a real vector
#define vGL_zvmags		vDSP_zvmags			// Square a complex vector
#define	vGL_sve			vDSP_sve
#define	vGL_vsadd		vDSP_vsadd
#define	vGL_vmul		vDSP_vmul
#define vGL_zvmul		vDSP_zvmul			// Complex vector multiply by complex vector
#define	vGL_vneg		vDSP_vneg
#define vGL_vabs		vDSP_vabs
#define vGL_zvabs		vDSP_zvabs
#define	vGL_svdiv		vDSP_svdiv
#define vGL_vdiv		vDSP_vdiv			// Vector divide
#define vGL_zvdiv		vDSP_zvdiv			// Complex vector divide
#define	vGL_vma			vDSP_vma
#define	vGL_vsma		vDSP_vsma

#define	vGL_vfrac		vDSP_vfrac
#define	vGL_vsub		vDSP_vsub
#define	vGL_vindex		vDSP_vindex
#define	vGL_vramp		vDSP_vramp
#define	vGL_vsbsm		vDSP_vsbsm
#define	vGL_vmmsb		vDSP_vmmsb
#define	vGL_vmma		vDSP_vmma

#define vGL_mmov		vDSP_mmov
#define vGL_mtrans		vDSP_mtrans			// Matrix transpose
#define vGL_vsorti		vDSP_vsorti			// Vector index in-place sort

#define	vGL_fft2d_zip			vDSP_fft2d_zip
#define	vGL_fft2d_zop			vDSP_fft2d_zop
#define	vGL_create_fftsetup		vDSP_create_fftsetup
#define	vGL_destroy_fftsetup	vDSP_destroy_fftsetup


#define	vGL_ctoz        vDSP_ctoz			// Interleaved Complex->Split Complex
#define	vGL_ztoc        vDSP_ztoc			// Split Complex->Interleaved Complex

#define vGL_vtrapz		vDSP_vtrapz			// Trapezoidal integration
#define	vGL_vsimps		vDSP_vsimps			// Simpson integration

#define vGL_mmul		vDSP_mmul			// Real matrix multiplication
#define vGL_zmmul		vDSP_zmmul			// Complex matrix multiplication

#define vGL_ggev		sggev_
#define vGL_geev		sgeev_
#define vGL_gesv		sgesv_
#define vGL_getrf		sgetrf_
#define vGL_getri		sgetri_

#endif

#endif

