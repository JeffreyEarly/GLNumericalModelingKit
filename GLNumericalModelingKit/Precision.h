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
//#include <libMatrixFFT/fftUtils.h>
//#include <libMatrixFFT/MatrixFFT.h>
//#include <libMatrixFFT/vdspUtils.h>
//#include <libMatrixFFT/complexBufUtils.h>

#ifndef _GLPRECISION_
#define _GLPRECISION_

// Must also change this in MatrixFFTConfig.h, line 71.
//#define DOUBLE_PREC

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
#define vGL_vvexp		vvexp               // Exponential
#define vGL_vvlog       vvlog               // Natural logarithm
#define vGL_vvlog10		vvlog10
#define vGL_vvpow		vvpow
#define vGL_vvsqrt		vvsqrt

#define vGL_dotpr		vDSP_dotprD
#define vGL_zdotpr		vDSP_zdotprD

#define vGL_vlint		vDSP_vlintD

#define vGL_vrvrs		vDSP_vrvrsD

#define vGL_zvzsml		vDSP_zvzsmlD		// Complex vector multiply by complex scalar
#define	vGL_vsmul		vDSP_vsmulD
#define	vGL_vadd		vDSP_vaddD
#define	vGL_vsmsa		vDSP_vsmsaD
#define	vGL_vsq			vDSP_vsqD
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
#define	vGL_vindex		vDSP_vindexD
#define	vGL_vramp		vDSP_vrampD
#define	vGL_vsbsm		vDSP_vsbsmD
#define	vGL_vmmsb		vDSP_vmmsbD
#define	vGL_vmma		vDSP_vmmaD

#define vGL_mmov	vDSP_mmovD
#define	vGL_fft2d_zip			vDSP_fft2d_zipD
#define	vGL_fft2d_zop			vDSP_fft2d_zopD
#define	vGL_create_fftsetup		vDSP_create_fftsetupD
#define	vGL_destroy_fftsetup	vDSP_destroy_fftsetupD
#define vGL_zvmags	vDSP_zvmagsD

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
#define vGL_vvatan2     vvatan2f			// Inverse tanget
#define vGL_vvexp		vvexpf              // Exponential
#define vGL_vvlog       vvlogf              // Natural logarithm
#define vGL_vvlog10		vvlog10f
#define vGL_vvpow		vvpowf
#define vGL_vvsqrt		vvsqrtf

#define vGL_dotpr		vDSP_dotpr
#define vGL_zdotpr		vDSP_zdotpr

#define vGL_vlint		vDSP_vlint

#define vGL_vrvrs		vDSP_vrvrs

#define vGL_zvzsml		vDSP_zvzsml			// Complex vector multiply by complex scalar
#define	vGL_vsmul		vDSP_vsmul
#define	vGL_vadd		vDSP_vadd
#define	vGL_vsmsa		vDSP_vsmsa
#define	vGL_vsq			vDSP_vsq
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

#define vGL_mmov	vDSP_mmov
#define	vGL_fft2d_zip			vDSP_fft2d_zip
#define	vGL_fft2d_zop			vDSP_fft2d_zop
#define	vGL_create_fftsetup		vDSP_create_fftsetup
#define	vGL_destroy_fftsetup	vDSP_destroy_fftsetup
#define vGL_zvmags	vDSP_zvmags


#endif

#endif

