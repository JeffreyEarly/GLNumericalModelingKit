!	File: callFFTReal.f90
!	
!	Description:
!		Example to demonstrate performing a 1-dimension real-signal FFT using 
!       MatrixFFT in Fortran.
!	
!	Copyright:
!		Copyright (C) 2009 Apple Inc.  All rights reserved.
!	
!	Disclaimer:
!		IMPORTANT:  This Apple software is supplied to you by Apple
!		Computer, Inc. ("Apple") in consideration of your agreement to
!		the following terms, and your use, installation, modification
!		or redistribution of this Apple software constitutes acceptance
!		of these terms.  If you do not agree with these terms, please
!		do not use, install, modify or redistribute this Apple
!		software.
!
!		In consideration of your agreement to abide by the following
!		terms, and subject to these terms, Apple grants you a personal,
!		non-exclusive license, under Apple’s copyrights in this
!		original Apple software (the "Apple Software"), to use,
!		reproduce, modify and redistribute the Apple Software, with or
!		without modifications, in source and/or binary forms; provided
!		that if you redistribute the Apple Software in its entirety and
!		without modifications, you must retain this notice and the
!		following text and disclaimers in all such redistributions of
!		the Apple Software.  Neither the name, trademarks, service
!		marks or logos of Apple Computer, Inc. may be used to endorse
!		or promote products derived from the Apple Software without
!		specific prior written permission from Apple.  Except as
!		expressly stated in this notice, no other rights or licenses,
!		express or implied, are granted by Apple herein, including but
!		not limited to any patent rights that may be infringed by your
!		derivative works or by other works in which the Apple Software
!		may be incorporated.
!
!		The Apple Software is provided by Apple on an "AS IS" basis.
!		APPLE MAKES NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
!		WITHOUT LIMITATION THE IMPLIED WARRANTIES OF NON-INFRINGEMENT,
!		MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, REGARDING
!		THE APPLE SOFTWARE OR ITS USE AND OPERATION ALONE OR IN
!		COMBINATION WITH YOUR PRODUCTS.
!
!		IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT,
!		INCIDENTAL OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
!		TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!		DATA, OR PROFITS; OR BUSINESS INTERRUPTION) ARISING IN ANY WAY
!		OUT OF THE USE, REPRODUCTION, MODIFICATION AND/OR DISTRIBUTION
!		OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER UNDER THEORY
!		OF CONTRACT, TORT (INCLUDING NEGLIGENCE), STRICT LIABILITY OR
!		OTHERWISE, EVEN IF APPLE HAS BEEN ADVISED OF THE POSSIBILITY OF
!		SUCH DAMAGE.
!

! Please see section "Data Packing for Real FFTs" in the following for instructions
! http://developer.apple.com/hardwaredrivers/ve/downloads/vDSP_Library.pdf

	program callC
	
	implicit none
	
	double precision realPart(8) 
	double precision imagPart(8) 
	
	integer val 
	integer isReal
	integer size
	integer i
	integer fwd
	integer dims
	integer res
    
	val = 1
	size = 16
	isReal = 1  ! call as real input data 
	dims = 1
	
    ! MatrixFFT expects real data as inteleaved between 2 parallel arrays.
    ! Initialize to incrementing real data.
	do i=1,size/2
		realPart(i) = dble( val - 1)
		imagPart(i) = dble( val )
		val = val + 2
	enddo
		
	print *, "original real", realPart
	print *, "original imag", imagPart 	
		
    ! Forward FFT
	fwd = 1
	call fft_wrapper(realPart,imagPart, fwd, isReal, size, dims, res )
	
	print *, "forward FFT result", res
	print *, "forward FFT real", realPart
	print *, "forward FFT imag", imagPart
		
    ! Inverse FFT
	fwd = 0
	call fft_wrapper(realPart,imagPart, fwd, isReal, size, dims, res )
	
	print *, "inverse FFT result", res
	print *, "inverse FFT real", realPart
	print *, "inverse FFT imag", imagPart
    
    ! Release resources allocated in fft_wrapper
	call fft_release( )
	return 
	end
	