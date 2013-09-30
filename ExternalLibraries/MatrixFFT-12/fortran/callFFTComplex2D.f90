!	File: callFFTComplex2D.f90
!	
!	Description:
!		Example to demonstrate performing a 2-dimension complex FFT using 
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

	program callC
	
	implicit none
	
	double precision realPart(8, 4) 
	double precision imagPart(8, 4) 
	
	double precision oneDReal(32)
	double precision oneDImag(32)
	
	integer i, j
	integer count
	
	integer size(2)  ! array of sizes, length 2 for 2D fft
	integer fwd
	integer isReal
	integer dims ! number of dimensions. 1 and 2 are supported. 
	integer res
    
	size(1) = 8
	size(2) = 4
	isReal = 0
	dims = 2
    
	count = 0 
	
	do i=1,size(1)
		do j =1,size(2)
			realPart(i,j) = dble( count )
			imagPart(i,j) = dble( count + 1)
			count = count + 2
		enddo
	enddo

		
	print *, "original real"
	call printMatrix(realPart, size(1), size(2) )
	print *, "original imag"
	call printMatrix(imagPart, size(1), size(2) )
	
	! copy input arrays to C style 1D arrays for input
	call copyTo1D(realPart, oneDReal, size(1), size(2))
	call copyTo1D(imagPart, oneDImag, size(1), size(2))					
			
    ! Forward FFT
	fwd = 1
	call fft_wrapper(oneDReal,oneDImag,fwd, isReal, size, dims, res )
	
	! copy back for output 
	! data to pass back to library is not modified and remains in 1D arrays
	call returnTo2D(realPart, oneDReal, size(1), size(2))
	call returnTo2D(imagPart, oneDImag, size(1), size(2))	
	
	
	print *, "forward FFT result", res
	print *, "forward FFT real"
	call printMatrix(realPart, size(1), size(2) )
	print *, "forward FFT imag"
	call printMatrix(imagPart, size(1), size(2) )
			
    ! Inverse FFT
	fwd = 0
	call fft_wrapper(oneDReal,oneDImag,fwd, isReal, size, dims, res )
	
	! copy back for output 
	! output data not modified
	call returnTo2D(realPart, oneDReal, size(1), size(2))
	call returnTo2D(imagPart, oneDImag, size(1), size(2))
	
	print *, "inverse FFT result", res
	print *, "inverse FFT real"
	call printMatrix(realPart, size(1), size(2) )
	print *, "inverse FFT imag"
	call printMatrix(imagPart, size(1), size(2) )
        
    ! Release resources allocated in fft_wrapper
	call fft_release( )
    
	return 
	end
	
	
	

! helper routines
	
	subroutine printMatrix( x, m, n )
		
		!outputs a double precision matrix of size m by n
		
		integer m, n
		integer i, j
		
		double precision x(m,n)
		
		do i = 1,m
			print *,(x(i,j), j=1,n)
		enddo
	
	return
	end
	

	
	subroutine copyTo1D( twoDim, oneDim, m, n)
	! copies a 2D array in column major order to a 1D array in row major order
	! puts standard 2D Fortran format into 1D format expected by C library

	! calling this routine can be avoided by careful placement in 1D array. 
	! this is advised for performance, as the current scheme wastes memory

		integer m, n
		integer i, j
		integer count
	
		double precision twoDim(m,n)
		double precision oneDim(m*n)
		
		count = 1 
	
		do i = 1,m
			do j = 1,n
				oneDim(count) = twoDim(i,j)
				count = count + 1
			enddo
		enddo
	
	return 
	end
	


	subroutine returnTo2D( twoDim, oneDim, m, n)
	! reverses the above copy
	! places a row-major one dimensional array in a column-major, Fortran style, 2D array

		integer m, n
		integer i, j
		integer count
	
		double precision twoDim(m,n)
		double precision oneDim(m*n)
		
		count = 1 
	
		do i = 1,m
			do j = 1,n
				twoDim(i,j) = oneDim(count)
				count = count + 1
			enddo
		enddo
	
	return 
	end
	
	