!Calculates various error metrics to compare results in cxgood/cxfin

Program gearproc
	Integer, parameter :: sp = Kind(0.0E0)
	Integer, parameter :: dp = Kind(0.0D0)
	Integer, parameter :: qp = Kind(0.0Q0)
	Integer, parameter :: r_kind = qp
	Integer, parameter :: Spcnum = 121
	Integer, parameter :: Cellnum = 6912
	Real, parameter :: Eps = 1.0e-16
	Real (r_kind), dimension(Cellnum,Spcnum) :: Cxf,Cxg,Logdiff !Test values Cxf, "good" values Cxg
	Real (r_kind) :: Mean, Std2
print*, "after declars"
	Cxf = 0
	Cxg = 0

	open (unit = 2, file = "cxfin")
	open (unit = 3, file = "cxgood")
	print*, "opened files"
	!Remove header trash, cx(itloop,IGAS) = x by y
	read(2,*)
	read(3,*)

	do i=1,Cellnum
		!Remove cell info trash, Cell: n Max/Min = ...
		read(2,*)
		read(3,*)
		do j=1,Spcnum
			read(2,*) Cxf(i,j)
			read(3,*) Cxg(i,j)
			!K: If less than double-eps, ignore in calculation
			if (Cxf(i,j) .lt. Eps) then
				Cxf(i,j) = Cxg(i,j)
			endif
		enddo
	enddo

	Logdiff = abs(Log10(Cxf) - Log10(Cxg))

	!Print and remove outliers
	do i=1,Cellnum
		do j=1,Spcnum
			if (Logdiff(i,j) .gt. 1) then
				Write(*,*) 'Cell/Spc = ', i, '/',j
				Write(*,*) '  Cxg = ', Cxg(i,j)
				Write(*,*) '  Cxf = ', Cxf(i,j)
				Cxf(i,j) = Cxg(i,j)
			endif
		enddo
	enddo
	Logdiff = abs(Log10(Cxf)-Log10(Cxg))

	Mean = Sum(Logdiff)/(Spcnum*Cellnum)
   Std2 = Sum( (Logdiff-Mean)**2.0 )/(Spcnum*Cellnum)
	do i=1,Cellnum
		do j=1,Spcnum
			if (Logdiff(i,j) .gt. 0.25) then
				Write(*,*) 'Cell/Spc = ', i, '/', j
				Write(*,*) '  Cxg = ', Cxg(i,j)
				Write(*,*) '  Cxf = ', Cxf(i,j)
			end if
		end do
	end do

	Write(*,*) 'X = |Log(C_ij) - Log(D_ij)|'
	Write(*,*) '  Max  = ', Maxval(Logdiff)
	Write(*,*) '  Mean = ', Mean
	Write(*,*) '  Std  = ', sqrt(Std2) 
	close(2)
	close(3)
end Program gearproc
