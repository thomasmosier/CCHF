cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,
     &  sfc_pressure,ht_windobs,windspd,ro_air,Cp,emiss_sfc,
     &  Stef_Boltz,gravity,xLs,xkappa,z_0,Tf,Qc)

      implicit none

      real Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,sfc_pressure,ht_windobs,
     &  windspd,ro_air,Cp,emiss_sfc,Stef_Boltz,gravity,xLs,xkappa,
     &  z_0,Tf,Qc,AAA,CCC,DDD,EEE,FFF,C1,C2,B1,B2,z_0_tmp

      AAA = ro_air * Cp * De_h
      CCC = 0.622 / sfc_pressure
      DDD = emiss_sfc * Stef_Boltz
      EEE = (1.0-albedo) * Qsi + Qli + Qc
      FFF = ro_air * xLs * De_h

c Compute the constants used in the stability coefficient
c   computations.
      z_0_tmp = min(0.25*ht_windobs,z_0)
      C1 = 5.3 * 9.4 * (xkappa/(log(ht_windobs/z_0_tmp)))**2 *
     &  sqrt(ht_windobs/z_0_tmp)
      C2 = gravity * ht_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      CALL SOLVE(Tsfc,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLVE(xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf)

      implicit none

      integer maxiter,i

      real tol,old,A,B,C,other1,other2,es0,dother1,dother2,xnew,
     &  Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,B3,stability,
     &  dstability,fprime1,fprime2,fprime3,fprime4,B8,funct,
     &  fprime

      tol = 1.0e-2
      maxiter = 20
      old = Tair

c Because I am interested in sublimation over snow during the
c   winter, do these calculations over ice.

c Over water.
c       A = 6.1121 * 100.0
c       B = 17.502
c       C = 240.97
c Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55

      do i=1,maxiter

c This section accounts for an increase in turbulent fluxes
c   under unstable conditions.
        other1 = AAA * (Tair - old)
        es0 = A * exp((B * (old - Tf))/(C + (old - Tf)))
        other2 = FFF*CCC*(ea-es0)

        dother1 = - AAA
        dother2 = (- FFF)*CCC*es0*B*C/((C + (old - Tf))**2)

      if (old.gt.Tair) then
c Unstable case.
        B3 = 1.0 + B2 * sqrt(old - Tair)
        stability = 1.0 + B1 * (old - Tair) / B3
        dstability = B1/B3 - (B1*B2*(old-Tair))/
     &    (2.0*B3*B3*sqrt(old-Tair))
        fprime1 = (- 4.0)*DDD*old**3
        fprime2 = stability * dother1 + other1 * dstability
        fprime3 = stability * dother2 + other2 * dstability
        fprime4 = - 0.0

      elseif (old.lt.Tair) then
c Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - old))**2)
        dstability = 2.0 * B8 / ((1.0 + B8 * (Tair - old))**3)
        fprime1 = (- 4.0)*DDD*old**3
        fprime2 = stability * dother1 + other1 * dstability
        fprime3 = stability * dother2 + other2 * dstability
        fprime4 = - 0.0

      else
c Neutrally stable case.
        stability = 1.0
        fprime1 = (- 4.0)*DDD*old**3
        fprime2 = dother1
        fprime3 = dother2
        fprime4 = - 0.0
      endif

        funct = EEE - DDD*old**4 + AAA*(Tair-old)*stability +
     &    FFF*CCC*(ea-es0)*stability +
     &    0.0
        fprime = fprime1 + fprime2 + fprime3 + fprime4

        xnew = old - funct/fprime

        if (abs(xnew - old).lt.tol) return
        old = xnew

      end do

c If the maximum iterations are exceeded, send a message and set
c   the surface temperature to the air temperature.
      write (*,102)
  102 format('max iteration exceeded when solving for Tsfc')
      xnew = Tair

      return
      end