c  IMPORTANT SNOWMODEL CODE (TRANSLATION OF MELT ENERGY INTO SNOW MELT FOR SINGLE LAYER SNOWPACK OPTION)



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE NSNOWDEN(ro_nsnow,Twb,Tf,dt)

      implicit none

      real Twgmax,Tf,Twb,ro_nsnow,scalefact,dt,wt

      Twgmax = Tf + 1.0
      if (Twb.ge.258.16 .and. Twb.le.Twgmax) then
        ro_nsnow = 50. + 1.7 * (Twb - 258.16)**1.5
      elseif (Twb.lt.258.16) then
        ro_nsnow = 50.0
      else
        ro_nsnow = 158.8
      endif

c For one day time steps, this equation gives a new snow density at
c   the end of the 24 hour period which is too low, by an approximate
c   factor of X.  Thus, for a daily time step, I scale the density by
c   X before returning it to the main program.

      scalefact = 1.0
      if (dt.eq.86400.0) then
        if (ro_nsnow.le.158.8) then
          wt = 1.0 + (50.0 - ro_nsnow) / 108.8
          ro_nsnow = wt * scalefact * ro_nsnow + ro_nsnow
          ro_nsnow = min(158.8,ro_nsnow)
        endif
      endif

      return
      end
	  
	  
	  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE DDENSITY(ro_snow_grid,swe_depth,Tf,Tsfc,dt,A1,A2,
     &  snow_depth,ro_water,ro_ice)

      implicit none

      real snow_depth,Tsg,Tf,Tsnow,Tsfc,ro_snow_grid,dt,A1,A2,
     &  swe_depth_star,ro_ice,ro_water,swe_depth

      if (snow_depth.gt.0.0) then

c Assume that the snow-ground interface temperature is -1.0 C.
        Tsg = Tf - 1.0
        Tsnow = 0.5 * (Tsg + Tsfc)
        swe_depth_star= 0.5 * swe_depth
        ro_snow_grid = ro_snow_grid + dt *
     &    (A1 * swe_depth_star * ro_snow_grid *
     &    exp((- 0.08)*(Tf-Tsnow)) * exp((- A2)*ro_snow_grid))
        ro_snow_grid = min(ro_ice,ro_snow_grid)
        snow_depth = ro_water * swe_depth / ro_snow_grid

      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SNOWPACK(swe_depth,snow_d,ro_snow_grid,
     &  prec,ro_water,ro_nsnow,runoff,Qm,xLf,dt,rain,sprec,
     &  sum_prec,sum_runoff,soft_snow_d,sum_sprec,ro_snow,
     &  snow_depth,sprec_grnd,vegtype,glacier_melt,sum_glacmelt,
     &  swemelt,canopy_unload,Qe,sfc_sublim_flag,sum_sfcsublim,
     &  sum_swemelt,corr_factor,icorr_factor_index,swesublim,
     &  ro_snowmax)

      implicit none

      real ro_snowmax,runoff,Qm,swe_depth,potmelt,swemelt,dt,
     &  ro_water,xLf,snow_depth,ro_snow_grid,snow_d_melt,dz_water,
     &  soft_snow_d,prec,rain,snow_d,sum_sprec,sum_prec,
     &  sum_runoff,ro_nsnow,sprec,ro_snow,snow_d_new,sprec_grnd,
     &  vegtype,glacier_melt,sum_glacmelt,canopy_unload,Qe,
     &  xLsublim,potsublim,swesublim,snow_d_sublim,sfc_sublim_flag,
     &  sum_sfcsublim,sum_swemelt,corr_factor,potmelt_tmp
      integer icorr_factor_index

      runoff = 0.0

c SURFACE SUBLIMATION.

c Whether static-surface (non-blowing snow) sublimation is included
c   in the model calculations is controlled by the sfc_sublim_flag.
c   I am waiting for the flux-tower data Matthew and I are collecting
c   in Alaska, to compare with the model simulations, before
c   including this part of the model in all simulations.

c If the sfc_sublim_flag is turned on, the latent heat flux (Qe)
c   calculated in ENBAL is used to add/remove snow from the snowpack.
c   xLsublim = xLf + xLv = 2.5104x10^6 J/kg + 3.334x10^5 J/kg, and
c   potsublim is in m swe.

      if (swe_depth.gt.0.0  .and.  sfc_sublim_flag.eq.1.0) then
        if (Qe.lt.0.0) then

c Compute the snow-surface sublimation (m, swe).
          xLsublim = 2.844e6
          potsublim = (- dt) * Qe / (ro_water * xLsublim)
          swesublim = min(potsublim,swe_depth)

c Save a summing array of the static surface snow sublimation.
          sum_sfcsublim = sum_sfcsublim + swesublim

c Compute the change in snow depth.  Assume that this sublimated
c   snow does not change the snow density and does not change the
c   soft snow depth.  It only reduces the snow depth and the
c   associated swe depth.
          swe_depth = swe_depth - swesublim
          if (swe_depth.eq.0.0) then
            snow_depth = 0.0
          else
            snow_d_sublim = swesublim * ro_water / ro_snow_grid
            snow_depth = snow_depth - snow_d_sublim
          endif
        else
          swesublim = 0.0
        endif
      else
        swesublim = 0.0
      endif

c MELTING.

c If melting occurs, decrease the snow depth, and place the melt
c   water in the 'runoff' variable.  Keep track of the liquid water
c   produced.

      if (Qm.gt.0.0) then

c Compute the snow melt (m).
        potmelt = dt * Qm / (ro_water * xLf)

c Account for any snowmelt data assimilation.
        if (icorr_factor_index.lt.0.0) then
          potmelt_tmp = potmelt * corr_factor
          swemelt = min(potmelt_tmp,swe_depth)
c Handle the case of no snowmelt data assimilation.
        else
          swemelt = min(potmelt,swe_depth)
        endif

c Compute any glacier or permanent snow-field melt (m water equiv.).
        if (vegtype.eq.20.0) then
          glacier_melt = potmelt - swemelt
        else
          glacier_melt = 0.0
        endif

c Save a summing array of the glacier melt.
        sum_glacmelt = sum_glacmelt + glacier_melt

c Save the runoff contribution.
        runoff = runoff + glacier_melt

c Save a summing array of the snow melt.
        sum_swemelt = sum_swemelt + swemelt

c Compute the change in snow depth.
        snow_d_melt = swemelt * ro_water / ro_snow_grid
        snow_depth = snow_depth - snow_d_melt
        snow_depth = max(0.0,snow_depth)

c Compute the changes in snow density resulting from the melt.
c   Assume that the melted snow is redistributed through the new
c   snow depth up to a maximum density.  Any additional melt water
c   is added to the runoff.
        if (snow_depth.eq.0.0) then
          ro_snow_grid = ro_snowmax
          runoff = runoff + swemelt
        else
          ro_snow_grid = swe_depth * ro_water / snow_depth
        endif

        if (ro_snow_grid.gt.ro_snowmax) then
          dz_water = snow_depth *
     &      (ro_snow_grid - ro_snowmax) / ro_water
          ro_snow_grid = ro_snowmax
          swe_depth = snow_depth * ro_snow_grid / ro_water
          runoff = runoff + dz_water
        else
          swe_depth = snow_depth * ro_snow_grid / ro_water
        endif

        soft_snow_d = 0.0

      else

c These prevent values from the previous time step from being
c   carried through to the next time step.
        swemelt = 0.0
        glacier_melt = 0.0

      endif

c PRECIPITATION.

c Precipitation falling as rain on snow contributes to a snow
c   density increase, precipitation falling as snow adds to the
c   snow depth, and rain falling on bare ground contributes to the
c   runoff.

c We have precipitation.
      if (prec.gt.0.0) then
        rain = prec - sprec

c We have rain.
        if (rain.gt.0.0) then

c Rain on snow.  Note that we can also have snow unloading here.
c   Assume this unloading is wet as rain.
          if (snow_depth.gt.0.0) then
            swe_depth = swe_depth + rain + canopy_unload
            ro_snow_grid = swe_depth * ro_water / snow_depth
            if (ro_snow_grid.gt.ro_snowmax) then
              dz_water = snow_depth * (ro_snow_grid - ro_snowmax) /
     &          ro_water
              ro_snow_grid = ro_snowmax
              swe_depth = snow_depth * ro_snow_grid / ro_water
              runoff = runoff + dz_water
            endif

c Rain on bare ground.  Assume any unloading is as wet as rain.
          else
            runoff = runoff + rain + canopy_unload
          endif

c We have snow precipitation (on either snow or bare ground).
        else
          swe_depth = swe_depth + sprec_grnd
          snow_d_new = ro_water / ro_nsnow * sprec_grnd
          snow_depth = snow_depth + snow_d_new
          ro_snow_grid = ro_water * swe_depth / snow_depth
        endif

c Here we handle the case where there is no precipitation, but
c   there is snow falling from the canopy to the snowpack.
      else
        rain = 0.0
        if (sprec_grnd.gt.0.0) then
          swe_depth = swe_depth + sprec_grnd
          snow_d_new = ro_water / ro_snow * sprec_grnd
          snow_depth = snow_depth + snow_d_new
          ro_snow_grid = ro_water * swe_depth / snow_depth
        endif
      endif

c The following are set up to be compatible with SnowTran-3D, and
c   are in snow-depth units.  The sum_sprec corrections are done
c   in the SnowTran-3D code.
      soft_snow_d = soft_snow_d + sprec_grnd * ro_water / ro_snow
      snow_d = swe_depth * ro_water / ro_snow
c     sum_sprec = sum_sprec + sprec_grnd * ro_water / ro_snow
      sum_sprec = sum_sprec + sprec_grnd

c The following are in swe-depth units.
      sum_prec = sum_prec + prec
      sum_runoff = sum_runoff + runoff

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc