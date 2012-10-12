pro build_wisesfrs, clobber=clobber, nmonte=nmonte, debug=debug
; jm11apr28ucsd - take the output of BUILD_WISESFRS_PARENT and compute
; everything we need to write a paper
    
    sfrspath = getenv('WISENYU_DATA')+'/wisesfrs/'

    kcorr = mrdfits(sfrspath+'wisesfrs_kcorr_all.fits.gz',1)
    wise = mrdfits(sfrspath+'wisesfrs_wise.fits.gz',1)
    mass = mrdfits(sfrspath+'wisesfrs_mpamass.fits.gz',1)
    ispec = mrdfits(sfrspath+'wisesfrs_mpaispec.fits.gz',1)
    witt = mrdfits(sfrspath+'wisesfrs_witt.fits.gz',1)
    ngal = n_elements(wise)

    if (n_elements(nmonte) eq 0) then nmonte = 0
    
; initialize the output data structure
    out = {$
      isagn:            0,$
      class:           '',$
      ewha:        -999.0,$
      ewha_err:    -999.0,$
      d4000:       -999.0,$
      d4000_err:   -999.0,$
      mpamass:     -999.0,$
      mpamass_err: -999.0,$
      mpaoh:       -999.0,$
      mpaoh_err:   -999.0,$
      mpasfr:      -999.0,$
      mpasfr_err:  -999.0,$
                 
      k_mass:              -999.0,$ ; from K-correct
      k_coeffs:         fltarr(5),$
      k_maggies:       fltarr(12),$
      k_ivarmaggies:   fltarr(12),$
      wise_maggies:     fltarr(4),$ ; aperture 8 fluxes
      wise_ivarmaggies: fltarr(4),$
      wise_maggies_mpro:     fltarr(4),$ ; profile fluxes
      wise_ivarmaggies_mpro: fltarr(4),$

      Mr:        -999.0,$       ; absolute magnitude
      umr:       -999.0,$ ; u-r (rest)
      rmj:       -999.0,$ ; r-J (rest)
      nuvmr:     -999.0,$ ; NUV-r (rest)

      fuvmnuv:     -999.0,$     ; FUV-NUV color (rest)
      fuvmnuv_err: -999.0,$
      l1500:     -999.0,$
      l1500_err: -999.0,$
      beta:      -999.0,$
      beta_err:  -999.0,$

      ch1mch2:   -999.0,$ ; [3.4]-[4.6]
      ch1mch3:   -999.0,$ ; [3.4]-[12]
      ch2mch3:   -999.0,$ ; [4.6]-[12]

      m22good:        0,$ ; significant 22-micron detection
      m22upper:       0,$ ; upper limit on the 22-micron flux
      lir:      -999.0,$ ; average!
      lir_err:  -999.0,$
      irx:      -999.0,$
      irx_err:  -999.0,$
      modelindx_ce01: -999.0,$
      modelindx_dh02: -999.0,$

      a1500:     -999.0,$
      l1500_cor: -999.0,$

      sfr1500_obs: -999.0,$ ; from L(1500)_obs
      sfr1500_cor: -999.0,$ ; from L(1500)_cor
      sfrtot:      -999.0}  ; from UV+IR
    out = replicate(out,ngal)    

    out = struct_addtags(struct_trimtags(kcorr,$
      select=['object_position','ra','dec','z']),out)

; get the spectroscopic properties
    out.ewha = ispec.h_alpha_ew[0]
    out.ewha_err = ispec.h_alpha_ew[1]
    out[where(out.ewha gt 1E4)].ewha = -999.0 ; crap
    out[where(out.ewha gt 1E4)].ewha_err = -999.0 
    out.d4000 = ispec.d4000_narrow_cor[0]
    out.d4000_err = ispec.d4000_narrow_cor[1]

    out.mpaoh = mass.oh_avg
    out.mpamass = mass.mass_avg
    out.mpasfr = mass.sfr_avg
    out.mpaoh_err = (mass.oh_p84-mass.oh_p16)/2.0
    out.mpamass_err = (mass.mass_p84-mass.mass_p16)/2.0
    out.mpasfr_err = (mass.sfr_p84-mass.sfr_p16)/2.0

    cc = iclassification(ispec,snrcut=2.0,ratios=ratios,doplot=doplot)
    out.class = strtrim(ratios.final_class,2)
    out.isagn = strmatch(out.class,'*AGN*')

; get the optical properties
    out.k_mass = kcorr.k_mass
    out.k_coeffs = kcorr.k_coeffs
    out.k_maggies = kcorr.k_maggies
    out.k_ivarmaggies = kcorr.k_ivarmaggies
    
    out.Mr = kcorr.k_ugrizjhk_absmag[2]
    out.umr = kcorr.k_ugrizjhk_absmag[0]-kcorr.k_ugrizjhk_absmag[2]
    out.rmj = kcorr.k_ugrizjhk_absmag[2]-kcorr.k_ugrizjhk_absmag[5]
    out.nuvmr = kcorr.k_galex_absmag[1]-kcorr.k_ugrizjhk_absmag[2]

; get the UV properties; to get the error in the UV slope Monte Carlo
; the K-corrections
    nuvweff = k_lambda_eff(filterlist='galex_NUV.par',band_shift=0.1)
    fuvweff = k_lambda_eff(filterlist='galex_FUV.par',band_shift=0.1)
    mbolsun = 4.74

    l1500factor = 1500D*4.0*!dpi*dluminosity(kcorr.z,/cm)^2/3.826D33
    out.l1500 = l1500factor*kcorr.k_uvflux[0]
    out.fuvmnuv = kcorr.k_galex_absmag[0]-kcorr.k_galex_absmag[1]

    if (nmonte gt 0) then begin
       t0 = systime(1)
       splog, 'Monte Carlo error on FUV-NUV and L(1500)', nmonte
       filterlist = [galex_filterlist(),sdss_filterlist(),$
         twomass_filterlist(),(wise_filterlist())[0:1]]
       nfilt = n_elements(filterlist)
       err = 1.0/sqrt(kcorr.k_ivarmaggies+(kcorr.k_ivarmaggies eq 0))*(kcorr.k_ivarmaggies ne 0)
       fuvmnuv_monte = fltarr(ngal,nmonte)
       l1500_monte = fltarr(ngal,nmonte)
       for ii = 0L, nmonte-1 do begin
          maggies = kcorr.k_maggies + randomn(seed,nfilt,ngal)*err
          kk = wisesfrs_do_kcorrect(out.z,maggies,kcorr.k_ivarmaggies,$
            filterlist=filterlist,/silent)
          fuvmnuv_monte[*,ii] = kk.k_galex_absmag[0]-kk.k_galex_absmag[1]
          l1500_monte[*,ii] = l1500factor*kk.k_uvflux[0]
       endfor
       for jj = 0L, ngal-1 do out[jj].fuvmnuv_err = djsig(fuvmnuv_monte[jj,*])
       for jj = 0L, ngal-1 do out[jj].l1500_err = djsig(l1500_monte[jj,*])/out.l1500/alog(10)
       splog, 'Total time = ', (systime(1)-t0)/60
    endif else begin
       out.fuvmnuv_err = 0.1
       out.l1500_err = 0.1 ; [dex]
    endelse
    out.l1500 = alog10(out.l1500) ; take the log

    out.beta = 2.32*out.fuvmnuv-2 ; from Hao+11, via Kong+04
    out.beta_err = 2.32*out.fuvmnuv_err

; take the uncertainty in L(1500) to be that of FUV_rest, which is at
; 1390 A, so not a bad approximation
;   out.fuvmnuv_err = 1.0/sqrt(kcorr.k_galex_absmag_ivar[0])
;   out.fuvmnuv_err = sqrt((1.0/sqrt(kcorr.k_galex_absmag_ivar[0]))^2+$
;     (1.0/sqrt(kcorr.k_galex_absmag_ivar[1]))^2)
;   out.l1500 = 10^(-0.4*(kcorr.k_galex_absmag[0]-4.74))
;   out.l1500_err = 0.4*alog(10)/sqrt(kcorr.k_galex_absmag_ivar[0])

; finally get the IR properties
    splog, 'Getting IR properties'
    wise_to_maggies, wise, maggies, ivarmaggies
    wise_to_maggies, wise, maggies_mpro, ivarmaggies_mpro, /mpro
    out.wise_maggies = maggies
    out.wise_ivarmaggies = ivarmaggies
    out.wise_maggies_mpro = maggies_mpro
    out.wise_ivarmaggies_mpro = ivarmaggies_mpro
    
; ch1-ch2 and ch1-ch3 colors    
    good12 = where((maggies[0,*] gt 0.0) and (maggies[1,*] gt 0.0))
    good23 = where((maggies[1,*] gt 0.0) and (maggies[2,*] gt 0.0))
    good13 = where((maggies[0,*] gt 0.0) and (maggies[2,*] gt 0.0))
    out[good12].ch1mch2 = reform(-2.5*alog10(maggies[0,good12]/maggies[1,good12]))
    out[good23].ch2mch3 = reform(-2.5*alog10(maggies[1,good23]/maggies[2,good23]))
    out[good13].ch1mch3 = reform(-2.5*alog10(maggies[0,good13]/maggies[2,good13]))

; get L(IR), dealing with upper limits    
    gd22 = where((maggies[3,*] gt 0.0),comp=up22)
    fitmaggies = maggies
    fitivarmaggies = ivarmaggies
    fitmaggies[3,up22] = 1.0/sqrt(ivarmaggies[3,up22])/2.0 ; 1-sigma
    fitivarmaggies[3,up22] = median(ivarmaggies[3,gd22])    ; typical value

    m22lim = 14.45
    out.m22good = reform((fitmaggies[3,*] ge 10^(-0.4*m22lim))) ; m(22)<14.45
    out.m22upper = out.m22good and reform(maggies[3,*] eq 0)    ; m(22) upper limits

    m22good = where(out.m22good,ngood)
    if (ngood ne 0L) then begin
       lir_ce01 = wise2lir(out[m22good].z,fitmaggies[*,m22good],$
         fitivarmaggies[*,m22good],/chary,/no12micron,debug=debug,$
         model_indx=indx_ce01)
       lir_dh02 = wise2lir(out[m22good].z,fitmaggies[*,m22good],$
         fitivarmaggies[*,m22good],/dale,/no12micron,debug=debug,$
         model_indx=indx_dh02)

; get the statistical errors on L(IR) using Monte Carlo
       lir = (lir_ce01+lir_dh02)/2.0
       lir_err = abs(lir_ce01-lir_dh02)
       out[m22good].lir = alog10(lir)
       out[m22good].lir_err = lir_err/lir/alog(10)
       out[m22good].modelindx_ce01 = indx_ce01[m22good]
       out[m22good].modelindx_dh02 = indx_dh02[m22good]
       
; compute the infrared excess and A(1500) from the Witt & Gordon dust
; models; also correct L(1500) for attenuation
       out[m22good].irx = out[m22good].lir-out[m22good].l1500
       out[m22good].irx_err = out[m22good].lir_err

       linterp, witt.irx, witt.a1500, 10^out[m22good].irx, a1500, missing=0.0
       out[m22good].a1500 = a1500
       out[m22good].l1500_cor = out[m22good].l1500 + 0.4*a1500

; calculate the SFR
       monofactor = 1500.0/im_light(/ang)
       l1500_obs = 10^out[m22good].l1500*3.826D33  ; [erg/s]
       l1500_cor = 10^out[m22good].l1500_cor*3.826D33 ; [erg/s]
       lir = 10^out[m22good].lir*3.826D33             ; [erg/s]
       
       out[m22good].sfr1500_obs = alog10(1.4D-28*monofactor*l1500_obs)
       out[m22good].sfr1500_cor = alog10(1.4D-28*monofactor*l1500_cor)
       out[m22good].sfrtot = alog10(9.8D-11*(lir+2.2*l1500_obs)/3.826D33)
    endif
       
; write out    
    outfile = sfrspath+'wisesfrs.fits'
    im_mwrfits, out, outfile, clobber=clobber
    
stop    

return
end
    
