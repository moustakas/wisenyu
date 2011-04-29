pro build_wisesfrs, clobber=clobber, debug=debug
; jm11apr28ucsd - take the output of BUILD_WISESFRS_PARENT and compute
; everything we need to write a paper
    
    sfrspath = getenv('WISENYU_DATA')+'/wisesfrs/'
    kcorr = mrdfits(sfrspath+'wisesfrs_kcorr.fits.gz',1)
    wise = mrdfits(sfrspath+'wisesfrs_wise.fits.gz',1)
    mass = mrdfits(sfrspath+'wisesfrs_mpamass.fits.gz',1)
    ispec = mrdfits(sfrspath+'wisesfrs_mpaispec.fits.gz',1)
    ngal = n_elements(wise)

; initialize the output data structure
    out = {$
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
                 
      fuvmnuv:   -999.0,$ ; FUV-NUV color (rest)
      nuvmr:     -999.0,$ ; NUV-r (rest)
      umr:       -999.0,$ ; u-r (rest)
      Mr:        -999.0,$ ; absolute magnitude
      l1500:     -999.0,$
      beta:      -999.0,$
      k_mass:    -999.0,$ ; from K-correct
      k_coeffs: fltarr(5),$

;     l22_ce01: -999.0,$
;     l22_dh02: -999.0,$
      lir_ce01: -999.0,$
      lir_dh02: -999.0,$
      irx_ce01: -999.0,$
      irx_dh02: -999.0}
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

    cc = iclassification(ispec,snrcut=2.0,ratios=ratios,/doplot)
    out.class = strtrim(ratios.final_class,2)

; get the UV/optical properties
    out.fuvmnuv = kcorr.k_galex_absmag[0]-kcorr.k_galex_absmag[1]
    out.nuvmr = kcorr.k_galex_absmag[1]-kcorr.k_ugrizjhk_absmag[2]
    out.umr = kcorr.k_ugrizjhk_absmag[0]-kcorr.k_ugrizjhk_absmag[2]
    out.Mr = kcorr.k_ugrizjhk_absmag[2]
    out.k_mass = kcorr.k_mass
    out.k_coeffs = kcorr.k_coeffs

    out.beta = 2.32*out.fuvmnuv-2 ; from Hao+11
    out.l1500 = alog10(1500D*kcorr.k_uvflux[0]*4.0*!dpi*3.085678D19^2/3.826D33)
    out.l1500 = alog10(1500D*kcorr.k_uvflux[0]*(dluminosity(kcorr.z,/cm)/3.085678D19)^2/3.826D33)
    out.l1500 = alog10(1500D*kcorr.k_uvflux[0]*4.0*!dpi*dluminosity(kcorr.z,/cm)^2/3.826D33)
    

stop    
    
; finally get the IR properties    
    splog, 'Getting IR properties'
    wise_to_maggies, wise, maggies, ivarmaggies

    out.lir_ce01 = alog10(wise2lir(out.z,maggies,ivarmaggies,/chary,/no12micron,debug=debug))
    out.lir_dh02 = alog10(wise2lir(out.z,maggies,ivarmaggies,/dale,/no12micron,debug=debug))
    out.irx_ce01 = out.lir_ce01-out.l1500
    out.irx_dh02 = out.lir_dh02-out.l1500
    
    outfile = sfrspath+'wisesfrs.fits'
    im_mwrfits, out, outfile, clobber=clobber
    
stop    

return
end
    
