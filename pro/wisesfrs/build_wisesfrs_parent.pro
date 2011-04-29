pro build_wisesfrs_parent, clobber=clobber
; jm11apr19ucsd - build the VAGC parent sample for the WISE SFR
; calibration paper

    sample = 'dr72' & letter = 'bsafe' & poststr = '0' ; everything!
    vagcpath = getenv('VAGC_REDUX')+'/'
    lsspath = getenv('LSS_REDUX')+'/'+sample+'/'+letter+'/'+poststr+'/'
    
    allpost = read_vagc_garching(sample=sample,$
      letter=letter,poststr=poststr,/postlss)
    allwise = mrdfits(vagcpath+'object_wise.fits.gz',$
      row=allpost.object_position,1)
    allgalex = mrdfits(vagcpath+'object_galex_gr6.fits.gz',$
      row=allpost.object_position,1)

; sample cuts: (1) in the VAGC at 0.0; (2) m22<14.5; (3) NUV=FUV<23
    zmin = 0.01
    zmax = 0.25
    m22cut = 14.5
    muvcut = 23.0
    
    wise_to_maggies, allwise, wmaggies, wivarmaggies
    im_galex_to_maggies, allgalex, gmaggies, givarmaggies

;   plot, allpost.z, -2.5*alog10(wmaggies[3,*]), psym=3, yr=[10,15]
;   plot, allpost.z, -2.5*alog10(gmaggies[0,*]), psym=3, yr=[10,30]
    keep = where((allpost.z gt zmin) and (allpost.z lt zmax) and $
      (allwise.wise_tag gt -1) and $ ; in the WISE footprint
      (allgalex.nuv_weight gt 1E3) and $
      (gmaggies[0,*] ge 10^(-0.4*muvcut)) and $
      (gmaggies[1,*] ge 10^(-0.4*muvcut)),ngal)
;   keep = where((allpost.z gt zmin) and (allpost.z lt zmax) and $
;     (wmaggies[3,*] ge 10^(-0.4*m22cut)) and $
;     (gmaggies[0,*] ge 10^(-0.4*muvcut)) and $
;     (gmaggies[1,*] ge 10^(-0.4*muvcut)),ngal)
;   keep2 = where((allpost.z gt zmin) and (allpost.z lt zmax) and $
;     (gmaggies[0,*] ge 10^(-0.4*muvcut)) and $
;     (gmaggies[1,*] ge 10^(-0.4*muvcut)))

;; test the selection - it's mostly due to the 22-micron cut    
;    im_plothist, allpost.z, bin=0.002                               
;    im_plothist, allpost[keep].z, bin=0.002, /over, norm=0.1, color='blue'
;    im_plothist, allpost[keep2].z, bin=0.002, /over, color='red'

; now read the rest of the data we need; need stellar masses at some
; point 
    post = allpost[keep]
    wise = allwise[keep]
    galex = allgalex[keep]

    sdss = mrdfits(vagcpath+'object_sdss_imaging.fits.gz',row=post.object_position,1)
    twomass = mrdfits(vagcpath+'object_twomass.fits.gz',row=post.object_position,1)
    mpamass = mrdfits(lsspath+'mpamassoh.dr72bsafe0.fits.gz',row=keep,1)
    mpaispec = mrdfits(lsspath+'ispecline.dr72bsafe0.fits.gz',row=keep,1)

; compute K-corrections
    im_galex_to_maggies, galex, gmaggies, givarmaggies
    sdss_to_maggies, calib=sdss, smaggies, sivarmaggies, flux='cmodel'
    twomass_to_maggies, twomass, tmaggies, tivarmaggies
    filterlist = [galex_filterlist(),sdss_filterlist(),twomass_filterlist()]

    kcorr = wisesfrs_do_kcorrect(post.z,[gmaggies,smaggies,tmaggies],$
      [givarmaggies,sivarmaggies,tivarmaggies],filterlist=filterlist)

; write out!
    sfrspath = getenv('WISENYU_DATA')+'/wisesfrs/'
    im_mwrfits, struct_addtags(post,kcorr), sfrspath+'wisesfrs_kcorr.fits', clobber=clobber
    im_mwrfits, wise, sfrspath+'wisesfrs_wise.fits', clobber=clobber
    im_mwrfits, mpamass, sfrspath+'wisesfrs_mpamass.fits', clobber=clobber
    im_mwrfits, mpaispec, sfrspath+'wisesfrs_mpaispec.fits', clobber=clobber
    
stop    

return
end
    
