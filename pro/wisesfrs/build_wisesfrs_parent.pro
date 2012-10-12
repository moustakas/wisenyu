pro build_wisesfrs_parent, clobber=clobber
; jm11apr19ucsd - build the parent sample for the WISE SFR calibration
; paper

; basic sample definitions    
    sfrspath = getenv('WISENYU_DATA')+'/wisesfrs/'
    splog, file=sfrspath+'build_wisesfrs_parent.log'

    sample = 'dr72' & letter = 'bsafe' & poststr = '0' ; everything!
    vagcpath = getenv('VAGC_REDUX')+'/'
    lsspath = getenv('LSS_REDUX')+'/'+sample+'/'+letter+'/'+poststr+'/'

    zmin = 0.01
    zmax = 0.25
    muvcut = 23.0
    m22lim = 14.45
;   m22lim = (wise_flimit())[3]

    raxis = range(14.5,17.6,50) ; for the QAplot
    
; apply the tighter redshift cuts up front    
    allpost = read_vagc_garching(sample=sample,$
      letter=letter,poststr=poststr,/postlss)
    zcut = where(allpost.z ge zmin and allpost.z le zmax)
    allpost = allpost[zcut]

    allwise = mrdfits(vagcpath+'object_wise.fits.gz',$
      row=allpost.object_position,1)
    allgalex = mrdfits(vagcpath+'object_galex_gr6.fits.gz',$
      row=allpost.object_position,1)
    nall = n_elements(allpost)
    
    im_galex_to_maggies, allgalex, gmaggies, givarmaggies
    wise_to_maggies, allwise, wmaggies, wivarmaggies

    good22 = where(wmaggies[3,*] gt 0.0,ngood22)
    upper22 = where((wmaggies[3,*] eq 0.0) and (wivarmaggies[3,*] gt 0.0),nupper22)
    m22 = fltarr(nall)
    m22[good22] = -2.5*alog10(wmaggies[3,good22])
    m22[upper22] = -2.5*alog10(1.0/sqrt(wivarmaggies[3,upper22])/2.0) ; 1-sigma
    m22upper = long(m22*0)
    m22upper[upper22] = 1
    
; build the sample piece by piece and log the output
    inwise = (allwise.wise_tag gt -1) and (wmaggies[3,*] ne -1.0) ; very soft cut on 22-microns!
    ingalex = (allwise.wise_tag gt -1) and (wmaggies[3,*] ne -1.0) and $
      (allgalex.nuv_weight gt 1E3) and (allgalex.fuv_weight gt 1E3)
    galexgood = (allwise.wise_tag gt -1) and (wmaggies[3,*] ne -1.0) and $
      (allgalex.nuv_weight gt 1E3) and (allgalex.fuv_weight gt 1E3) and $
      (gmaggies[0,*] ge 10^(-0.4*muvcut)) and (gmaggies[1,*] ge 10^(-0.4*muvcut))
    m22good = (allwise.wise_tag gt -1) and (wmaggies[3,*] ne -1.0) and $
      (allgalex.nuv_weight gt 1E3) and (allgalex.fuv_weight gt 1E3) and $
      (gmaggies[0,*] ge 10^(-0.4*muvcut)) and (gmaggies[1,*] ge 10^(-0.4*muvcut)) and $
      (m22 gt 0) and (m22 le m22lim)
    m22detect = (allwise.wise_tag gt -1) and (wmaggies[3,*] ne -1.0) and $
      (allgalex.nuv_weight gt 1E3) and (allgalex.fuv_weight gt 1E3) and $
      (gmaggies[0,*] ge 10^(-0.4*muvcut)) and (gmaggies[1,*] ge 10^(-0.4*muvcut)) and $
      (m22 gt 0) and (m22 le m22lim) and (m22upper eq 0)
    keep = where(galexgood,ngal)
    
    splog, 'VAGC = '+sample+letter+poststr
    splog, '  Redshift limits = ', zmin, zmax
    splog, '  Magnitude limits = 14.5 < r < 17.6'
    splog, '  Ngal = ', nall
    splog, ' '
    splog, 'WISE footprint ', total(inwise), float(nall), total(inwise)/float(nall)
    splog, 'WISE+GALEX footprint ', total(ingalex), total(inwise), total(ingalex)/total(inwise)
    splog, '  NUV=FUV<23 ', total(galexgood), total(ingalex), total(galexgood)/total(ingalex)
    splog, '  m(22)<14.5', total(m22good), total(galexgood), total(m22good)/total(galexgood)
    splog, '  m(22)<14.5 detections', total(m22detect), total(galexgood), $
      total(m22detect)/total(galexgood), total(m22detect), total(m22good), $
      total(m22detect)/total(m22good)
    
; show the effect of our UV flux cuts 
;   ww = where(cut2)
;   plot, allpost[ww].m, -2.5*alog10(gmaggies[0,ww])-allpost[ww].m, $
;     psym=3, xsty=3, ysty=3, yrange=[0,10]
;   djs_oplot, raxis, muvcut-raxis, line=0
;   djs_oplot, raxis, 25-raxis, line=5
    
;; establish the flux limits    
;    keep = where((allpost.z gt zmin) and (allpost.z lt zmax) and $
;      (allwise.wise_tag gt -1) and $ ; in the WISE footprint
;      (allgalex.nuv_weight gt 1E3) and $
;      (wmaggies[3,*] ne -1.0) and $
;      (gmaggies[0,*] ge 10^(-0.4*muvcut)) and $
;      (gmaggies[1,*] ge 10^(-0.4*muvcut)))
;    det = where(wmaggies[3,keep] gt 0.0)
;    up = where(wmaggies[3,keep] eq 0.0 and wivarmaggies[3,keep] gt 0)
;    m22_det = -2.5*alog10(wmaggies[3,keep[det]])
;    m22_up = -2.5*alog10(1.0/sqrt(wivarmaggies[3,keep[up]])/2.0)
;    mr = allpost[keep].m
;    
;    raxis = range(14.5,17.6,50)
;
;    djs_plot, mr[det], mr[det]-m22_det, ps=3, xsty=3, ysty=3, yr=[0,6]
;    djs_oplot, mr[up], mr[up]-m22_up, ps=3, color='red'
;;   djs_oplot, mr[up], mr[up]-(m22_up+1), ps=3, color='blue'
;    djs_oplot, mr[det], mr[det]-m22_det, ps=6, sym=0.05, color='yellow'
;    djs_oplot, raxis, raxis-lim[3], line=0

;   plot, allpost.z, -2.5*alog10(wmaggies[3,*]), psym=3, yr=[10,15]
;   plot, allpost.z, -2.5*alog10(gmaggies[0,*]), psym=3, yr=[10,30]

; read the rest of the data structures we need
    post = allpost[keep]
    wise = allwise[keep]
    galex = allgalex[keep]

    sdss = mrdfits(vagcpath+'object_sdss_imaging.fits.gz',row=post.object_position,1)
    twomass = mrdfits(vagcpath+'object_twomass.fits.gz',row=post.object_position,1)
    mpamass = mrdfits(lsspath+'mpamassoh.dr72bsafe0.fits.gz',row=zcut[keep],1)
    mpaispec = mrdfits(lsspath+'ispecline.dr72bsafe0.fits.gz',row=zcut[keep],1)

; compute K-corrections using various bands; kcorr1=all bands;
; kcorr2=no wise; kcorr3=nowise,no 2mass
    im_galex_to_maggies, galex, gmaggies, givarmaggies
    sdss_to_maggies, calib=sdss, smaggies, sivarmaggies, flux='cmodel'
    twomass_to_maggies, twomass, tmaggies, tivarmaggies
    wise_to_maggies, wise, wmaggies, wivarmaggies
    filterlist = [galex_filterlist(),sdss_filterlist(),twomass_filterlist(),$
      (wise_filterlist())[0:1]]

    kcorr1 = wisesfrs_do_kcorrect(post.z,[gmaggies,smaggies,tmaggies,wmaggies[0:1,*]],$
      [givarmaggies,sivarmaggies,tivarmaggies,wivarmaggies[0:1,*]],filterlist=filterlist)
    kcorr2 = wisesfrs_do_kcorrect(post.z,[gmaggies,smaggies,tmaggies,wmaggies[0:1,*]],$
      [givarmaggies,sivarmaggies,tivarmaggies,0*wivarmaggies[0:1,*]],filterlist=filterlist)
    kcorr3 = wisesfrs_do_kcorrect(post.z,[gmaggies,smaggies,tmaggies,wmaggies[0:1,*]],$
      [givarmaggies,sivarmaggies,0*tivarmaggies,0*wivarmaggies[0:1,*]],filterlist=filterlist)

; write out!
    im_mwrfits, struct_addtags(post,kcorr1), sfrspath+'wisesfrs_kcorr_all.fits', clobber=clobber
    im_mwrfits, struct_addtags(post,kcorr2), sfrspath+'wisesfrs_kcorr_nowise.fits', clobber=clobber
    im_mwrfits, struct_addtags(post,kcorr3), sfrspath+'wisesfrs_kcorr_no2mass.fits', clobber=clobber
    im_mwrfits, wise, sfrspath+'wisesfrs_wise.fits', clobber=clobber
    im_mwrfits, mpamass, sfrspath+'wisesfrs_mpamass.fits', clobber=clobber
    im_mwrfits, mpaispec, sfrspath+'wisesfrs_mpaispec.fits', clobber=clobber
    splog, /close
    
stop    

return
end
    
