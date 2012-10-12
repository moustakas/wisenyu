function wisesfrs_restore, info
; internal function to rebuild the best-fitting SED in the observed 
; frame  

    vname = 'default.nolines'
    light = 2.99792458D18       ; speed of light [A/s]
    ngal = n_elements(info)

    dale = wisenyu_read_02dale()
    chary = wisenyu_read_01chary()

    k_load_vmatrix, vmatrix, lambda, vfile=vfile, $
      lfile=lfile, vpath=vpath, vname=vname
    wave = k_lambda_to_centers(lambda)

    kkeep = where(wave lt 1D5)
    dkeep = where(dale.wave gt 3D4)
    ckeep = where(chary.wave gt 3D4)
    
    model = {wave: wave[kkeep], flux: wave[kkeep]*0.0, $
      dale_wave: dale.wave[dkeep], dale_flux: dale.wave[dkeep]*0.0, $
      chary_wave: chary.wave[ckeep], chary_flux: chary.wave[ckeep]*0.0}
    model = replicate(temporary(model),ngal)

    flux = vmatrix#info.k_coeffs
    model.flux = flux[kkeep,*]
    dlum = dluminosity(info.z,/cm)
    for igal = 0L, ngal-1L do begin
       zplus1 = (1.0+info[igal].z)
       model[igal].wave = model[igal].wave*zplus1 ; [A]
       model[igal].flux = model[igal].flux/zplus1

       model[igal].dale_wave = model[igal].dale_wave*zplus1
       model[igal].chary_wave = model[igal].chary_wave*zplus1

       model[igal].dale_flux = interpolate(dale.flux[dkeep,*],$
         info[igal].modelindx_dh02)/(4.0*!dpi*dlum[igal]^2)/zplus1
       model[igal].chary_flux = interpolate(chary.flux[ckeep,*],$
         info[igal].modelindx_ce01)/(4.0*!dpi*dlum[igal]^2)/zplus1
    endfor

    model.flux = -2.5*alog10((model.flux*model.wave^2/light)>1D-50)-48.6
    model.dale_flux = -2.5*alog10((model.dale_flux*model.dale_wave^2/light)>1D-50)-48.6
    model.chary_flux = -2.5*alog10((model.chary_flux*model.chary_wave^2/light)>1D-50)-48.6
    
return, model
end

pro qaplot_wisesfrs, info, psfile=psfile, clobber=clobber
    
    ngal = n_elements(info)
    if (ngal eq 0L) then begin
       doc_library, 'qaplot_wisesfrs'
       return
    endif

; restore the filters and best-fitting models
    in_filterlist = [galex_filterlist(),sdss_filterlist(),$
      twomass_filterlist(),wise_filterlist(),(wise_filterlist())[0:1]]
    in_filtinfo = im_filterspecs(filterlist=in_filterlist)
    in_nfilt = n_elements(in_filterlist)

    model = wisesfrs_restore(info)

    if tag_exist(info,'galaxy') then galaxy = info.galaxy else begin
       fmt = '(I'+string(5L,format='(I0)')+'.'+string(5L,format='(I0)')+')'
       galaxy = 'Galaxy_'+string(lindgen(ngal),format=fmt)
    endelse

; generate the QA-plot
    xtitle1 = 'Observed Wavelength (\AA)'
    xtitle2 = 'Rest Wavelength (\AA)'
    ytitle1 = 'm_{AB}'

    if (n_elements(psfile) ne 0) then begin
       if file_test(psfile+'.gz',/regular) and $
         (keyword_set(clobber) eq 0) then begin
          splog, 'Output file '+psfile+' exists; use /CLOBBER'
          return
       endif
    endif

    im_plotconfig, 8, pos, psfile=psfile, ymargin=[0.8,1.1]
    for igal = 0L, ngal-1L do begin
       if ((igal mod 10) eq 0) then print, igal+1L, ngal, string(13b), $
         format='("Building QAplot for galaxy ",I0,"/",I0,A10,$)'
       
; redshift and best-fit model
       zobj = info[igal].z
       wave = model[igal].wave    ; [A]
       flux = model[igal].flux    ; [AB]

       maggies = [info[igal].k_maggies,info[igal].wise_maggies[2:3]]
       ivarmaggies = [info[igal].k_ivarmaggies,info[igal].wise_ivarmaggies[2:3]]
       
       dale_wave = model[igal].dale_wave ; [A]
       dale_flux = model[igal].dale_flux ; [AB]
       chary_wave = model[igal].chary_wave   ; [A]
       chary_flux = model[igal].chary_flux   ; [AB]

; make the plot       
       xrange1 = [915,5D6]
;      xrange1 = [min(in_filtinfo.weff-1.3*in_filtinfo.fwhm),$
;        max(in_filtinfo.weff+2.1*in_filtinfo.fwhm)]
;      xrange1[0] = xrange1[0]>(500.0*(1+zobj))
;      xrange1[1] = xrange1[1]<max(wave)

       get_element, wave, xrange1, xx
       get_element, dale_wave, xrange1, dale_xx
       yrange = [max(flux[xx[0]:xx[1]]),min(dale_flux[dale_xx[0]:dale_xx[1]])*0.95]
;      bestmab = -2.5*alog10(info[igal].k_bestmaggies)

       djs_plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
         xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, /xlog, $
         position=pos           ; ymargin=[4,3]       ;, ;, xtickinterval=1000
       axis, /xaxis, xsty=1, xtitle=textoidl(xtitle2), xrange=xrange2, /xlog

       djs_oplot, wave, flux, line=0, color='grey'
       djs_oplot, dale_wave, dale_flux, line=0, color='blue'
       djs_oplot, chary_wave, chary_flux, line=0, color='red'
;      djs_oplot, in_filtinfo.weff, bestmab, $
;        psym=symcat(6,thick=6), symsize=2.5
       legend, ['K-correct','DH02','CE01'], /right, /bottom, box=0, $
         margin=0, pspacing=1.8, line=0, thick=4, $
         color=djs_icolor(['grey','blue','red'])
       
; overplot the data; distinguish between three different cases, based
; on the input photometry
       used = where((maggies gt 0.0) and (ivarmaggies gt 0.0),nused)
       notused = where((maggies gt 0.0) and (ivarmaggies eq 0.0),nnotused)
       nodata = where((maggies eq 0.0) and (ivarmaggies eq 0.0),nnodata)
       upper = where((maggies le 0.0) and (ivarmaggies gt 0.0),nupper)

       if (nused ne 0L) then begin
          mab = maggies2mag(maggies[used],$
            ivar=ivarmaggies[used],magerr=mab_err)
          oploterror, in_filtinfo[used].weff, mab, in_filtinfo[used].fwhm/2.0, $
            mab_err, psym=symcat(16), symsize=2.0, color=djs_icolor('dark green'), $
            errcolor=djs_icolor('dark green'), errthick=!p.thick
       endif

       if (nnotused ne 0L) then begin
          mab = maggies2mag(maggies[notused])
          oploterror, in_filtinfo[notused].weff, mab, in_filtinfo[notused].fwhm/2.0, $
            mab*0.0, psym=symcat(4,thick=6.0), symsize=3.0, color=djs_icolor('red'), $
            errcolor=djs_icolor('red'), errthick=!p.thick
       endif

       if (nupper ne 0) then begin
          mab = maggies2mag(1.0/sqrt(ivarmaggies[upper])/2.0) ; 1-sigma
          oploterror, in_filtinfo[upper].weff, mab, in_filtinfo[upper].fwhm/2.0, $
            mab*0.0, psym=symcat(18), symsize=3.0, color=djs_icolor('blue'), $
            errcolor=djs_icolor('blue'), errthick=!p.thick
       endif

; overlay the observed *output* photometry
       if tag_exist(info,'outmaggies_obs') then begin
          outmab = maggies2mag(info[igal].outmaggies_obs,$
            ivar=info[igal].outivarmaggies_obs,magerr=outmab_err)
          oploterror, out_filtinfo.weff, outmab, out_filtinfo.fwhm/2.0, $
            outmab_err, psym=symcat(14), symsize=2.5, color=djs_icolor('orange'), $
            errcolor=djs_icolor('orange'), errthick=!p.thick
       endif
       
; overlay the synthesized output photometry
       if tag_exist(info,'synth_outmaggies_obs') then begin
          synth_outmab = -2.5*alog10(info[igal].synth_outmaggies_obs)
          djs_oplot, out_filtinfo.weff, synth_outmab, $
            psym=symcat(5,thick=6), symsize=2.5, color='blue'
       endif

;; legend
;       label = textoidl([strtrim(repstr(galaxy[igal],'_',' '),2),$
;         'z = '+string(zobj,format='(F6.4)'),'\chi^{2} = '+$
;         strtrim(string(abs(info[igal].k_chi2),format='(F12.2)'),2),$
;         'log (M/M'+sunsymbol()+') = '+$
;         strtrim(string(info[igal].k_mass,format='(F12.2)'),2)])
;       legend, label, /left, /top, box=0, spacing=1.5, charsize=1.4

       if (n_elements(psfile) eq 0) then cc = get_kbrd(1)
    endfor       
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
return
end
