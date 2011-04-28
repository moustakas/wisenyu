;+
; NAME:
;   WISE2LIR()
;
; PURPOSE:
;   Given a redshift and the four bands of *observed* WISE photometry,
;   compute the total infrared luminosity, L(IR)=L(8-1100), based on
;   various infrared infrared SED models.  See COMMENTS for details.
;
; INPUTS: 
;   redshift - redshift for each object [NGAL]
;   maggies - observed 3.4, 4.6, 12, and 22 micron fluxes, as output
;     by WISE_TO_MAGGIES [4,NGAL]
;   ivarmaggies - corresponding inverse variances [4,NGAL]
;
; OPTIONAL INPUTS: 
;   h100 - Hubble constant divided by 100 (default 0.7)
;   omega0 - matter density (default 0.3)
;   omegal - cosmological constant density (default 0.7)
;
; KEYWORD PARAMETERS: 
;   no12micron - do not use the 12-micron flux to constrain the fit 
;   chary - use the Chary & Elbaz (2001) models (default)
;   dale - use the Dale & Helou (2002) models
;   debug - make a QAplot on the screen and wait for a keystroke 
;   zlog - use a logarithmically spaced fiducial redshift interval
;     (default is linear); useful when your sample includes very low
;     redshift objects
;
; OUTPUTS: 
;   lir - infrared luminosity for each object [L_sun, NGAL]
;
; OPTIONAL OUTPUTS:
;   chi2 - chi^2 of the fit [NGAL]
;   model_indx - fractional index number of the best-fitting model 
;     [NGAL] 

; COMMENTS:
;   Note that the 3.6- and 4.6-micron fluxes are *never* used to
;   constrain the SED since they are affected significantly by the
;   stellar light.  The 12-micron flux is used by default unless
;   /NO12MICRON is set.  Also note that the 22-micron flux point is
;   forced to have the minimum inverse variance.
;
;   Could use more error checking on the inputs!
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2011 Apr 19, UCSD
;-

function wise2lir_model_maggies, model, zref, dlum=dlum
; compute model photometry on a fiducial redshift grid    
    nmodel = n_elements(model.lir)
    nzref = n_elements(zref)
; k_projection_table can't deal with the units so we have to
; loop and call k_project_filters    
;   k_projection_table, rmatrix, model.flux, model.wave, $
;     zref, (wise_filterlist())[2:3]
    model_maggies = dblarr(4,nmodel,nzref)
    for ii = 0, nmodel-1 do begin
       for jj = 0, nzref-1 do begin
          model_maggies[*,ii,jj] = k_project_filters(k_lambda_to_edges($
            model.wave*(1.0+zref[jj])),model.flux[*,ii]/(4.0*!dpi*dlum[jj]^2)/$
            (1.0+zref[jj]),filterlist=wise_filterlist())
       endfor
    endfor
return, model_maggies
end

function wise2lir, redshift, maggies, ivarmaggies, h100=h100, $
  omega0=omega0, omegal=omegal, no12micron=no12micron, chi2=chi2, $
  model_indx=model_indx, chary=chary, dale=dale, debug=debug, $
  zlog=zlog

    common com_wise2lir, chary_maggies, dale_maggies
    
    ngal = n_elements(redshift)
    if (ngal eq 0L) then begin
       doc_library, 'wise2lir'
       return, -1
    endif

; cosmology
    if (n_elements(h100) eq 0) then h100 = 0.7
    if (n_elements(omega0) eq 0) then omega0 = 0.3
    if (n_elements(omegal) eq 0) then omegal = 0.7
    
; define the fiducial redshift grid
    zrefmin = min(redshift)
    zrefmax = max(redshift)
    nzz = (ceil((zrefmax-zrefmin)/0.01)>50)<200
    zref = wisenyu_range(zrefmin,zrefmax,nzz,log=zlog)
    dlum = lumdist(zref,H0=h100*100.0,omega_m=omega0,lambda0=omegal,/silent)*3.086D24
    
; read the relevant set of models and build the fiducial table of
; photometry; use Chary & Elbaz (2001) by default
    if (keyword_set(chary) eq 0) and (keyword_set(dale) eq 0) then chary = 1
    
; Chary & Elbaz (2001)    
    if keyword_set(chary) then begin
       model = wisenyu_read_01chary()
       if (n_elements(chary_maggies) eq 0) then $
         chary_maggies = wise2lir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = chary_maggies
    endif
; Dale & Helou (2002)
    if keyword_set(dale) then begin
       model = wisenyu_read_02dale()
       if (n_elements(dale_maggies) eq 0) then $
         dale_maggies = wise2lir_model_maggies(model,zref,dlum=dlum)
       model_maggies_grid = dale_maggies
    endif
    nmodel = n_elements(model.lir)
    npix = n_elements(model.wave)

    filt = wise_filterlist()
    weff = k_lambda_eff(filterlist=filt)
    nfilt = n_elements(filt)
    
; interpolate the models at the redshifts of interest
    zindx = wisenyu_findex(zref,redshift)
    model_maggies = interpolate(model_maggies_grid,zindx)
    idlum = interpolate(dlum,zindx)

; ignore the 3.3- and 4.6-micron fluxes and give the 12- and 22-micron
; bands equal statistical weight
    fitivarmaggies = ivarmaggies
    fitivarmaggies[0:1,*] = 0.0
    fitivarmaggies[2,*] = max(ivarmaggies,dim=1)*(keyword_set(no12micron) eq 0)
    fitivarmaggies[3,*] = max(ivarmaggies,dim=1)

; need to loop, unfortunately...
    chi2floor = 100
    nfine = 20 ; chi^2 oversampling factor

    model_indx = fltarr(ngal)-1.0
    lir = model_indx*0.0-1.0
    chi2 = model_indx*0.0-1.0

    for igal = 0L, ngal-1 do begin
       nmaggies = abs(maggies[*,igal]*1.0D) ; need absolute value to deal with negative fluxes correctly
       nivarmaggies = fitivarmaggies[*,igal]*1.0D
       if (total(nivarmaggies gt 0) ge 1) then begin

; compute and minimize chi2
          vmodelmaggies = model_maggies[*,*,igal]
          vmaggies = rebin(reform(nmaggies,nfilt,1),nfilt,nmodel)
          vivarmaggies = rebin(reform(nivarmaggies,nfilt,1),nfilt,nmodel)
          vchi2 = total(vivarmaggies*(vmaggies-vmodelmaggies)^2,1,/double)

; oversample
          xx = findgen(nmodel)
          minx = min(xx,max=maxx)
          y2 = spl_init(xx,vchi2)
          x2 = findgen(nmodel*nfine+1)/(nmodel*nfine)*(maxx-minx)+minx

; we have to add a floor to chi^2 to prevent negative minima, which
; can occur because of the coarse sampling of the models; use a
; lower-order interpolation scheme if the chi^2 goes negative anyway 
          vchi2 = vchi2+chi2floor 
          chi2fit = spl_interp(xx,vchi2,y2,x2)
          if (min(chi2fit) lt 0.0) then chi2fit = interpol(vchi2,xx,x2,/quad)
          if (min(chi2fit) lt 0.0) then chi2fit = interpol(vchi2,xx,x2)
          minchi2 = min(chi2fit,fitplace)
          if (minchi2 le 0) then message, 'Negative chi^2'

          modelindx1 = x2[fitplace] 
          model_indx[igal] = modelindx1 ; fractional index number
          lir[igal] = interpolate(model.lir,modelindx1)

          minchi2 = minchi2-chi2floor
          chi2[igal] = minchi2
          
; debugging QAplots
          if keyword_set(debug) then begin
             splog, minchi2, modelindx1, lir[igal]

             djs_plot, x2, chi2fit, psym=-6, sym=0.5, xsty=3, ysty=3, /ylog, $
               xrange=modelindx1+5*[-1,1], xtitle='Model Number', ytitle='\chi^{2}'
             djs_oplot, xx, vchi2, psym=-6, sym=3, color='red'
             djs_oplot, modelindx1*[1,1], 10^!y.crange, color='dark green'
             djs_oplot, !x.crange, minchi2*[1,1], color='dark green'
             cc = get_kbrd(1)
             
             modelflam = interpolate(model.flux,modelindx1)/(4.0*!dpi*idlum[igal]^2)/$
               (1.0+redshift[igal])
             good = where(modelflam gt 0.0,ngood)
             modelwave = model.wave[good]*(1.0+redshift[igal])
             modelmab = -2.5*alog10(modelflam[good]*rebin(reform(modelwave,ngood,1),$
               ngood,nmodel)^2/2.9979246D18)-48.6

             djs_plot, weff*(1.0+redshift[igal])/1D4, -2.5*alog10(maggies[*,igal]), $
               psym=7, color='yellow', sym=3, xr=[1,500], /xlog, yr=[max(modelmab),min(modelmab)], $
               xsty=3, ysty=3, xtitle='Observed Wavelength (\AA)', ytitle='Observed AB Magnitude'
             djs_oplot, weff*(1.0+redshift[igal])/1D4, -2.5*alog10(interpolate(vmodelmaggies,modelindx1)), $
               psym=6, sym=3
             djs_oplot, modelwave/1D4, modelmab
             cc = get_kbrd(1)
          endif
       endif 
    endfor

return, lir
end
