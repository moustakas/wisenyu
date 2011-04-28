;+
; NAME:
;   WISENYU_READ_01CHARY()
;
; PURPOSE:
;   Read the Chary & Elbaz (2001) infrared models and put them into a
;   standard format.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;   data - structure containing the following fields:
;     lir - total infrared luminosity (8-1000 micron) [NMODEL, L_sun] 
;     l24 - Spitzer/24-micron luminosity [NMODEL, L_sun]
;     l22 - WISE/22-micron luminosity [NMODEL, L_sun]
;     wave - wavelength vector [NPIX, A]
;     flux - flux density [NPIX,NMODEL, erg/s/A]
;
; COMMENTS:
;   Both the Chary & Elbaz models *and* the Dale & Helou (2002) models
;   were downloaded from http://www.its.caltech.edu/~rchary.
; 
; MODIFICATION HISTORY:
;   J. Moustakas, 2010 Feb 17, UCSD
;-

function chary_qpint1d_func, x, wave=wave, flux=flux
return, interpol(flux,wave,x)
end
 
function wisenyu_read_01chary
    
    path = getenv('WISENYU_DIR')+'/data/'
    if (file_test(path,/dir) eq 0) then begin
       splog, 'Data path '+path+' does not exist!'
       return, -1
    endif

; restore Ranga-Ram's IDL save sets; rrcsedsforevol.save
; contains the SEDs themselves, while templatelumatirbands.save
; contains Ranga-Ram's estimates of L(IR), etc.; my values
; match his pretty well except at low infrared luminosity where the
; discrepancies are up to 0.1 dex
    restore, path+'rrcsedsforevol.save'
;   restore, path+'templatelumatirbands.save' 

    filt24 = 'spitzer_mips_24.par'
    filt22 = 'wise_w4.par'
    weff24 = (k_lambda_eff(filterlist=filt24))[0]
    weff22 = (k_lambda_eff(filterlist=filt22))[0]
    
    lsun = 3.826D33
    light = 2.9979246D18

    sz = size(nulnuinlsun,/dim)
    nmodel = sz[0]
    npix = sz[1]

    wave = lambda*1D4 ; [A]
    data = {$
      lir:  dblarr(nmodel),$
      l24:  dblarr(nmodel),$ ; spitzer
      l22:  dblarr(nmodel),$ ; wise
      wave:           wave,$
      flux: dblarr(npix,nmodel)}
    
; convert the models from [nu*L_nu/L_sun] to [erg/s/A]
    wave_edges = k_lambda_to_edges(wave)
    for ii = 0, nmodel-1 do begin
       fnu = lsun*nulnuinlsun[ii,*]*(wave/light) ; [erg/s/Hz]
       data.flux[*,ii] = fnu*(light/wave^2)      ; [erg/s/A]
; integrate to get L(IR)=L(8-1000)
       data.lir[ii] = wisenyu_qpint1d('chary_qpint1d_func',8D4,1000D4,$
         functargs={wave:wave,flux:data.flux[*,ii]})/lsun
; convolve the 22- and 24-micron filter curves to get L(22) L(24); we
; need the factor of 10^40 to trick k_project_filters() from running
; into numerical problems
       f22 = k_project_filters(wave_edges,data.flux[*,ii]/1D40,filterlist=filt22)
       f24 = k_project_filters(wave_edges,data.flux[*,ii]/1D40,filterlist=filt24)
       data.l22[ii] = f22*10.0^(-0.4*48.6)*1D40*(light/weff22^2)*weff22/lsun
       data.l24[ii] = f24*10.0^(-0.4*48.6)*1D40*(light/weff24^2)*weff24/lsun
    endfor

return, data
end
