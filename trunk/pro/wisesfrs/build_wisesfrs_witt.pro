pro build_wisesfrs_witt, debug=debug
; jm11apr29ucsd - build the A(lambda) vs lambda look-up table using
; the Witt & Gordon (2000) dust models and the flux-ratio method of
; Gordon et al. (2000)

; choose three metallicities four ages, and five a_d values 
    read_many_fit, dust_info, sf_type='c' ; continuous star formation
    metal = [-0.4,0.0,0.4]
    ages = [10.0,50.0,100.0,1000.0]
    ad = [0.01,0.25,0.5,0.75,0.95]
    ndust = n_elements(metal)*n_elements(ages)*n_elements(ad)

    irx = range(0.1,1D6,1D3,/log)
    nirx = n_elements(irx)
    witt_dust = {irx: float(irx), a1500: float(irx*0.0), $
      a1500_grid: fltarr(nirx,ndust), params: strarr(ndust)}

    counter = 0
    for im = 0, n_elements(metal)-1 do begin
       for ia = 0, n_elements(ages)-1 do begin
          for id = 0, n_elements(ad)-1 do begin
             witt_dust.params[counter] = 'Z'+repstr(string(metal[im],format='(F4.1)'),' ','+')+$
               '_Age'+string(ages[ia],format='(I4.4)')+'_ad'+string(ad[id],format='(F4.2)')
             for jj = 0L, nirx-1L do witt_dust.a1500_grid[jj,counter] = fr2atten(irx[jj],$
               '',1500.0,ages[ia],metal[im],ad[id],dust_info,/silent)
             counter++
          endfor
       endfor
    endfor

; compute the media A(1500)-IRX relation and write out
    witt_dust.a1500 = total(witt_dust.a1500_grid,2)/float(ndust)
;      for jj = 0L, nirx-1L do witt_dust.a1500[jj] = djs_median(witt_dust.a1500_grid[jj,*])

    outfile = getenv('WISENYU_DATA')+'/wisesfrs/wisesfrs_witt.fits'
    im_mwrfits, witt_dust, outfile, /clobber
    
    if keyword_set(debug) then begin
       djs_plot, 1.0+irx, witt_dust.a1500_grid[*,0], xr=[1,1D3], yr=[0,6], /xlog, ps=3
       for jj = 1, ndust-1 do djs_oplot, 1.0+irx, witt_dust.a1500_grid[*,jj], ps=3
       djs_oplot, 1.0+irx, witt_dust.a1500, color='red', thick=3
    endif

return
end
