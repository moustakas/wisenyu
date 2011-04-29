pro wiseplots_lir
; jm11apr27ucsd - build the paper plots

    sfrspath = wise_path(/wisesfrs)
    lir = mrdfits(sfrspath+'wisesfrs_lir.fits.gz',1)
    ised = mrdfits(sfrspath+'wisesfrs_isedfit.fits.gz',1)

    psfile = sfrspath+'wisesfrs.ps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=5.0

; IRX vs FUV-NUV
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[-0.2,2], yrange=[-1,3], xtitle='FUV - NUV', $
      ytitle='IRX = log [L(IR) / L(1500 \AA)]'
    djs_oplot, lir.fmnuv, lir.irx_dh02, psym=6, sym=0.3, color=im_color('dodger blue',101)
    legend, 'Hao+11', /right, /bottom, line=0, box=0, pspacing=1.9
    
    fmnuvaxis = range(0.0,2.0,100)
    sfuv = 3.83
    afuv = 0.46
    fmnuv_int = 0.022
    
    djs_oplot, fmnuvaxis, alog10(10^(0.4*sfuv*(fmnuvaxis-fmnuv_int))-1)-alog10(afuv), $
      thick=6
    
;; IRX vs mass    
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      xrange=[8.0,12.0], yrange=[-1,3], xtitle='log (M/M_{'+sunsymbol()+'})', $
;      ytitle='IRX'
;    djs_oplot, ised.mass_avg, lir.irx_dh02, psym=6, sym=0.3;, color=im_color('dodger blue',101)
        
;; L(IR) vs redshift
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      xrange=[0.0,0.053], yrange=[8,12.0], xtitle='Redshift', $
;      ytitle='log (L_{IR}/L_{'+sunsymbol()+'})'
;    djs_oplot, lir.zdist, alog10(lir.lir_dh02), psym=6, sym=0.3;, color=im_color('dodger blue',101)
;    
;; L(IR) vs stellar mass
;;   hogg_scatterplot, ised.mass_avg, alog10(lir.lir_dh02), position=pos, $
;;     xsty=1, ysty=1, xrange=[8.0,12], yrange=[8,12.2], xtitle=mfplot_masstitle(), $
;;     ytitle='log (L_{IR}/L_{'+sunsymbol()+'})', /internal, /outlier, $
;;     outcolor=djs_icolor('grey'), nxpix=30, nypix=30, levels=[0.1,0.25,0.5,0.75,0.9]
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      xrange=[8.0,12], yrange=[8,12.0], xtitle='log (M/M_{'+sunsymbol()+'})', $
;      ytitle='log (L_{IR}/L_{'+sunsymbol()+'})'
;    djs_oplot, ised.mass_avg, alog10(lir.lir_dh02), psym=6, sym=0.3;, color=im_color('dodger blue',101)
;;   djs_oplot, ised.mass_avg, alog10(lir.lir_ce01), psym=6, sym=0.1, color=im_color('firebrick',101)

    im_plotconfig, psfile=psfile, /psclose, /pdf;gzip, /


return
end
    
