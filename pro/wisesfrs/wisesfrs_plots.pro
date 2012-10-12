pro wisesfrs_plots
; jm11apr27ucsd - build the paper plots

    sfrspath = getenv('WISENYU_DATA')+'/wisesfrs/'
    texpath = getenv('WISENYU_DIR')+'/tex/wisesfrs/'

; read the output of BUILD_WISESFRS    
    wise = mrdfits(sfrspath+'wisesfrs.fits.gz',1)
    witt = mrdfits(sfrspath+'wisesfrs_witt.fits.gz',1)
    ngal = n_elements(wise)

; ---------------------------------------------------------------------------
; UV-optical vs mass selection diagram; show the number contours of
; the parent sample and the fraction of 22-micron sources in each bin
; of color and mass
    m22all = where(wise.m22good)
    m22good = where(wise.m22good and (wise.m22upper eq 0))

    allmass = wise.k_mass
    allnuvmr = wise.nuvmr
    mass = allmass[m22all]
    nuvmr = allnuvmr[m22all]
    
    nbins = [40,40]
    hmin = [8,0.5]
    hmax = [12.0,7.0]
    binsz = (hmax-hmin)/float(nbins-1)
    
    nall = hist_nd(transpose([[allmass],[allnuvmr]]),nbins=nbins,min=hmin,max=hmax)
    nm22 = hist_nd(transpose([[mass],[nuvmr]]),nbins=nbins,min=hmin,max=hmax)
    m22frac = float(nm22)/(float(nall)+(nall eq 0))*(nall ne 0)
    low = where(nall lt 10) & m22frac[low] = -1.0
    
    cumindex = reverse(sort(nall))
    cumimage = fltarr(nbins)
    cumimage[cumindex] = total((nall)[cumindex],/cumu)
    fracobj = cumimage/total(nall)

    massaxis = findgen(nbins[0])*binsz[0]+hmin[0]
    coloraxis = findgen(nbins[1])*binsz[1]+hmin[1]

    nmock = 300
    mockmass = randomu(seed,nmock,nmock)*(hmax[0]-hmin[0])+hmin[0]
    mockcolor = randomu(seed,nmock,nmock)*(hmax[1]-hmin[1])+hmin[1]
    mockfrac = interpolate(m22frac,findex(massaxis,mockmass),$
      findex(coloraxis,mockcolor),missing=-1.0)

    mass2daxis = massaxis#(coloraxis*0.0+1.0)
    color2daxis = transpose(coloraxis#(massaxis*0.0+1.0))
    
    gd = where(mockfrac gt 0)

    mincol = 50
    levels = [0.10,0.25,0.5,0.75,0.9,0.95]
    cannotation = string(levels,format='(F4.2)')

    psfile = texpath+'nuvmr_mass.eps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=5.0
    
    loadct, 3 ; 15
    wise_hogg_scatterplot, mockmass[gd], mockcolor[gd], weight=mockfrac[gd], position=pos, $
      xsty=1, ysty=1, xrange=[hmin[0],hmax[0]], yrange=[hmin[1],hmax[1]], xnpix=nbins[0], $
      ynpix=nbins[1], exp=2.0, levels=levels, /nocontour, $
      xtitle=textoidl('log (M / M_{'+sunsymbol()+'})'), /internal, $
      ytitle=textoidl('^{0.1}(NUV - r)'), darkest=mincol
    contour, fracobj, mass2daxis, color2daxis, levels=levels, $
      /overplot, c_annotation=cannotation
    loadct, 0
    
    im_plotconfig, psfile=psfile, /psclose

; ---------------------------------------------------------------------------
; IRX vs FUV-NUV
    psfile = texpath+'irx.eps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.8,1.1], $
      xmargin=[1.3,1.2], width=6.0, height=4.7

    agn = where(wise.isagn eq 1 and wise.m22good and (wise.m22upper eq 0))
    sf = where(wise.isagn eq 0 and wise.m22good and (wise.m22upper eq 0))
    hogg_scatterplot, wise[sf].fuvmnuv, wise[sf].irx, position=pos, $
      xsty=9, ysty=9, xrange=[-0.3,1.8], yrange=[-0.5,3.5], $
      xtitle=textoidl('^{0.1}(FUV - NUV)'), $
      ytitle=textoidl('IRX = log [L(IR) / L(1500 \AA)]'), $ 
      /outliers, /internal, levels=[0.25,0.5,0.75,0.9], $
      xnpix=30, ynpix=30, /nogrey, outcolor=im_color('dodger blue',101)
    axis, /xaxis, xsty=1, xrange=2.32*!x.crange-2, xtitle=textoidl('\beta')
    linterp, witt.irx, witt.a1500, 10^!y.crange, yrange2
    axis, /yaxis, ysty=1, yrange=yrange2, ytitle=textoidl('A(1500) (mag)')
    
    legend, 'Hao+11', /left, /top, line=0, box=0, pspacing=1.9, charsize=1.4
;   legend, ['SF','AGN'], /right, /bottom, 

    plotsym, 1, 1.0
;   djs_oplot, wise[sf2].fuvmnuv, wise[sf2].irx, psym=8
;   djs_oplot, wise[sf].fuvmnuv, wise[sf].irx, psym=symcat(6,thick=5), $
;     sym=0.3, color=im_color('dodger blue',101)
;   djs_oplot, wise[agn].fuvmnuv, wise[agn].irx, psym=symcat(9,thick=5), sym=0.9
;   djs_oplot, wise[agn].fuvmnuv, wise[agn].irx, psym=symcat(16), $
;     sym=0.3, color=im_color('forest green',101)
    
    fmnuvaxis = range(0.05,2.0,100)
    sfuv = 3.83
    afuv = 0.46
    fmnuv_int = 0.022
    
    djs_oplot, fmnuvaxis, alog10(10^(0.4*sfuv*(fmnuvaxis-fmnuv_int))-1)-alog10(afuv), $
     thick=6

    resid = wise[sf].irx-alog10(10^(0.4*sfuv*(wise[sf].fuvmnuv-fmnuv_int))-1)-alog10(afuv)
    
    im_plotconfig, psfile=psfile, /psclose

stop    
    
; ---------------------------------------------------------------------------
; UV-optical-IR color-color diagram
    psfile = texpath+'nuvmr_rmj.eps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=5.0

    xrange = [-0.5,2.1] & yrange=[0.3,7]
    hogg_scatterplot, wise.rmj, wise.nuvmr, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl('^{0.1}(r - J)'), ytitle=textoidl('^{0.1}(NUV - r)'), $
      levels=[0.25,0.5,0.75,0.9,0.975], xnpix=40, ynpix=40, $
      outliers=1, outcolor=im_color('grey',101), /internal, $
      /nogrey;, ccolor=im_color('dodger blue',101)

    djs_oplot, wise.rmj, wise.nuvmr, psym=6, sym=0.2, color='grey'
    w1 = where(wise.irx gt 2.0 and wise.irupper eq 0)
    w2 = where(wise.irx lt 2.0 and wise.irx gt 1.0 and wise.irupper eq 0)
    w3 = where(wise.irx lt 1.0 and wise.irupper eq 0)
    djs_oplot, wise[w2].rmj, wise[w2].nuvmr, psym=6, sym=0.2, color='dark green'
    djs_oplot, wise[w3].rmj, wise[w3].nuvmr, psym=6, sym=0.2, color='blue'
    djs_oplot, wise[w1].rmj, wise[w1].nuvmr, psym=6, sym=0.2, color='red'
    
;   hogg_scatterplot, wise[cut2].ch1mch3, wise[cut2].nuvmr, /noerase, position=pos, $
;     xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;     /nogrey, ccolor=im_color('firebrick',101), cthick=5, $
;     levels=[0.25,0.5,0.75,0.9], /outliers

    im_plotconfig, psfile=psfile, /psclose

; ---------------------------------------------------------------------------
; UV-optical-IR color-color diagram
    psfile = texpath+'nuvmr_ch1mch3.eps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=5.0

    xrange = [-3.0,3.0] & yrange=[0.3,7]
    hogg_scatterplot, wise[cut1].ch1mch3, wise[cut1].nuvmr, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl('[3.4] - [12]'), ytitle=textoidl('^{0.1}(NUV - r)'), $
      levels=[0.25,0.5,0.75,0.9,0.975], xnpix=40, ynpix=40, $
      outliers=1, outcolor=im_color('grey',101), /internal, $
      /nogrey, ccolor=im_color('dodger blue',101)
    hogg_scatterplot, wise[cut2].ch1mch3, wise[cut2].nuvmr, /noerase, position=pos, $
      xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      /nogrey, ccolor=im_color('firebrick',101), cthick=5, $
      levels=[0.25,0.5,0.75,0.9], /outliers

    im_plotconfig, psfile=psfile, /psclose

; ---------------------------------------------------------------------------
; NUV-r vs Mr color-magnitude diagram
    psfile = texpath+'nuvmr.eps'
    im_plotconfig, 0, pos, psfile=psfile, ymargin=[0.6,1.1], $
      xmargin=[1.3,0.4], width=6.8, height=5.0

    xrange = [-15.5,-24.5] & yrange=[0.3,7]
    hogg_scatterplot, wise.mr, wise.nuvmr, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xtitle=textoidl('M_{0.1r}'), ytitle=textoidl('^{0.1}(NUV - r)'), $
      levels=[0.25,0.5,0.75,0.9,0.975], xnpix=50, ynpix=50, $
      outliers=1, outcolor=im_color('grey',101), /internal, $
      /nogrey, ccolor=im_color('dodger blue',101)
    hogg_scatterplot, wise[ir].mr, wise[ir].nuvmr, /noerase, position=pos, $
      xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      /nogrey, ccolor=im_color('firebrick',101), cthick=5, $
      levels=[0.25,0.5,0.75,0.9], /outliers

; inset redshift distribution
    bin = 0.002
    im_plothist, wise.z, bin=bin, xall, yall, omin=0.0, omax=0.25, /noplot
    im_plothist, wise[ir].z, bin=bin, xir, yir, omin=0.0, omax=0.25, /noplot
    yrange = [0.0,max(yall)*1.05]

    im_plothist, wise.z, bin=bin, omin=0.0, omax=0.25, $
      /noerase, position=[0.24,0.68,0.48,0.88], $
      xsty=1, ysty=1, xrange=[0.0,0.2], yrange=yrange, $
      charsize=1.0, xtitle='Redshift', ytitle='Number', $
      xtickinterval=0.1, ytickinterval=500
    im_plothist, wise[ir].z, bin=bin, omin=0.0, omax=0.25, /overplot, $
      /fill, fcolor=djs_icolor('grey')

    im_plotconfig, psfile=psfile, /psclose

stop    
    

return
end
    
