pro wise_qaplot_templates

;; SSP models
indir = wisenyu_path()+'/templates/'
outdir = wisenyu_path()+'/qaplots/'
infile = indir + 'bc03_all_with_maggies.fits'

final = mrdfits(infile, 1)
iage = where(final.age ge 5.E9 and final.age le 15.E9 $
         and final.met ge 0.39 and final.met le 1.01, nage)
iz = where(final[0].zvals gt 0.2 and final[0].zvals le 0.4, nz)

zz = 0.00

;; 10 Gyr
tmp_min = min(abs(final[iage].age-10.0E9), imin)
iunique = iage[imin]
wave = final[iunique].lambda/1.E4*(1.+zz)
flux = final[iunique].flux
;; 10 Gyr and 0.4 Z_solar 
print, final[iunique].age, final[iunique].met

filters = (strtrim(final[0].filters,2))[12:15]
k_load_filters, filters, filter_nlambda, filter_lambda, filter_pass
weff =  (primus_lambda_eff(filterlist=filters))/1.E4
tmp_min = min(abs(wave-weff[0]), imin)

flux = -2.5*alog10(flux[imin]*wave[imin]^2)+2.5*alog10(flux*wave^2)

;; maggies: 5 SDSS, 3 2MASS, 4 IRAC, 4 WISE
ch1 = -2.5*alog10(final[iage].maggies[iz,12])+2.5*alog10(final[iage].maggies[iz,12])
ch2 = -2.5*alog10(final[iage].maggies[iz,12])+2.5*alog10(final[iage].maggies[iz,13])
ch3 = -2.5*alog10(final[iage].maggies[iz,12])+2.5*alog10(final[iage].maggies[iz,14])
ch4 = -2.5*alog10(final[iage].maggies[iz,12])+2.5*alog10(final[iage].maggies[iz,15])
npoint = n_elements(ch1)
ch2_mean = (moment(ch2,sdev=ch2_sdev))[0]
ch3_mean = (moment(ch3,sdev=ch3_sdev))[0]
ch4_mean = (moment(ch4,sdev=ch4_sdev))[0]
weff_tmp = weff[0]

;; pah demarcation
pah = wise_pahdem_ch13()
pah_wave = pah.lambda/1.E4*(1.+zz)
pah_flux = pah.flux
tmp_min = min(abs(pah_wave-weff_tmp), imin)
pah_flux = -2.5*alog10(pah_flux[imin]*pah_wave[imin]^2)+2.5*alog10(pah_flux*pah_wave^2)
;stop

;; Rieke's templates
tfile = wisenyu_path()+'/templates/rieke_all_with_maggies.fits'
riekeall = mrdfits(tfile, 1)
object = wise_riekeobject()
ii = where(strmatch(riekeall.object, object+'*', /fold_case) eq 1, nn)
ngc = riekeall[ii]

ngc_wave = ngc.lambda/1.E4*(1.+zz)
tmp_min = min(abs(ngc_wave-weff_tmp), imin)
ngc_flux = -2.5*alog10(ngc.flux[imin]*ngc_wave[imin]^2)+2.5*alog10(ngc.flux*ngc_wave^2)


;; AGN
agn = wise_agnpl(index=0.5)
agn_wave = agn.lambda/1.E4*(1.+zz)
tmp_min = min(abs(agn_wave-weff_tmp), imin)
agn_flux = -2.5*alog10(agn.flux[imin]*agn_wave[imin]^2)+2.5*alog10(agn.flux*agn_wave^2)

xthick=8
ythick=8
thick=8
charthick=3.0
charsize=1.3
psym=6
symsize=0.2

if (not keyword_set(outfile)) then $
    outfile = outdir+'wise_templates.ps'

k_print, filename=outfile, axis_char_scale=1.4, xsize=10., ysize=10.
   xtitle = 'Wavelength (\mum)'
   ytitle = '[3.4] - mag (AB)'
   plotsym, 8, 2, /fill
   djs_plot, weff, [ch1[0], ch2_mean, ch3_mean, ch4_mean], $
        psym=8, yra=[-4,5.99], xra=[1.8000, 18.0000], $
        xthick=xthick, ythick=ythick, thick=thick, $
        xtitle=xtitle, ytitle=ytitle, charsize=charsize, $
        charthick=charthick, /nodata;, /xlog

   djs_oplot, filter_lambda[*,0]/1.E4, filter_pass[*,0]/max(filter_pass[*,0])*2.0-4, $
        thick=thick, color='red'
   djs_xyouts, 2.90, $
               !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), $
               'ch1', charsize=charsize+0.5, charthick=charthick+0.6
   djs_oplot, filter_lambda[*,1]/1.E4, filter_pass[*,1]/max(filter_pass[*,1])*2.0-4, $
        thick=thick, color='red'
   djs_xyouts, 4.1, $
               !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), $
               'ch2', charsize=charsize+0.5, charthick=charthick+0.6
   djs_oplot, filter_lambda[*,2]/1.E4, filter_pass[*,2]/max(filter_pass[*,2])*2.0-4, $
        thick=thick, color='red'
   djs_xyouts, 11.0, $
               !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), $
               'ch3', charsize=charsize+0.5, charthick=charthick+0.6
   djs_oplot, filter_lambda[*,3]/1.E4, filter_pass[*,3]/max(filter_pass[*,3])*2.0-4, $
        thick=thick, color='red'
;  djs_xyouts, 7.6, $
;              !y.crange[0]+0.10*(!y.crange[1]-!y.crange[0]), $
;              'ch4', charsize=charsize+0.5, charthick=charthick+0.6
   djs_oplot, wave, flux, thick=thick+5, color='magenta'
   djs_oplot, pah_wave, pah_flux, thick=thick+5, color='brown', linestyle=2
   djs_oplot, ngc_wave, ngc_flux+0.5, thick=thick+5, color='blue', linestyle=1
   djs_oplot, agn_wave, agn_flux, thick=thick+5, color='dark green', linestyle=3

   temp = min(abs(ngc_wave-3.3*(1.+zz)), imin)
   temp = min(abs(ngc_wave-3.1*(1.+zz)), imin0)
   temp = min(abs(ngc_wave-3.0*(1.+zz)), imin1)
   djs_xyouts, ngc_wave[imin0]-0.3, ngc_flux[imin]+0.35+0.5, $
               'PAH', charsize=charsize+0.0, charthick=charthick+0.0
   djs_xyouts, ngc_wave[imin1]-0.3, ngc_flux[imin]+0.12+0.5, $
               '3.3 \mum', charsize=charsize+0.0, charthick=charthick+0.0

   temp = min(abs(ngc_wave-6.2*(1.+zz)), imin)
   temp = min(abs(ngc_wave-6.0*(1.+zz)), imin0)
   temp = min(abs(ngc_wave-5.9*(1.+zz)), imin1)
   djs_xyouts, ngc_wave[imin0]-0.3, ngc_flux[imin]+0.35+0.5, $
               'PAH', charsize=charsize+0.0, charthick=charthick+0.0
   djs_xyouts, ngc_wave[imin1]-0.3, ngc_flux[imin]+0.12+0.5, $
               '6.2 \mum', charsize=charsize+0.0, charthick=charthick+0.0

   temp = min(abs(ngc_wave-7.7*(1.+zz)), imin)
   temp = min(abs(ngc_wave-7.5*(1.+zz)), imin0)
   temp = min(abs(ngc_wave-7.2*(1.+zz)), imin1)
   djs_xyouts, ngc_wave[imin0]-0.3, ngc_flux[imin]+0.35+0.5, $
               'PAH', charsize=charsize+0.0, charthick=charthick+0.0
   djs_xyouts, ngc_wave[imin1]-0.3, ngc_flux[imin]+0.12+0.5, $
               '7.7+8.6 \mum', charsize=charsize+0.0, charthick=charthick+0.0

   items = ['IRAS12112+0305', 'Mixed', 'SSP10G04Z', textoidl('AGN (\alpha=0.5)')]
   colors = [djs_icolor('blue'), djs_icolor('brown'), djs_icolor('magenta'), djs_icolor('dark green')]
   legend, items, colors=colors, linestyle=[1,2,0,3], pspacing=1.6, $
        box=0, thick=thick+5, charsize=charsize+0.4, charthick=charthick

;  djs_xyouts, !x.crange[0]+0.07*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.75*(!y.crange[1]-!y.crange[0]), $
;              'z='+string(zz, format='(f4.2)'), charsize=charsize+0.8, charthick=charthick+0.5

   djs_oplot, weff, [ch1[0], ch2_mean, ch3_mean, ch4_mean], $
        psym=8, yra=[-4,5], xra=[1.8000, 12.0000], $
        thick=thick

k_end_print
print, ch2_mean, ch2_sdev
print, ch3_mean, ch3_sdev
print, ch4_mean, ch4_sdev

stop

end
