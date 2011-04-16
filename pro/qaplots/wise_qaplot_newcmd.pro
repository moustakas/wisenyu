;; to-do
;; AGN fraction as a function of color
;; Define green valley consistently
;; Define transitioning objects consistently
;; 15-Apr-2011, Guangtun Zhu, Started, NYU
pro wise_qaplot_newcmd, outfile=outfile

if (not keyword_set(outdir)) then $
    outdir = wisenyu_path()+'/qaplots/'
indir = wisenyu_path()
sdss_indir = '~/mrb_dimage/data/atlas/'

zpoint = wise_vega2ab()
flimit = wise_flimit()

;; readin the sample
wise0 = mrdfits(indir+'atlas_wise_cat.fits', 1)
wise_kcorrect0 = mrdfits(indir+'atlas_wise_absmag.fits', 1)
atlas0 = mrdfits(sdss_indir+'atlas.fits',1)

ii = where(atlas0.isdss ge 0 and wise0.w1mag_8 lt flimit[0] and wise0.w3mag_8 lt flimit[2] $
and wise0.w1mag_8 gt 0 and wise0.w3mag_8 gt 0, nn)
atlas = atlas0[ii]
wise = wise0[ii]
wise_kcorrect = wise_kcorrect0[ii]

sdss0 = mrdfits(sdss_indir+'sdss_atlas.fits',1)
sdss_line0 = mrdfits(sdss_indir+'sdssline_atlas.fits',1)
sdss = sdss0[atlas.isdss]
sdss_line = sdss_line0[atlas.isdss]
sdss_to_maggies, maggies, ivar, calibobj=sdss
mr = -2.5*alog10(maggies[2, *])

mlimit = 2.5*alog10(0.7) - lf_distmod(0.05) + 17.77 - 0.2
;jj = where(mr lt 18.0 and wise_kcorrect.absmag[2] lt mlimit)
jj = where(mr lt 17.77)
atlas = atlas[jj]
wise = wise[jj]
wise_kcorrect = wise_kcorrect[jj]
sdss = sdss[jj]
sdss_line = sdss_line[jj]

;; temporary
;wise_rsbcgv, wise_kcorrect, ired=ired, iblue=iblue, igreen=igreen
;; 2.6-0.074(x+20), red
;; 1.5-0.074(x+19), blue
;; 2.3-0.074(x+20) and 1.8-0.08(x+19), green
mag = wise_kcorrect.absmag[2]
color = wise_kcorrect.absmag[0] - wise_kcorrect.absmag[2]
ircolor = wise.w1mag_8 + zpoint[0] - wise.w3mag_8 - zpoint[2]
redc = 2.45
bluec = 1.90
ired = where(color gt redc-0.08*(mag+20))
iblue = where(color lt bluec-0.08*(mag+19))
igreen = where(color ge bluec-0.08*(mag+19) and color le redc-0.08*(mag+20))

inoblue = [ired, igreen]
xagn = alog10(sdss_line0.n2lum/sdss_line0.halum) 
yagn = alog10(sdss_line0.o3lum/sdss_line0.hblum)
iagn = where(xagn[inoblue] gt alog10(0.6) and yagn[inoblue] gt alog10(3.0))


;; all the objects with redshifts of good quality
;infile = indir+version+'.fits.gz'
;esam = mrdfits(infile, 1)
;esam = primus_sfrs_select(esam)


;; Apply redshift cut, and separate samples to swire and scosmos subsamples
;;esam = primus_sfrs_select(esam, zmin=0.20, zmax=0.40, swire=swire, scosmos=scosmos)

;; Apply flux limit cut in i band and ch1
;primus_sfrs_flimit_cut, swire, scosmos, out_swire=swire_tmp, out_scosmos=scosmos_tmp
;swire = swire_tmp
;scosmos = scosmos_tmp

;; Separate sample into detections and non-detections in ch4
;; Note this changes the flux of ch4 for non-detections
;primus_sfrs_ch4sep, swire, scosmos, $
;           out_swire=swire_tmp, out_scosmos=scosmos_tmp, $
;           swire_noch4=swire_noch4, scosmos_noch4=scosmos_noch4

;; want both
;sample = [swire_tmp, swire_noch4, scosmos_tmp, scosmos_noch4]
;
;primus_sfrs_rsbcgv, sample, ired=ired, iblue=iblue, igreen=igreen
;sample_green = sample[igreen]
;sample_blue = sample[iblue]
;sample_red = sample[ired]

;ch1 = -2.5*alog10(sample.maggies_ap2[0])
;ch4 = -2.5*alog10(sample.maggies_ap2[3])
;inz = sample.zprimus
;inch14 = ch1-ch4
;outch14 = inch14
;outz = 0.3
;for i=0, n_elements(inz)-1L do begin
;    outch14[i] = primus_sfrs_simplekcorrect_ch14(ch1[i], ch4[i], inz[i], outz)
;endfor

;coeffs = primus_sfrs_pahdem_ch14()
pah = wise_pahdem_ch13(fracsf=0.02)
pahcolor = pah.color[0]

xthick=8
ythick=8
thick=8
charthick=3
charsize=1.5
psym=6 
symsize=0.7

if (not keyword_set(outfile)) then $
    outfile = outdir+'wise_newcmd.ps'

;psfile1 = repstr(outfile, '.ps', '_pg1.ps')
psfile1 = outfile
k_print, filename=psfile1, axis_char_scale=1.3, xsize=15., ysize=10.
    !p.multi=[0,2,2]
    !x.margin=0
    !y.margin=0
    xtitle = textoidl('!8M_{r}-5!6log_{10}!8h!6')
    ytitle = textoidl('!8M_{u} - M_{r}!6')
    xra = [-14.5, -23.6]
    yra = [0.51, 3.49]

    plotsym, 0, 1, /fill
    flevels = [0.40, 0.60, 0.95]
    hogg_scatterplot, mag, color, xthick=xthick, ythick=ythick,thick=thick, $
              xra=xra, yra=yra, levels=flevels, exponent=0.5, $
              xnpix=40L, ynpix=40L, charthick=charthick, charsize=charsize, $
              xtickformat='(A1)', ytitle=ytitle, title=title, outliers=outliers
    djs_oplot, mag[inoblue[iagn]], color[inoblue[iagn]], psym=6, symsize=0.2, color='green', thick=thick
    legend = 'Optical color'
    djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.8, charthick=5., color='blue'
    legend = 'Red sequence'
    djs_xyouts, !x.crange[0]+0.69*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.70*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.7, charthick=6., color='red'
;   legend = 'DUSTY SF / TRANSITION'
    legend = 'Green valley'
;   djs_xyouts, !x.crange[0]+0.45*(!x.crange[1]-!x.crange[0]), $
    djs_xyouts, !x.crange[0]+0.69*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.53*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.7, charthick=6., color='dark green'
    legend = 'Blue cloud'
;   djs_xyouts, !x.crange[0]+0.32*(!x.crange[1]-!x.crange[0]), $
    djs_xyouts, !x.crange[0]+0.69*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.30*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.7, charthick=6., color='blue'
    djs_oplot, !x.crange, bluec-0.074*(!x.crange+19), thick=thick
    djs_oplot, !x.crange, redc-0.074*(!x.crange+20), thick=thick


;   djs_plot, mag[ired], color[ired], xthick=xthick, ythick=ythick,thick=thick, $
    djs_plot, mag[igreen], color[igreen], xthick=xthick, ythick=ythick,thick=thick, $
              xra=xra, yra=yra, charthick=charthick, charsize=charsize, $
              xtickformat='(A1)', ytickformat='(A1)', psym=8, symsize=0.5, $
              color='dark green'
    djs_oplot, !x.crange, bluec-0.074*(!x.crange+19), thick=thick
    djs_oplot, !x.crange, redc-0.074*(!x.crange+20), thick=thick

;   legend = 'Red sequence only'
    legend = 'Green valley only'
    djs_xyouts, !x.crange[0]+0.50*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.10*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.8, charthick=5., color='black'

;   ytitle = textoidl('!8!u0.3!n([8.0] - [3.6])!6')
    ytitle = textoidl('!8M_{[12.0]} - M_{[3.4]}!6')
    yra = [-5.49+2.5, 0.99+2.5]
    hogg_scatterplot, mag, -ircolor, xthick=xthick, ythick=ythick,thick=thick, $
              xra=xra, yra=yra, levels=flevels, exponent=0.5, $
              xnpix=40L, ynpix=40L, charthick=charthick, charsize=charsize, $
              xtitle=xtitle, ytitle=ytitle, title=title, outliers=outliers
    djs_oplot, mag[inoblue[iagn]], -ircolor[inoblue[iagn]], psym=6, symsize=0.2, color='green', thick=thick

    legend = 'IR color'
    djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.8, charthick=5., color='blue'

    legend = 'Quiescent'
    djs_xyouts, !x.crange[0]+0.69*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.70*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.7, charthick=6., color='red'
;   legend = 'DUSTY SF / TRANSITION'
    legend = 'Transition'
;   djs_xyouts, !x.crange[0]+0.45*(!x.crange[1]-!x.crange[0]), $
    djs_xyouts, !x.crange[0]+0.69*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.50*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.7, charthick=6., color='dark green'
    legend = 'Star-forming'
;   djs_xyouts, !x.crange[0]+0.32*(!x.crange[1]-!x.crange[0]), $
    djs_xyouts, !x.crange[0]+0.69*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.30*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.7, charthick=6., color='blue'
    legend = 'SSFR'
    djs_xyouts, !x.crange[0]+0.02*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.9, charthick=5., color='blue'
    tmp_x = !x.crange[0]+0.072*(!x.crange[1]-!x.crange[0])
    tmp_y = !y.crange[0]+0.78*(!y.crange[1]-!y.crange[0])
    one_arrow, tmp_x, tmp_y, 90, '', color=djs_icolor('blue'), $
        charsize=2., arrowsize=[240.0, 50.0, 10.0], thick=10., /data
;   djs_oplot, !x.crange, [-pahcolor, -pahcolor], thick=thick


;   djs_plot, mag[ired], -ircolor[ired], xthick=xthick, ythick=ythick,thick=thick, $
    djs_plot, mag[igreen], -ircolor[igreen], xthick=xthick, ythick=ythick,thick=thick, $
              xra=xra, yra=yra, charthick=charthick, charsize=charsize, $
              ytickformat='(A1)', xtitle=xtitle, psym=8, symsize=0.5, $
              color='dark green'
;   legend = 'Red sequence only'
    legend = 'Green valley only'
    djs_xyouts, !x.crange[0]+0.50*(!x.crange[1]-!x.crange[0]), $
                !y.crange[0]+0.10*(!y.crange[1]-!y.crange[0]), $
                legend, charsize=1.8, charthick=5., color='black'


k_end_print

spawn, 'convert '+psfile1+' '+repstr(psfile1, '.ps', '.png')
stop

end
