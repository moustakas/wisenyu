pro wise_pahdem_template 

indir = wisenyu_path()+'/templates/'

file = indir+'rieke_all_with_maggies.fits'
all = mrdfits(file, 1)
;; the object we like
object = wise_riekeobject()
ii = where(strmatch(all.object, object+'*', /fold_case) eq 1, nn)
ngc = all[ii]

sspfile = indir+'bc03_all_with_maggies.fits'
ssp = mrdfits(sspfile, 1)
tmp_min = min(abs(ssp.age-10.0E9)+abs(ssp.met-0.4), imin)
;issp = where(ssp.age eq 10.0E9 and ssp.met eq 0.4, nssp)
redssp = ssp[imin]
red_wave = redssp.lambda
red_flux = redssp.flux


tmp_wave = findgen(1000L)*100.+10E+4
tmp_flux = interpol(red_flux, red_wave, tmp_wave)

ii = where(red_wave lt 10E+4)
red_wave = [red_wave[ii], tmp_wave]
red_flux = [red_flux[ii], tmp_flux]

ngc_flux = interpol(ngc.flux, ngc.lambda, red_wave)

;normalize
tmp_min = min(abs(red_wave - 2E4), imin) 
red_flux = red_flux/red_flux[imin];*red_wave/red_wave[imin]
ngc_flux = ngc_flux/ngc_flux[imin];*red_wave/red_wave[imin]

djs_plot, red_wave, red_flux, xra=[2., 10.]*1E4, yra=[0,2.8]
djs_oplot, red_wave, ngc_flux
stop

nbin = 101
fracsf = findgen(nbin)/float(nbin-1)

outfilters = [primus_sdss_filterlist(), $
              primus_twomass_filterlist(), $
              primus_irac_filterlist(), $
              wise_filterlist()]
nfilter = n_elements(outfilters)

zbin = 0.025
zmin = 0.00
zmax = 1.20
nz = (zmax-zmin)/zbin + 1.
nz = fix(nz)
zvals = findgen(nz)*zbin+zmin

strtmp = {fracsf: 0.0, $
          lambda: red_wave, $
          flux: red_flux, $
          filters: outfilters, $
          zvals: fltarr(nz), $
          maggies: fltarr(nz, nfilter)}
allstr = replicate(strtmp, nbin)

olambda = k_lambda_to_edges(red_wave)
for ibin=0L, nbin-1L do begin
    oflux = red_flux*(1.-fracsf[ibin]) + ngc_flux*fracsf[ibin]
    k_projection_table, maggies, oflux, olambda, zvals, outfilters
    allstr[ibin].maggies = maggies
    allstr[ibin].zvals = zvals
    allstr[ibin].fracsf = fracsf[ibin]
    allstr[ibin].flux = oflux
endfor

outfile = indir+'wise_pahdem_template_with_maggies.fits'
mwrfits, allstr, outfile, /create
stop

x = allstr[0].zvals
i=0 & y = -2.5*alog10(allstr[i].maggies[*, 12])+2.5*alog10(allstr[i].maggies[*, 14]) & djs_plot, x, y, xra=[0., 0.7], yra=[-3,3] & print, allstr[i].fracsf
i++ & y = -2.5*alog10(allstr[i].maggies[*, 12])+2.5*alog10(allstr[i].maggies[*, 14]) & djs_oplot, x, y & print, allstr[i].fracsf

end
