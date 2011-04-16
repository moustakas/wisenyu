pro wise_rieke_maggies

indir = wisenyu_path()+'/templates/'

infile = indir+'rieke_all.fits'

if (file_test(infile) eq 0) then begin

    in_obj = mrdfits(indir+'/apj295357t3_mrt.fits', 1)
    nobj = n_elements(in_obj)
    strtmp = {object: ' ' $
            , lambda_ori: fltarr(870) $
            , flux_ori: fltarr(870) $
            , lambda: fltarr(8691) $
            , flux: fltarr(8691) $
            }
    outstr = replicate(strtmp, nobj)
    outstr.object = in_obj.name
    outstr.lambda_ori = in_obj.wave
    outstr.flux_ori = in_obj.flux

    for i=0L, nobj-1L do begin
        temp_lambda = in_obj[i].wave*1.E4
        temp_flux = in_obj[i].flux/in_obj[i].wave/in_obj[i].wave
        ii = where(temp_lambda gt 4500. and temp_lambda lt 5500., nn)
        temp_flux = temp_flux/mean(temp_flux[ii])*1.E-17
        ;; resampling:
        temp_lambda1 = (rebin(temp_lambda, 8700))[0:8690]
        temp_flux1 = interpol(temp_flux, temp_lambda, temp_lambda1)
        outstr[i].lambda = temp_lambda1
        outstr[i].flux = temp_flux1
    endfor

    mwrfits, outstr, infile, /create
endif

if (n_elements(outstr) eq 0) then $
    outstr = mrdfits(infile, 1)
nobj = n_elements(outstr)

outfilters = [primus_sdss_filterlist(), $
              primus_twomass_filterlist(), $
              primus_irac_filterlist(), $
              wise_filterlist()]
nfilter = n_elements(outfilters)

zbin = 0.025
zmin = 0.00
zmax = 0.50
nz = (zmax-zmin)/zbin + 1.
nz = fix(nz)
zvals = findgen(nz)*zbin+zmin

allstr = {filters: outfilters, $
          zvals: fltarr(nz), $
          maggies: fltarr(nz, nfilter)}
allstr = replicate(allstr, nobj)

for i=0L, nobj-1L do begin
    olambda = k_lambda_to_edges(outstr[i].lambda)
    oflux = outstr[i].flux
    k_projection_table, maggies, oflux, olambda, zvals, outfilters
    allstr[i].maggies = maggies
    allstr[i].zvals = zvals
endfor

final = struct_addtags(outstr, allstr)

outfile = indir + 'rieke_all_with_maggies.fits'
mwrfits, final, outfile, /create

end
