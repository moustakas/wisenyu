pro wise_redssp_maggies

indir = wisenyu_path()+'/templates/'

infile = indir+'bc03_all.fits'
if (file_test(infile) eq 0) then begin

    for imet=0,5 do begin

        case imet of
            0:  met = 0.0001/0.02
            1:  met = 0.0004/0.02
            2:  met = 0.004/0.02
            3:  met = 0.008/0.02
            4:  met = 0.02/0.02
            5:  met = 0.05/0.02
        endcase

        bc03 = k_im_read_bc03(bc03_extra=bc03_extra, metallicity=imet)
        all_age = where(bc03.age gt 1.0E9 and bc03.age lt 17.0E9, nage)
        if (n_elements(outstr) eq 0) then begin 
            strtmp = {age: ' ' $
                    , met: 0.0 $
                    , mstar: ' ' $
                    , lambda: bc03.wave $
                    , flux: fltarr(n_elements(bc03.wave)) $
                  }
            outstr = replicate(strtmp, nage*6)
        endif 

        ibegin = imet*nage
        for i=0, nage-1L do begin
            outstr[i+ibegin].age = bc03.age[all_age[i]]
            outstr[i+ibegin].met = met
            outstr[i+ibegin].mstar = bc03_extra[all_age[i]-1].m_  ;; stellar mass in this SSP
            outstr[i+ibegin].flux = bc03.flux[*,all_age[i]]
        endfor

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

outfile = indir + 'bc03_all_with_maggies.fits'
mwrfits, final, outfile, /create

end
