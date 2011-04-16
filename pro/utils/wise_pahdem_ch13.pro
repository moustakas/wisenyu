function wise_pahdem_ch13, fracsf=fracsf

if (n_elements(fracsf) eq 0) then fracsf=0.10
indir = wisenyu_path()+'/templates/'
infile = indir+'wise_pahdem_template_with_maggies.fits'
allstr = mrdfits(infile, 1)

tmp_min = min(abs(allstr.fracsf-fracsf), ii)

outstr = {zvals: allstr[ii].zvals, $
          color: -2.5*alog10(allstr[ii].maggies[*, 12])+2.5*alog10(allstr[ii].maggies[*, 14]), $
          lambda: allstr[ii].lambda, $
          flux: allstr[ii].flux $
         }

return, outstr
end
