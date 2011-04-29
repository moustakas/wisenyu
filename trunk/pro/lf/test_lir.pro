pro test_lir

atlas= read_atlas()
wise= mrdfits('atlas_wise_cat.fits',1)
wise= wise[atlas.nsaid]

wise_to_maggies, wise, mgy, ivar

m22= -2.5*alog10(mgy[3,*])

icat= where(m22 gt 8. and m22 lt 13.5)
atlas= atlas[icat]
wise= wise[icat]
mgy= mgy[*,icat]
ivar= ivar[*,icat]
m22= m22[icat]

ir= wise2lir(atlas.zdist, mgy, ivar)

end
