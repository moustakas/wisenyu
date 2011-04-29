;+
; NAME: 
;   wise_22mu_lf
; PURPOSE: 
;   Calculate the 22 micron LF
; CALLING SEQUENCE: 
;   wise_22mu_lf
; REVISION HISTORY:
;   14-Apr-2011 MRB NYU
;-
pro wise_22mu_lf

m22lim= [8., 13.5]
mrlim= [10., 17.]
zlim= [0.01, 0.05]

absmag=mrdfits('atlas_wise_absmag.fits',1)
atlas= read_atlas()
absmag= absmag[atlas.nsaid]

ilss= get_sdss_ilss(atlas.ra, atlas.dec, sample='dr72', icomb=icomb)
lss= mrdfits(getenv('LSS_REDUX')+'/dr72/lss_geometry.dr72.fits',1)

euler, atlas.ra, atlas.dec, elon, elat, 3

indx= where(absmag.mag[2] gt mrlim[0] AND $
            absmag.mag[2] lt mrlim[1] AND $
            absmag.mag[8] gt m22lim[0] AND $
            absmag.mag[8] lt m22lim[1] AND $
            atlas.zdist gt zlim[0] AND $
            atlas.zdist lt zlim[1] AND $
            ((elon gt 30.5 and elon lt 130.5) OR $
             (elon gt 204.5 and elon lt 270.)) AND $
            ilss ge 0 and lss[ilss].fgotmain gt 0.5 and $
            lss[ilss].mmax ne 0, count)

vmax= lf_simple_vmax(absmag[indx].absmag[8], m22lim, zlim)
vmax.vmax= vmax.vmax*(3210.*!DPI^2/180.^2)/(4.*!DPI)

mwrfits, atlas[indx], 'atlas_wise_22mu_sample.fits', /create
mwrfits, absmag[indx], 'atlas_wise_22mu_sample.fits'
mwrfits, vmax, 'atlas_wise_22mu_sample.fits'


k_print,filename='lf_wise_22mu.ps'

hogg_plothist,absmag[indx].absmag[8], weight=1./vmax.vmax, npix=50, $
  xra=[-19., -26.], xvec=xvec, xtitle=textoidl('M_{22}'), $
  ytitle=textoidl('\Phi(M_{22})'), hist=phi, err=phi_err, $ 
  /ylog, yra=[ 1.e-6, 2.e-2]

k_end_print

END 
