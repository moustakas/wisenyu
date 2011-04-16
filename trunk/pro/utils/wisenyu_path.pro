;+
; NAME:
;   wisenyu_path
; PURPOSE:
;   given a type of data, return the path (directory).
; CALLING SEQUENCE:
; OPTIONAL INPUTS:
; OUTPUTS:
; COMMENTS:
; REVISION HISTORY:
;   14-Apr-2011: Started, Guangtun Zhu, NYU
;-

function wisenyu_path, lowz=lowz, dimage=dimage

   path = '/mount/hercules5/sdss/atlas/v0/wise-cat/'

   if (keyword_set(lowz)) then path = '/global/data/lowz-sdss/'
   if (keyword_set(dimage)) then path = getenv('DIMAGE_DIR')+'/data/atlas/'

   return, path
end
