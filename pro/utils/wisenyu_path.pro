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

function wise_path, lowz=lowz, atlas=atlas

   path = '/mount/hercules5/sdss/atlas/v0/wise-cat/'

   if (keyword_set(lowz)) then path = '/global/data/lowz-sdss/'
;  if (keyword_set(atlas)) then path = '/global/data/lowz-sdss'

   return, path
end
