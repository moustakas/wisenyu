;+
; NAME:
;  wise_filterlist
; PURPOSE:
;   return list of WISE broadband filter names
; CALLING SEQUENCE:
;   filterlist=wise_filterlist()
; OUTPUTS:
;   filterlist - [4] list of filter names
; COMMENTS:
;   See the documentation of wise_to_maggies.pro
; REVISION HISTORY:
;   14-Apr-2011  Guangtun Zhu, NYU
;-
;------------------------------------------------------------------------------
function wise_filterlist, bandnames=bandnames

    bandnames = ['w1',$
      'w2',$
      'w3',$
      'w4']
    filterlist='wise_'+bandnames+'.par'

return, filterlist
end
