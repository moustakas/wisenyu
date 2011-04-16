;+
; NAME:
;   wise_agnpl
; PURPOSE:
; CALLING SEQUENCE:
; OPTIONAL INPUTS:
; OUTPUTS:
; COMMENTS:
; REVISION HISTORY:
;   14-Apr-2011: Started, Guangtun Zhu, NYU
;-

function wise_agnpl, index=index

   if (n_elements(index) eq 0) then index=0.

   wave = findgen(20000L)*10.+1000.
;; index: f(v) = v^(-index), f(lambda) = lambda^(index-2)
   flux = (wave/10000.)^(index-2.)

   strtmp = {lambda: wave, flux:flux}

   return, strtmp
end
