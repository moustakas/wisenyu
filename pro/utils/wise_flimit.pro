;+
;; [10, 10, 43, 40] \muJy SWIRE
;; [0.9, 1.7, 11.3, 14.6] \muJy SCOSMOS
;-
function wise_flimit

   mjy = [0.08, 0.11, 1.0, 6.0]

   flimit = -2.5*alog10(mjy*1.E-3/3631.)

   return, flimit
end
