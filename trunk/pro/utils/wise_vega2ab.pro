;;http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html#WISEZMA
function wise_vega2ab
;; m(AB) = m(vega)+zpoint
;; 3.6, 4.5, 5.8, 8.0
    zpoint = [2.683, 3.319, 5.242, 6.604]
    return, zpoint
end
