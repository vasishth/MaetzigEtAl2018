;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Modelling of SR/OR, SR/OR-REFL and SR/OR-PRON sentence types from 
;;; Caplan et al. (2015).

;; Starting ACT-R, loading the model etc.
ccl
(load "../actr6/load-act-r-6.lisp")
(load "sp-lv05.lisp")

; ggf. ACT-R recompile
(push :act-r-recompile *features*)

(print-params)
(print-interface-params)

; reload
(rl)

; delete output files
(delete-output)
;(delete-trace)

; testing the changed functions -- goal: complete trace with relations printed to file
(search-param-space-em caplan15-full 10 '((:dat 0.05 0.1 0.05)))

; testing a more complicated paramspace
(search-param-space-em caplan15-full 300 '((:ga 0.25 1 0.25)
                                    (:dat 0.05 0.1 0.01)
                                    (:ans 0.15 0.45 0.1)))

(quit)
