;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Modelling response accuracy data via attachments in ACT-R

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
(delete-trace)

; testing the changed functions -- goal: complete trace with relations printed to file
; run several times and compare the results
(search-param-space-em gg-exp1 1000 '((:ans 0.15 0.50 0.35)))

(quit)
