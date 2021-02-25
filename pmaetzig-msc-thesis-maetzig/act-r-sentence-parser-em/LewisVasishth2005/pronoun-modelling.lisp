
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

;; 2017-05-26
;; test simulations for fixing the script
(search-param-space-em caplan15-sorefl 500 '((:ga 0.25 1 0.25)
                                             (:dat 0.05 0.1 0.01)
                                             (:ans 0.15 0.45 0.1)))

;; 2017-05-26
;; large scale simulations
(search-param-space-em caplan15-sorefl 1000 '((:ga 0.25 1 0.25)
                                             (:dat 0.05 0.1 0.01)
                                             (:ans 0.15 0.45 0.1)))

(quit)
