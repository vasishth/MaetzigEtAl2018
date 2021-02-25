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
(search-param-space-em gg-exp1 10 '((:ga 0.5 1 0.25)))

; testing a more complicated paramspace
(search-param-space-em gg-exp1 10 '((:ga 0.25 1 0.75)
                                    (:dat 0.05 0.1 0.05)
                                    (:ans 0.15 0.45 0.15)))

; 2017-01-02: doing a realistic test run to see trace file size
(search-param-space-em gg-exp1 500 '((:ga 0.5 1.5 0.1)))

; 2017-01-03: testing fewer iterations, but more parameter combinations for sanity check
(search-param-space-em gg-exp1 100 '((:ga 0.2 1 0.2)
                                     (:dat 0.05 0.1 0.01)
                                     (:ans 0.15 0.45 0.05)))

; 2017-01-06: adding slightly wider pspace and more iterations to run on
;             different computer in the background
(search-param-space-em gg-exp1 1000 '((:ga 0.2 1 0.1)
                                     (:dat 0.05 0.1 0.01)
                                     (:ans 0.15 0.45 0.05)))

(quit)
