ccl
(load "../actr6/load-act-r-6.lisp")
(load "sp-lv05.lisp")

; ggf. ACT-R recompile
(push :act-r-recompile *features*)

(print-params)

; reload
(rl)

; delete output files
(delete-output)
(delete-trace)

;; GETTING PARSER PREDICTIONS
;;; two values of :ga
;;; Date: 2016-10-12
;;; :mas is set to 3.5 according to Felix' diss
;;; :mp is not yet set to NIL, try this in another step
(setf i 1)
(sgp :mas)

;; Redoing the simulations, 2016-12-07 --  2016-12-12
;; 
;; This is, firstly, to check wether taking the mean of the four simulations
;; of the last time is valid and yields a similar result as just doing 500
;; simulations.
(search-param-space-em gg-exp1 500 '((:ga 0.5 1.5 0.25)))

(loop for i from 1 to 4
      do (search-param-space-em gg-exp1 100 '((:ga 0.5 1.5 0.25))) )
(rl)
(loop for i from 1 to 4
      do (search-param-space-em gg-exp1 100 '((:dat 0.05 0.1 0.05))) )
(rl)
(loop for i from 1 to 4
      do (search-param-space-em gg-exp1 100 '((:ga 0.5 1.5 0.25)
                                              (:dat 0.05 0.1 0.025))) )
(rl)
(loop for i from 1 to 4
      do (search-param-space-em gg-exp1 100 '((:ga 0.5 1.5 0.25)
                                              (:ans 0.15 0.25 0.02))) )
(rl)
(loop for i from 1 to 4
      do (search-param-space-em gg-exp1 100 '((:dat 0.05 0.1 0.025)
                                              (:ans 0.15 0.25 0.02))) )

(loop for i from 1 to 4
      do (search-param-space-em gg-exp1 100 '((:ga 0.2 0.8 0.2)
                                              (:dat 0.075 0.1 0.025)
                                              (:ans 0.16 0.2 0.01))) )

;; from 2016-12-16 on
;; this will take at least a full day to run, should do it on another machine
(search-param-space-em gg-exp1 1000 '((:ga 0.2 1.5 0.1)
                                      (:ans 0.15 0.5 0.01)))

;; trying to replicate the behaviour of ga deviating from 1.0
(search-param-space-em gg-exp1 200 '((:ga 0.5 1.5 0.25)
                                     (:ans 0.15 0.25 0.02)))

;;; GETTING ATTACHMENT DATA
; move standard-output to file in ./output/
(setf *log* (open "/Users/paul/ownCloud/potsdam_msc_ling/msc-thesis-maetzig/act-r-sentence-parser-em/LewisVasishth2005/output/trace.txt" :direction :output :if-exists :append :if-does-not-exist :create))
(setf *standard-output* *log*)
(force-output *log*)

(ps "the dog bit the boy *")

(quit)
