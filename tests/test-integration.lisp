(in-package :cl-numlib)




(trapezoidal-closed (map-a-b #'identity 0 10) 1e-4)

;; exponential distribution, should integrate to 1
(let ((l 0.5))
  (trapezoidal-closed (map-a-inf (lambda (x) (* l (exp (- (* l x))))) 0)
		      1e-7))

;; shifted exponential distribution
(let ((l 0.2)
      (a -9))
  (trapezoidal-closed (map-a-inf (lambda (x) (* l (exp (- (* l (- x a)))))) a)
		      1e-5))

(let ((l 0.2)
      (a -9))
  (integrate (lambda (x) (* l (exp (- (* l (- x a))))))
	     a 'inf))

;; exponential distribution on -x
(let ((l 0.2)
      (b -99))
  (integrate (lambda (x) (* l (exp (* l (- x b)))))
	     'inf b))

(funcall (map-a-inf (lambda (x) (* 0.1 (exp (- (* 0.1 x))))) 0) .1)
