(in-package :cl-numlib)

(defparameter *p* #(1 2 3))
(defun count-draws-for (probs N)
  (let ((count (make-array (length probs) :element-type 'fixnum :initial-element 0))
	(sum (reduce #'+ probs)))
    (dotimes (i N)
      (incf (aref count (weighted-draw probs))))
    (values
     ;; empirical frequency
     (map 'vector #'(lambda (c) (coerce (/ c N) 'double-float)) count)
     ;; expected frequency
     (map 'vector #'(lambda (p) (coerce (/ p sum) 'double-float)) probs))))


(defun test-irs (f a b x &optional (total 10000))
  "Test inverse random sampler.  Two values returned should be about
the same."
  (let ((emp (/ (array-count 
		 (replicate (make-empirical-inverse-transform-sampler
			     f a b) total 'double-float)
		 (lambda (v) (< v x)))
		total))
	(th (funcall f x)))
    (values (coerce emp 'double-float) th)))

(test-irs #'identity 0d0 1d0 0.2d0)
(test-irs (lambda (x) (- 1 (exp (-* 3d0 x)))) 0d0 :infinity 0.2d0)
