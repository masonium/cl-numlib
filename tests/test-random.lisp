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
