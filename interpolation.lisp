(in-package :cl-numlib)

(defun make-interpolating-functions (x &rest ys)
  "Return a list of univariate interpolating functions.

For each y in ys, x -> y is interpolated linearly.  x needs to be
strictly increasing, and all the ys need to have the same length as
x.

For values outside x, the value at the nearest edge is returned (ie
constant extrapolation)."
  (let ((binning-function (make-binning-function x))
	(n (length x)))
    (iter
      (for y :in ys)
      (assert (and (vectorp y) (= (length y) n)))
      (let ((y (array-copy y)))
	(collect (lambda (v)
		   (let ((bin (funcall binning-function v)))
		     (cond
		       ((<= bin 0) (aref y 0))
		       ((<= n bin) (aref y (1- n)))
		       (t (let ((x-left (aref x (1- bin)))
				(x-right (aref x bin))
				(y-left (aref y (1- bin)))
				(y-right (aref y bin)))
			    (+ y-left (* (- v x-left)
					 (/ (- y-right y-left)
					    (- x-right x-left))))))))))))))
 
(defun make-interpolating-function (x y)
  "Makes a single interpolating function.  See
make-interpolating-functions for details."
  (first (make-interpolating-functions x y)))
