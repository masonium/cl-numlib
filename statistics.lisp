(in-package :cl-numlib)

(defun quantiles (data quantiles)
  "Calculate and return the empirical quantiles of data, both of which
are vectors.  Quantiles are interpolated between data points.

The data is copied and then sorted, so it is wiser to call this
function once with many quantiles on the same data."
  (assert (and (vectorp data) (vectorp quantiles)))
  (let ((copy (array-copy data))
	(n (length data)))
    (sort copy #'<)
    (flet ((quantile (q)
	     (assert (<= 0 q 1))
	     (bind ((position (* q (1- n)))
		    ((:values int-part frac-part) (floor position)))
	       (if (zerop frac-part)
		   (aref copy int-part)
		   (convex-combination (aref copy (1+ int-part))
				       (aref copy int-part) frac-part)))))
      (array-map #'quantile quantiles))))

(defun quantiles-weighted (data quantiles key weight-key)
  "Calculate and return the empirical quantiles of a weighted data
series, where functions `key' and `weight' are used to extract the
value and the weight from elements of `data'.  Return a single vector
of double-floats.

Quantiles are interpolated between data points.  The data is copied
and then sorted, so it is wiser to call this function once with many
quantiles on the same data."
  (assert (and (vectorp data) (vectorp quantiles)))
  (let ((copy (array-copy data))
	(n (length data)))
    (sort copy #'< :key key)
    (let ((cumsum-weights (make-array n :element-type 'double-float))
	  (extracted-values (map '(simple-array * (*)) key copy))
	  (sum 0d0))
      ;; cumulative sum
      (iter
	(for elt :in-vector copy)
	(for i :from 0)
	(for w := (funcall weight-key elt))
	(assert (plusp w))
	(setf (aref cumsum-weights i) sum)
	(incf sum w))
      ;; interpolation of the inverse
      (let ((last-weight (aref cumsum-weights (1- (length cumsum-weights))))
	    (interp (car (make-interpolating-functions cumsum-weights 
						       extracted-values))))
	;; calculating quantiles
	(array-map (lambda (x)
;		     (assert (<= 0d0 x 1d0))
		     (coerce (funcall interp (* x last-weight)) 'double-float))
		   quantiles 'double-float)))))

;; (defparameter *data* #((1 1) (0 2) (3 1)))

;; (quantiles-weighted *data* (num-sequence :from 0 :to 1 :by 1/5 :type 'rational) #'first #'second)

    
;; (defun mean-variance-weighted-old (data key weight-key)
;;   "Calculate and return the mean and variance of a weighted data vector,
;; where functions `key' and `weight' are used to extract the value and
;; the weight from elements of `data'.  Return (values mean variance),
;; both are double-float."
;;   (let ((weighted-sum 0d0)
;; 	(weighted-square-sum 0d0)
;; 	(cumulative-weights 0d0))
;;     (dotimes (i (length data))
;;       (let* ((elt (aref data i))
;; 	     (value (funcall key elt))
;; 	     (weight (funcall weight-key elt)))
;; 	(incf weighted-sum (* value weight))
;; 	(incf weighted-square-sum (* value (square weight)))
;; 	(incf cumulative-weights weight)))
;;     (let ((mean (/ weighted-sum cumulative-weights))
;; 	  (mean-square (/ weighted-square-sum cumulative-weights)))
;;       (values mean (- mean-square (square mean))))))

;; (declaim (optimize (speed 3)))

(defun mean-variance-weighted (function n)
  "Calculate and return the mean and variance of a weighted data vector.
 (function i) should return a vector of 2 doubles, a value and its
 weight, for 0 <= i < n, both are double-float.  Return (values mean
 variance), both are double-float.  i is a fixnum."
  (let ((weighted-sum 0d0)
	(weighted-square-sum 0d0)
	(cumulative-weights 0d0))
    (declare ; ((function (fixnum) (simple-array double-float (2))) function)
	     (fixnum n)
	     (double-float weighted-sum weighted-square-sum
			   cumulative-weights))
    (dotimes (i n)
      (declare (fixnum i))
;;;       (let* ((value-weight (funcall function i))
;;; 	     (value (aref value-weight 0))
;;; 	     (weight (aref value-weight 1)))
      (bind ((#(value weight) (funcall function i)))
	(declare ;((simple-array double-float (2)) value-weight)
		 (double-float value weight))
	(incf weighted-sum (* value weight))
	(incf weighted-square-sum (* value value weight))
	(incf cumulative-weights weight)))
    (if (plusp cumulative-weights)
	(let ((mean (/ weighted-sum cumulative-weights))
	      (mean-square (/ weighted-square-sum cumulative-weights)))
	  (values mean (- mean-square (* mean mean))))
	(values nil nil))))
 
