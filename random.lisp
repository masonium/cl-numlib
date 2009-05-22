(in-package :cl-numlib)

(defun weighted-draw (probs &optional (sum (reduce #'+ probs)))
  "Draws a nonnegative integer in [0,(length probs)) using the weights
given, where probs is a vector.  If the weights don't sum to one, they
are automatically normalized, but you can give the sum explicitly if
known."
  (let ((x (random sum)))
    (iter
      (for n :from 0)
      (for p :in-vector probs)
      (assert (<= 0 p))
      (if (< x p)
	  (return-from weighted-draw n)
	  (decf x p)))
    ;; should never get this far, only in case of numerical
    ;; inaccuracies, then we return the first index with nonzero
    ;; probability.  Mathematically, this has ever happening has zero
    ;; probability, but floating point arithmetic may be inaccurate
    (let ((first-nonzero (position-if-not #'zerop probs)))
      (unless first-nonzero
	(error "All the probabilities are zero!"))
      first-nonzero)))

(defun round-probabilities-to-count (probs n)
  "Approximate probabilities with proportional integers that add up to
n."
  ;; truncate and keep remainder
  (let* ((truncated-pairs (map 'vector (lambda (p)
					(assert (plusp p))
					(multiple-value-list (floor (* p n))))
			      (normalize-probability probs)))
	 (sum (array-sum truncated-pairs #'first))
	 (need-to-round-up (- n sum)))
    (cond
      ((zerop need-to-round-up) (array-map #'first truncated-pairs))
      ((minusp need-to-round-up)
       (error "this is not supposed to happen, check algorithm"))
      (t 
       ;; append an index so that we can keep track of permutations
       (dotimes (i (length truncated-pairs))
	 (setf (aref truncated-pairs i) (cons i (aref truncated-pairs i))))
       ;; sort based on remainder
       (setf truncated-pairs (stable-sort truncated-pairs #'< :key #'third))
       (let ((number-of-probs (length probs)))
	 ;; increase last need-to-round-up counts
	 (iter
	   (for i :from (- number-of-probs need-to-round-up)
		:below number-of-probs)
	   (incf (second (aref truncated-pairs i)))))
       ;; sort and collect
       (array-map #'second (sort truncated-pairs #'< :key #'first))))))

(defun fill-with-weighted-draw (count)
  "Create random a vector of fixnums that contains (aref count i) of
i."
  (let* ((needed (array-copy count))
	 (length (array-sum count))
	 (result (make-array length :element-type 'fixnum)))
    (dotimes (i length result)
      (let ((j (weighted-draw needed)))
	(decf (aref needed j))
	(setf (aref result i) j)))))

(defun replicate (function dimensions type)
  "Create an array with element-type type and given dimensions, fill
it with values obtained from repeated evaluation of function."
  (let ((array (make-array dimensions :element-type type)))
    (dotimes (i (array-total-size array))
      (setf (row-major-aref array i) (funcall function)))
    array))

(defun make-empirical-inverse-transform-sampler (F a b &key 
						 (F-maximum 1d0))
  "Return a function which will return random draws from the
probability distribution F, which has support [a,b], where `a' can
be :minusinfinity, `b' can be :inifinity.  Optionally, you can specify
F(b) as F-maximum (so you don't have to normalize functions), F(0) has
to be 0, and F has to be strictly monotone increasing (the latter two
conditions are of course not verified).

Notes: uses inverse transform sampling.  Uses rootfinding, should only
used when the inverse of F is not readily available."
  (assert (and (or (realp a) (eq a :minusinfinity))
	       (or (realp b) (eq b :infinity))))
  (let ((positive-rule (if (realp a)
			   (make-contracting-rule a 0.5d0) ; contact to a
			   (make-expanding-rule -1d0 2d0))) ; expand to -Inf
	(negative-rule (if (realp b)
			   (make-contracting-rule b 0.5d0) ; contact to b
			   (make-expanding-rule 1d0 2d0))) ; expand to Inf
	(x0 (cond
	      ((and (realp a) (realp b)) 0d0) ; not used
	      ((realp b) (- b 1d0))
	      ((realp a) (+ a 1d0))
	      (t 0d0))))
    (if (and (realp a) (realp b))
	(lambda ()
	  (let ((y (random F-maximum)))
	    (flet ((f (x) (- (funcall f x) y)))
	      (root-ridders #'f a b))))
	(lambda ()
	  (let ((y (random F-maximum)))
	    (flet ((f (x) (- (funcall f x) y)))
	      (root-autobracket #'f x0 negative-rule positive-rule)))))))

;; random draws
(defun exponential-random (rate &optional (random-state *random-state*))
  "Generate a random draw from an exponential distribution with given rate."
  (let ((u (random 1d0 random-state)))
    (/ (- (log (- 1d0 u))) rate)))

;; (defparameter *foo* (make-array 10000 :element-type 'double-float))
;; (dotimes (i 10000)
;;   (setf (aref *foo* i) (exponential-random 2d0)))
;; (array-mean *foo*)

(defun first-exponential-random (rates &key (random-state *random-state*)
				 (sum (array-sum rates)))
  "Imagine n independent Poisson processes with the given rates.  The
event is defined as the one occuring first, the function will return
the time to the event and the index of the event as (values time
index).  `rates' needs to be a vector."
  (values (exponential-random sum random-state)
          (weighted-draw rates sum)))

(defun make-poisson-random (lambda &optional (epsilon 1d-10))
  "Return a function which takes no arguments and returns a Poisson
random variate with mean lambda.  The table is created up to the point
where the probability of drawing that element is smaller than
epsilon."
  (assert (plusp epsilon))
  ;; fill table
  (let* ((sum 0d0)
	 (lambda (coerce lambda 'double-float))
	 (p (exp (- lambda)))
	 (table (iter
		  (for i :from 0)
		  (when (< p epsilon)
		    (finish))
		  (incf sum p)
		  (collect sum :result-type '(simple-array double-float (*)))
		  (multf p (/ lambda (1+ i)))))
	 (n (1- (length table))))
    ;; function that generates random variate
    (lambda ()
      (let ((x (random sum)))
	(if (< x (aref table 0))
	    0
	    (interval-binary-search table x 0 n))))))

;; (coerce (array-mean (replicate (make-poisson-random 1d0) 100000 'fixnum))
;; 		    'double-float)
    
