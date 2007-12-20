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

(defun normalize-probability (probs)
  "Normalize probs so that they sum to 1."
  (let ((sum (array-sum probs)))
    (assert (plusp sum))
    (array-map #'(lambda (p)
		   (assert (>= p 0))
		   (/ p sum)) probs)))

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
	 (result (make-ffa length 'fixnum)))
    (dotimes (i length result)
      (let ((j (weighted-draw needed)))
	(decf (aref needed j))
	(setf (aref result i) j)))))

(defun replicate (function length type)
  "Create a vector with element-type type, fill it with values
obtained from repeated evaluation of function."
  (let ((vector (make-array length :element-type type)))
    (dotimes (i length)
      (setf (aref vector i) (funcall function)))
    vector))

