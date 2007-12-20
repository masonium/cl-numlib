(in-package :cl-numlib)

(defun int-sequence (a &optional b)
  "If b is given, return a vector with integers from a to
b (inclusive), otherwise from 1 to a (inclusive)."
  (unless b
    (setf b a
	  a 1))
  (assert (<= a b))
  (let ((result (make-array (1+ (- b a)) :element-type 'fixnum)))
    (iter
      (for i :from 0)
      (for v :from a :to b)
      (setf (aref result i) v))
    result))

(defun num-sequence (&key (from nil from-given-p) (to nil to-given-p)
		     (by nil by-given-p) (length nil length-given-p)
		     (type 'double-float))
  "Create a one-dimensional array which contains a sequence with the
given parameters.  The parameters are checked for consistency.  Two
values are returned, first the sequence (simple-array with type), then
the fourth parameter which was not given.  Note: the sequence stops at
or before reaching \"to\", eg (num-sequence :from 0 :to 10 :by pi)
will stop end at (* 3 pi)."
  (unless (= 3 (reduce #'+
		       (mapcar #'(lambda (x) (if x 1 0))
			       (list from-given-p to-given-p
				     by-given-p length-given-p))))
    (error "Exactly 3 of from-given-p, to-given-p, by-given-p and length-given-p
 need to be supplied."))
  (flet ((seq (from by length)
	   "Create a sequence with given parameters."
	   (let* ((from (coerce from type))
		  (by (coerce by type))
		  (result (make-array length :element-type type)))
	     (dotimes (i length)
	       (setf (aref result i) (+ from (* i by))))
	     result)))
    (cond
      ;; from, to, by => sequence, length
      ((not length-given-p)
       (let* ((range (- to from))
	      (length (1+ (floor (/ range by)))))
	 (unless (plusp (* by range))
	   (error "by needs to have the same sign as (- to from)."))
	 (values (seq from by length) length)))
      ;; from, to, length => sequence, by
      ((not by-given-p)
       (cond
	 ((and (= from to) (= length 1))
	  (values (make-array 1 :element-type type :initial-element from) 0))
	 ((and (/= from to) (> length 1))
	  (let ((by (/ (- to from) (1- length))))
	    (values (seq from by length) (coerce by type))))
	 (t (error "Mismatching from, to and length."))))
      ;; from, by, length => sequence, to
      ((not to-given-p)
       (let ((result (seq from by length)))
	 (values result (aref result (1- length)))))
      ;; to, by, length => sequence, from
      (t 
       (let ((from (- to (* by length))))
	 (values (seq from by length) (coerce from type)))))))

(define-modify-macro maxf (&rest values) max
		     "Replace with largest of the place and the given values")

(define-modify-macro minf (&rest values) min
		     "Replace with smallest of the place and the given values")

(define-modify-macro multf (&rest values) * "Multiply by the arguments")

(defun square (x)
  "Return the square of the argument."
  (* x x))

(defun positive-part (x)
  "Return x if (plusp x), 0 otherwise."
  (if (plusp x)
      x
      (coerce 0 (type-of x))))

(defun negative-part (x)
  "Return (- x) if (minusp x), 0 otherwise."
  (if (minusp x)
      (- x)
      (coerce 0 (type-of x))))
  
(defun twice (x)
  "Return two times the argument."
  (* x 2))

(defun half (x)
  "Return half of the argument."
  (/ x 2))

(defun mean (&rest xs)
  "Return the mean of the arguments."
  (if xs
      (iter
	(for x :in xs)
	(for n :from 1)
	(summing x :into sum)
	(finally 
	 (return (/ sum n))))
      (error "we need at least one argument")))

(defun arrays-range (&rest arrays)
  "Find the range of the arrays combined, returning two values, min
and max."
  (let (min max)
    (iter
      (for array :in arrays)
      (for (values array-min array-max) := (array-range array))
      (if (first-time-p)
	  (setf min array-min
		max array-max)
	  (setf min (min min array-min)
		max (max max array-max))))
    (values min max)))

(defun matrix-multiply (A B &optional (element-type (array-element-type A)))
  "Multiply matrices A and B.  If not given, the resulting element
type is the type of A."
  (assert (= 2 (array-rank A) (array-rank B)))
  (bind (((nrow-A ncol-A) (array-dimensions A))
	 ((nrow-B ncol-B) (array-dimensions B))
	 (result (make-ffa (list nrow-A ncol-B) element-type)))
    (assert (= ncol-A nrow-B))
    (dotimes (row nrow-A)
      (dotimes (col ncol-B)
	(setf (aref result row col)
	      (iter
		(for j :from 0 :below ncol-A)
		(summing (* (aref A row j) (aref B j col)))))))
    result))

(defun matrix-vector-multiply (A B &optional (element-type (array-element-type A)))
  "Perform matrix multiplication on A and B.  If A is a vector, it is
used as a column vector, if B is a vector, it is used as a column
vector.  If not given, the resulting element type is the type of A.
The result is a scalar, a vector, or a matrix, depending on the rank
of the arguments."
  (let ((rank-A (array-rank A))
	(rank-B (array-rank B)))
    (cond
      ((= 1 rank-A rank-B) (dot-product A B))
      ((= 2 rank-A rank-B) (matrix-multiply A B element-type))
      ((and (= 1 rank-A) (= 2 rank-B))
       (find-or-displace-to-flat-array 
	(matrix-multiply (displace-array A (list 1 (array-dimension A 0)) 0)
			 B element-type)))
      ((and (= 2 rank-A) (= 1 rank-B))
       (find-or-displace-to-flat-array 
	(matrix-multiply A (displace-array B (list (array-dimension B 0) 1) 0)
			element-type)))
      (t (error "A and B can only be vectors or matrices.")))))

(defun identity-matrix (n &optional (one 1d0))
  "Create an identity matrix, with one on the diagonal.  The element
type will be the type of one."
  (let* ((type (type-of one))
	 (id (make-ffa (list n n) type 
		      :initial-element (coerce 0 type))))
    (dotimes (j n id)
      (setf (aref id j j) one))))

(defun sum-of-squares (v)
  "Sum of squared elements of v."
  (array-sum v #'square))

(defun l2-norm (v)
  "L2 (aka Euclidean) norm."
  (sqrt (sum-of-squares v)))

(defun convex-combination (x y alpha)
  "Return alpha*x+(1-alpha)*y."
  (+ (* alpha x) (* (- 1 alpha) y)))

(defun find-point-approaching-boundary (boundary x predicate &key
					(alpha 0.5) (nmax 100))
  "Find the first x such that (predicate x), gradually approaching the
boundary making the new x the convex combination of the boundary and
the old one.  A maximum of nmax iterations are allowed."
  (iter
    (for i :from 0 :below nmax)
    (when (funcall predicate x)
      (return-from find-point-approaching-boundary x))
    (setf x (convex-combination boundary x alpha)))
  (error "reached maximum number of iterations")
  #| !!!! raising appropriate condition would be nicer |#)

;; (find-point-approaching-boundary 2 10 (lambda (x) (< x 3))) => 2.5
   
(defun quadratic-roots (a b c &key
			refuse-identical-roots refuse-complex-roots)
  "Calculate the roots of the quadratic polynomial a x^2 + b x + c.
You can request that an error is signalled for identical or complex
roots.

Will return the two roots and the type of the roots (real, complex,
identical).  The roots are sorted in increasing
order (lexicographically in real and imaginary components)."
  (assert (not (zerop a)))
  (let* ((discriminant (- (square b) (* a c 4)))
	 (type (cond
		 ((zerop discriminant)
		  (assert (not refuse-identical-roots))
		  'identical)
		 ((minusp discriminant)
		  (assert (not refuse-complex-roots))
		  'complex)
		 (t 'real)))
	 (sqrt-discriminant (sqrt discriminant))
	 (r1 (/ (- 0 b sqrt-discriminant) a 2))
	 (r2 (/ (- sqrt-discriminant b) a 2)))
    (if (plusp a)
	(values r1 r2 type)
	(values r2 r1 type))))
