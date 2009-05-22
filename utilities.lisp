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
	 ((and (= from to))
	  (values (make-array length :element-type type :initial-element 
			      (coerce from type)) 0))
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

(define-modify-macro maxf* (&rest values) (lambda (&rest values)
					    (reduce (lambda (x y)
						      (cond
							((and x y) (max x y))
							(x x)
							(t y)))
						    values))
		     "Replace with largest of the place and the given
		     values, ignoring nil's")

(define-modify-macro minf* (&rest values) (lambda (&rest values)
					    (reduce (lambda (x y)
						    (cond
						      ((and x y) (min x y))
						      (x x)
						      (t y)))
						  values))
		     "Replace with largest of the place and the given
		     values, ignoring nil's")


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

(defun absolute-difference (a b)
  "Return |a-b|."
  (abs (- a b)))

(defun arrays-range (&rest arrays)
  "Find the range of the arrays combined, returning two values, min
and max."
  (let (min max)
    (iter
      (for array :in arrays)
      (for (array-min array-max) := (array-range array))
      (if (first-time-p)
	  (setf min array-min
		max array-max)
	  (setf min (min min array-min)
		max (max max array-max))))
    (values min max)))

(defun matrix-transpose (matrix)
  "Tranpose a matrix."
  ;; !!!! efficiency note: could walk one of the matrices as a flat
  ;; !!!! vector
  (assert (= 2 (array-rank matrix)))
  (bind (((n m) (array-dimensions matrix))
	 (result (make-array (list m n) :element-type (array-element-type matrix))))
    (dotimes (i n)
      (dotimes (j m)
	(setf (aref result j i) (aref matrix i j))))
    result))

(defun matrix-multiply (A B &optional (element-type (array-element-type A)))
  "Multiply matrices A and B.  If not given, the resulting element
type is the type of A."
  (assert (= 2 (array-rank A) (array-rank B)))
  (bind (((nrow-A ncol-A) (array-dimensions A))
	 ((nrow-B ncol-B) (array-dimensions B))
	 (result (make-array (list nrow-A ncol-B) :element-type element-type)))
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
used as a row vector, if B is a vector, it is used as a column vector.
If not given, the resulting element type is the type of A.  The result
is a scalar, a vector, or a matrix, depending on the rank of the
arguments."
  (let ((rank-A (array-rank A))
	(rank-B (array-rank B)))
    (cond
      ((= 1 rank-A rank-B) (dot-product A B))
      ((= 2 rank-A rank-B) (matrix-multiply A B element-type))
      ((and (= 1 rank-A) (= 2 rank-B))
       (flatten-array 
	(matrix-multiply (displace-array A (list 1 (array-dimension A 0)) 0)
			 B element-type)))
      ((and (= 2 rank-A) (= 1 rank-B))
       (flatten-array 
	(matrix-multiply A (displace-array B (list (array-dimension B 0) 1) 0)
			element-type)))
      (t (error "A and B can only be vectors or matrices.")))))


(defun quadratic-form (x A)
  "Calculate the quadratic form x^T A x.  The symmetry of A is not assumed."
  ;; Note: this the simplest algorithm, unoptimized, you are welcome
  ;; to improve it.
  (assert (vectorp x))
  (let ((n (length x))
	(sum 0))
    (assert (and (typep A 'array) (equal (array-dimensions A) (list n n))))
    (dotimes (i n)
      (dotimes (j n)
	(incf sum (* (aref x i) (aref x j) (aref A i j)))))
    sum))

(defun identity-matrix (n &optional (one 1d0))
  "Create an identity matrix, with one on the diagonal.  The element
type will be the type of one."
  (let* ((type (type-of one))
	 (id (make-array (list n n) :element-type type 
			 :initial-element (coerce 0 type))))
    (dotimes (j n id)
      (setf (aref id j j) one))))

(defun sum-of-squares (v)
  "Sum of squared elements of v."
  (array-sum v :key #'square))

(defun l1-norm (v)
  "L1 (sum of absolute values) norm."
  (array-sum v :key #'abs))

(defun l2-norm (v)
  "L2 (aka Euclidean) norm."
  (sqrt (sum-of-squares v)))

(defun sup-norm (v)
  "Supremum (maximum absolute value) norm."
  (array-max v :key #'abs))


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
identical) as (values r1 r2 type).  The roots are sorted in increasing
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

(defun solve-2x2-system (X z)
  "Solve the 2x2 linear system Xy=z, returning y.  If X is singular,
an error is signalled."
  (declare ((array * (2 2)) X) ((array * (2)) z))
  (bind ((#2A((a b)
	      (c d)) X)
	 (#(e f) z)
	 (determinant (- (* a d) (* b c))))
    (when (zerop determinant)
      (error "The linear system is singular."))
    (vector (/ (- (* d e) (* b f)) determinant)
	    (/ (- (* a f) (* c e)) determinant))))
  
(defmacro iter-indexing-vector ((vector var index) &body body)
  "An iter form that goes assigns the values in vector to var, also
running the index from 0."
  `(iter
     (for ,var :in-vector ,vector)
     (for ,index :from 0)
     ,@body))

(declaim (inline -+ -*))

(defun -+ (&rest args)
  "Shorthand for (- (+ ...))."
  (- (apply #'+ args)))

(defun -* (&rest args)
  "Shorthand for (- (* ...))."
  (- (apply #'* args)))

(defun displace-matrix-row (matrix row)
  "Return a the given row of a matrix as a displaced array."
  (assert (typep matrix '(array * (* *))))
  (bind (((nrow ncol) (array-dimensions matrix)))
    (assert (< -1 row nrow))
    (displace-array matrix ncol (* row ncol))))

(defun slice-matrix-rows (matrix)
  "Return matrix as a list of rows."
  (assert (typep matrix '(array * (* *))))
  (bind (((nrow ncol) (array-dimensions matrix)))
    (iter
      (for r :from 0 :below nrow)
      (collect (displace-array matrix ncol (* r ncol))))))

(defun normalize-vector (vector sum)
  "Divide each element of a nonnegative vector so that they will sum
to `sum', return the resulting vector."
  (assert (vectorp vector))
  (let ((cumsum 0d0))
    (iter
      (for v :in-vector vector)
      (assert (<= 0 v))
      (incf cumsum v))
    (when (zerop cumsum)
      (error "vector sums to zero"))
    (array-scalar/ vector (/ cumsum sum))))

(defun normalize-probability (vector)
  "Divide each element of a nonnegative vector so that they will sum
to 1, return the resulting vector."
  (normalize-vector vector 1d0))

(defun first-difference (vector &optional (coeff 1))
  "First-difference a vector a[0],a[1],... and return the result
a[0]-a[1],a[1]-a[2],..., multiplied by coeff (will be one element
shorter than the original, of course).  Setting coeff -1 allows first
differencing in the other direction, resulting in
a[1]-a[0],a[2]-a[1],..."
  (assert (and (vectorp vector) (< 1 (length vector))))
  (let* ((n (length vector))
	 (result (make-array (1- n) :element-type (array-element-type vector))))
    (dotimes (i (1- n))
      (setf (aref result i)
	    (* coeff (- (aref vector i) (aref vector (1+ i))))))
    result))

(defmacro conditions-to-alternative-value ((value conditions) &body body)
  "Evaluate body and return its value.  In case a condition is
signalled, and it is one of the conditions given as the argument,
value is returned."
  (with-unique-names (condition)
    `(handler-case (progn
		     ,@body)
     (t (,condition)
       (if (member ,condition ,conditions
		   :test (function typep))
	   ,value
	   (error ,condition))))))

;; (defun array-finite-difference (array coordinate delta
;; 				&key (alternative-value nil)
;; 				(result-element-type t)
;; 				(conditions-to-ignore '(simple-type-error)))
;;   "Calculate finite differences for array along the given coordinate,
;; by stepping `delta' (eg 1,-1,etc) elements.  For elements where we
;; would run outside the array boundaries, or where one of
;; `conditions-to-ignore' is signalled, `alternative-value' is used.  The
;; resulting array conforms to the dimensions of the original one and has
;; `result-element-type' as its element-type.  Please make sure that
;; `alternative-value' and `result-element-type' are compatible."
;;   (let* ((dimensions (array-dimensions array))
;;          (upper-bound (nth coordinate dimensions))
;;          (indexer (row-major-indexing:make-rm-indexer dimensions))
;;          (reindexer (row-major-indexing:make-rm-reindexer dimensions))
;;          (total-dimension (reduce #'* dimensions))
;;          (result (make-array dimensions :element-type result-element-type))
;; 	 (conditions (cons 'row-major-indexing:modified-index-outside-range-error
;; 			   conditions-to-ignore)))
;;     (dotimes (i total-dimension)
;;       (let* ((indexes (funcall reindexer i))
;;              (difference
;;               (conditions-to-alternative-value (alternative-value conditions)
;; 		(- (row-major-aref array 
;; 				   (funcall indexer 
;; 					    (row-major-indexing:modify-nth-index
;; 							      indexes
;; 							      coordinate
;; 							      delta
;; 							      upper-bound)))
;; 		   (row-major-aref array i)))))
;;         (setf (row-major-aref result i) difference)))
;;     result))

(defun extend-vector (vector direction &optional (initial-element 0))
  "If `direction' is positive, append that many elements at the end of
the vector, of it is negative, append (- direction) elements at the
start of the vector, copy the original elements, fill the rest with
initial-element, and return the resulting vector."
  (assert (and (vectorp vector) (integerp direction)))
  (let ((n (length vector))
	(element-type (array-element-type vector)))
    (cond
      ((plusp direction)
       (let ((result (make-array (+ n direction) 
				 :element-type element-type
				 :initial-element initial-element)))
	 (iter
	   (for v :in-vector vector)
	   (for i :from 0)
	   (setf (aref result i) v))
	 result))
      ((minusp direction)
       (let ((result (make-array (- n direction) 
				 :element-type element-type
				 :initial-element initial-element)))
	 (iter
	   (for v :in-vector vector)
	   (for i :from (- direction))
	   (setf (aref result i) v))
	 result))
      (t (error "direction has to be nonzero")))))

(defun pretty-axis-positions (lower upper n &key 
			      (round-to '(1 2 5))
			      (include-lower-p nil)
			      (include-upper-p nil))
  "Supply suitably rounded numbers between `lower' and `upper', with
at most `n' (or, in rare cases, n+1) numbers.  Returns multiple
values, a list with the positions and an exponent which can be used
for rounding these numbers.  When lower=upper, (values (lower) nil) is
returned.  round-to should be a list of numbers that divide 10,
mantissas are rounded to these values.

Arguments `include-lower-p' and `include-upper-p' can be used to
specify whether to include endpoints.

Examples:
 (pretty-axis-positions 1 2.1 2) => (1 3/2 2 5/2), -1
 (pretty-axis-positions 417 972 3) => (400, 500, 600, 700, 800, 900, 1000), 2
"
  ;; convert endpoints to rational numbers
  (let ((lower (rationalize lower))
	(upper (rationalize upper)))
    ;; check input
    (assert (<= 0 n) (n) "n=~A is not positive" n)
    (assert (<= lower upper) (lower upper)
	    "~A <= ~A does not hold." lower upper)
    ;; trivial axis: return single point
    (when (= lower upper)
      (return-from pretty-axis-positions (values (list lower) nil)))
    ;; non-trivial axis
    (let* ((raw-step (/ (- upper lower) (1+ n)))
	   (exponent (floor (log raw-step 10)))
	   (mantissa (/ raw-step (expt 10 exponent)))
	   (step-first-digit (or (find-if (lambda (s) (<= mantissa s))
				 round-to) 10))
	   (step (* step-first-digit
		    (expt 10 exponent)))
	   (start (* (if include-lower-p
			 (floor lower step)
			 (ceiling lower step))
		     step))
	   (end (* (if include-upper-p
		       (ceiling upper step)
		       (floor upper step))
		   step)))
      ;; correct exponent if step-first-digit is 10
      (if (= step-first-digit 10)
	  (incf exponent))
      ;; collect marks
      (iter
	(for label-pos :from start :to end :by step)
	(collecting label-pos :into positions :result-type 'vector)
	(finally
	 (return (values positions exponent)))))))

(defun vector-last (vector &optional (count 1))
  "Return length-count'th (last if count=1) element of a vector."
  (declare ((vector) vector))
  (aref vector (- (length vector) count)))

(defun make-gamma-transform (gamma)
  "Return a function that performs a gamma correction with the given
coefficient on its argument, on the interval [0,1]."
  (if (= 1 gamma)
      #'identity
      (progn
	(assert (plusp gamma))
	(lambda (x)
	  (assert (<= 0 x 1))
	    (expt x gamma)))))
  
