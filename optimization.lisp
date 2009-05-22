(in-package :cl-numlib)

(defun golden-section-combination (a b)
  "Return the convex combination (1-G)*a+G*b, where G is the
inverse of the golden ratio."
;;  (declare (double-float a b))
;;   (declare (optimize (speed 3) (safety 1)
;;                      (compilation-speed 0)
;;                      (space 0) (debug 0)))
  (+ (* #.(- 1d0 (/ (- 3d0 (sqrt 5d0)) 2d0)) a)
     (* #.(/ (- 3d0 (sqrt 5d0)) 2d0) b)))

(defun golden-section-minimize (f a b tol &optional (nmax 100))
  "Find a local minimum of f in the [a,b] interval.  The algorithm
terminates when the minimum is bracketed in an interval smaller than
tol.  Since the algorithm is slow, tol should not be chosen smaller
then necessary.  The algorithm will also find the local minimum at the
endpoints, and if f is unimodal, it will find the global minimum.
nmax is there for terminating the algorithm, in case tolerance is zero
or too small.  All values should be double-float, and f should also
return double-floats.

Note: when f is constant on a range, golden-section-minimize ``pulls
to the left'', ie will keep picking smaller values."
  (declare (double-float a b tol)
	   (fixnum nmax)
	   (type (function (double-float) double-float) f)
	   (inline golden-section-combination)
	   (optimize (speed 3) (safety 1)
                     (compilation-speed 0)
                     (space 0) (debug 0)))
  ;; reorder a and b if necessary
  (when (> a b)
    (rotatef a b))
  ;; start iteration with golden ratio inner points
  (let* ((m1 (golden-section-combination a b))
	 (m2 (golden-section-combination b a))
	 (f1 (funcall f m1))
	 (f2 (funcall f m2))
	 (n 0))
    (declare (double-float m1 m2 f1 f2))
    (do ()
	((or (<= (abs (- b a)) tol)
	     (>= n nmax))
	 (if (< f1 f2)			; change < to maximize
	     (values m1 f1)
	     (values m2 f2)))
;;;; uncomment below for debugging
;;       (incf n)
;;       (format t "bracket is a=~a~%m1=f(~a)=~a~%m2=f(~a)=~a~%b=~a~%"
;;  	      a m1 f1 m2 f2 b)
      (if (<= f1 f2)			; change < to maximize
	  (progn 
	    ;; new bracket is (a,m1,m2)
	    (shiftf b m2 m1 (golden-section-combination m1 a))
	    (shiftf f2 f1 (funcall f m1)))
	  (progn
	    ;; new bracket is (m1,m2,b)
	    (shiftf a m1 m2 (golden-section-combination m2 b))
	    (shiftf f1 f2 (funcall f m2)))))))

(define-condition univariate-outside-boundaries (error)
  ((x :initarg :x)
   (a :initarg :a)
   (b :initarg :b))
  (:report (lambda (condition stream)
	     (format stream "~a is outside [~a,~a]"
		     (slot-value condition 'x)
		     (slot-value condition 'a)
		     (slot-value condition 'b)))))

(define-condition reached-maximum-iterations (error)
  ((nmax :initarg :nmax)
   (x :initarg :x))
  ;;; !!! write report function
)

(define-condition univariate-found-local-maximum (error)
  ((x :initarg :x)
   (fp :initarg :fp)
   (fpp :initarg :fpp))
  (:report (lambda (condition stream)
	     (format stream "~a is not a local minimum: f''=~a, f'=~a"
		     (slot-value condition 'x)
		     (slot-value condition 'fpp)
		     (slot-value condition 'fp)))))

(defun newton-minimize (fp x a b &optional (tol 1d-5) (epsilon 1d-5) (nmax 100))
  "Use the Newton-Raphson method for finding the minimum of a
function.  Iteration stops when the change in x is below tol, or the
slope of f is below epsilon in absolute value.  The algorithm allows a
maximum of nmax iterations in order to avoid going in circles.  f is
supposed to return 2 values: the first and second derivatives. 

If you don't care about boundaries, set a and/or b to nil.

Returns two values, x and f''."
  (let ((n 0))
    (tagbody 
     top
       (when (or (and a (< x a)) (and b (> x b)))
	 (error 'univariate-outside-boundaries :x x :a a :b b))
       (when (>= n nmax)
	 (error 'reached-maximum-iterations :nmax nmax :x x))
       (multiple-value-bind (fp fpp) (funcall fp x)
	 (when (zerop fpp)
	   ;; we recover from f''=0 by using f''=-1
	   (setf fpp -1))
	 (let ((delta (/ fp fpp)))
	   (when (and (<= (abs delta) tol)
		      (<= (abs fp) epsilon))
	     (if (plusp fpp)
		 (return-from newton-minimize (values x fpp))
		 (error 'univariate-found-local-maximum :x x :fp fp :fpp fpp)))
	   (decf x delta)))
       (incf n)
       (go top))))

(defun negate-univariate-function (f)
  "Return the function (-f).  Useful for calling
golden-section-minimize for maximization."
  (lambda (x)
    (- (funcall f x))))

(defun negate-univariate-fp (fp)
  "Return a function that returns values -f'(0),-f''(0)  where fp returned 
f'(0) and f''(0).  Useful for calling newton-minimize for maximization."
  (lambda (x)
    (multiple-value-bind (fp fpp) (funcall fp x)
      (values (- fp) (- fpp)))))

