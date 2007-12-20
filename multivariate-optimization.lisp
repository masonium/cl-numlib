(in-package :cl-numlib)

;; multivariate optimization 

;; Calling convention: (f x &optional derivatives-p).  If
;; derivatives-p, return 2 values: a scalar f(x), and a vector f'(x).
;; If derivatives are not requested, then only the first value may be
;; expected.

;; Line search functions: return an alpha for scaling direction

(declaim (optimize (debug 3) (safety 3)))

(defun linesearch (f x0 direction &key
		   (alpha 1d0) (a1 0.1d0) (a2 1.5d0)
		   (m 0.0001) (m1 0.9) (maxiter 50) f0 fd0)
  "Perform line search on a h(alpha)=f(x0+alpha*direction).  You can
supply f0=f(x0) and fd0=f'(x0) for speed, otherwise they will be
obtained by calling f.  The algorithm terminates when both of the
conditions below (known as weak Wolfe conditions) hold, or the number
of steps reaches maxiter:

h(alpha) <= h(0) + m h'(0) alpha
h'(alpha) >= m1 h'(0)

Parameters:

- alpha: starting alpha, usually 1 for quasi-Newton methods
- a1: when shrinking finite intervals, determines padding on sides, 0 < a1 < 0.5.
- a2: extrapolation multiplier for shrinking non-finite intervals, a2 > 1
- m, m1: see above, the smoother the function, the larger m1 should be,
  0 < m < m1 < 1.  Defaults are usually good for BFGS.
- maxiter: if conditions are not satisfied in this many iterations, 
  a warning is given

Return values:

- alpha
- x=x0+direction*alpha
- f(x)
- f'(x)

References:

LEMARECHAL, C. 1981. A view of line-searches. In Optimization and
Optimal Control, A.  Auslender, W. Oettli, and J. Steer, Eds. Lecture
Notes in Control and Information Science, vol.  30. Springer-Verlag,
Heidelberg, 59-78."
  ;; check consistency of parameters
  (assert (< 0 a1 0.5))
  (assert (< 1 a2))
  (assert (< 0 m m1 1))
  (assert (plusp alpha))
  ;; if not given, calculate f and f' at x0
  (unless (and f0 fd0)
    (setf (values f0 fd0) (funcall f x0 t)))
  ;; 
  (let* ((n (length x0))
	 (alpha-left 0d0)		; left boundary
	 (alpha-right nil)		; right boundary, nil means infinity
	 (hd0 (dot-product fd0 direction)) ; h'(0)
	 alpha-prev			; previous value of alpha
	 (x (make-ffa n :double))	; x=x0+alpha*direction, we will reuse it
	 fdx				; f'(x) goes here
	 h-alpha hd-alpha h-alpha-prev hd-alpha-prev) ; h, h'
    (assert (minusp hd0))
    (flet ((interpolate (left right)
;; 	     (format t "interpolating on [~a,~a], alpha=~a h-alpha=~a,
;;  hd-alpha=~a, alpha-prev=~a h-alpha-prev=~a hd-alpha-prev=~a~%" left
;;  right alpha h-alpha hd-alpha alpha-prev h-alpha-prev hd-alpha-prev)
	     (assert (< left right))
	     (let ((padding (* a1 (- right left)))
		   (k (- h-alpha f0 (* hd0 alpha))))
	       (max (+ left padding)
		    (min (- right padding)
			 (flet ((quadratic-or-midpoint ()
;;				  (format t "quadratic-or-midpoint, k=~a~%" k)
				  (if (zerop k)
				      (half (+ left right)) ; midpoint
				      (/ (* hd0 (square alpha)) -2d0 k))))
;;			   (format t "f(alpha)=~a f(alpha-prev)=~a~%"
;; 				   (funcall f (array+
;; 					       x0 (array-scalar* direction alpha)))
;; 				   (if alpha-prev
;; 				       (funcall f (array+
;; 						   x0 (array-scalar* direction
;; 								     alpha-prev)))))
			   (if alpha-prev
			       ;; will try cubic
			       (let* ((k-prev (- h-alpha-prev f0 (* hd0 alpha-prev)))
				      (alpha-square (square alpha))
				      (alpha-prev-square (square alpha-prev))
				      (denominator (* alpha-square alpha-prev-square
						      (- alpha alpha-prev)))
				      (a (/ (- (* alpha-prev-square k)
					       (* alpha-square k-prev))
					   denominator))
				      (b (/ (- (* alpha-square alpha k-prev)
					       (* alpha-prev-square alpha-prev k))))
				      (discriminant (- (square b) (* 3 a hd0))))
;;				 (format t "trying cubic, a=~a~%" a)
				 (if (or (zerop a) (minusp discriminant))
				     (quadratic-or-midpoint)
				     (/ (- (sqrt discriminant) b)
					(* 3 a))))
			       ;; can only do quadratic on a single point
			       (quadratic-or-midpoint)))))))
	   (return-serious ()
	     ;; return normally
	     (return-from linesearch (values alpha x h-alpha fdx)))
	   (update-and-save-prev (alpha-next)
	     ;; save previous values
	     (setf alpha-prev alpha
		   h-alpha-prev h-alpha
		   hd-alpha-prev hd-alpha
		   alpha alpha-next)))
      (dotimes (i maxiter)
	;; calculate x = x0+alpha*direction, evaluate f(x), f'(x),
	;; h(alpha), h'(alpha)
	(dotimes (j n)
	  (setf (aref x j) (+ (aref x0 j) (* alpha (aref direction j)))))
	(setf (values h-alpha fdx)
	      (funcall f x t))
	(setf hd-alpha (dot-product direction fdx))
	;; shrink interval
;;	(format t "alpha=~a h(alpha)=~a h'(alpha)=~a~%" alpha h-alpha hd-alpha)
	(cond
	  ((< h-alpha (+ f0 (* hd0 m alpha)))
	   (when (>= hd-alpha (* hd0 m1)) ; serious conditions hold
		  (return-serious))
	   ;; too large
;;	   (format t "too large~%")
	   (setf alpha-left alpha)
	   (update-and-save-prev
	    (if alpha-right
		(interpolate alpha-left alpha-right) ; finite, interpolate
		(* alpha-left a2))))		     ; infinite, extrapolate
	  (t				; too small
;;	   (format t "too small~%")
	   (setf alpha-right alpha)
	   (update-and-save-prev (interpolate alpha-left alpha-right)))))
      (warn "reached maximum number of iterations")
      (return-serious))))

	  
(defun bfgs-minimize (f x0 &key (maxiter 100)
		      (epsilon 1d-8) (delta 1d-5)
		      (H-scaling 1d0 H-scaling-p)
		      (linesearch-a1 0.1d0) (linesearch-a2 1.5d0)
		      (linesearch-m 0.0001) (linesearch-m1 0.9)
		      (linesearch-maxiter 50))
  "Use the BGFS method to find the minimum of f, starting from x0.
The following parameters are optional:

- maxiter: maximum number of iterations
- epsilon: don't stop before change in x is below this
- delta: don't stop before L2 norm of gradient is below this

- linesearch-a1, linesearch-a2, linesearch-m, linesearch-m1,
  linesearch-maxiter: see linesearch

- H-scaling: a value for the diagonal of the initial scale of the
  inverse of the Hessian.  If not given, a heuristic is used that
  performs quite well in most cases.

Return the following values:
- x, the optimum
- f(x)
- number of iterations

References:

Nocedal, Jorge and Wright, Stephen J.  Numerical Optimization.
Springer, 1999.  Chapter 8.  Heuristic for initializing H is from p.
200."
;;;  (declare (ignore perturb-initial-critical-point))
  (bind ((n (length x0))
	 (H (identity-matrix n H-scaling))
	 (x x0)
	 ((values fx fdx) (funcall f x t)))
    ;; we may not believe that the critical point is an optimum
    (when (< (l2-norm fdx) (* delta (1+ fx)))
      (return-from bfgs-minimize (values x fx 0)))
    ;; iterate to convergence
    (dotimes (iteration maxiter)
      (bind ((direction (array-negate (matrix-vector-multiply H fdx)))
	     ((values alpha x-next fx-next fdx-next)
	      (linesearch f x direction
	       :a1 linesearch-a1 :a2 linesearch-a2 :m linesearch-m 
	       :m1 linesearch-m1 :maxiter linesearch-maxiter 
	       :f0 fx :fd0 fdx))
	     (y (array-scalar* direction alpha)) ; change in x
	     (s (array- fdx-next fdx))		 ; change in f'(x)
	     (ys (dot-product y s)))		 ; <y,s>
	;; when ys is very small, H is reset to identity
	(if (<= (abs ys) (* (l2-norm y) (l2-norm s)
			    least-positive-double-float))
	    (setf H (identity-matrix n))
	    ;; otherwise, we just proceed with updating A
	    (let* ((w (array- y (matrix-vector-multiply H s)))
		   (H-update (array-scalar* (outer-product y y)
					    (- (/ (dot-product w s)
						  ys)))))
	      ;; subtract w y^T and its transpose from H-update
	      (dotimes (i n)
		(dotimes (j n)
		  (incf (aref H-update i j)
			(+ (* (aref w i) (aref y j))
			   (* (aref w j) (aref y i))))))
	      ;; use heuristic if needed
	      (when (and (zerop iteration) (not  H-scaling-p))
		(setf H (identity-matrix n (/ ys (sum-of-squares y)))))
	      ;; update 
	      (setf H (array+ H (array-scalar/ H-update ys)))))
	;; check if the iteration is stuck
	(if (>= fx-next fx)
	    (warn "the quasi-Newton algorithm appears to be stuck."))
	;; now we update
	(setf x x-next
	      fx fx-next
	      fdx fdx-next)
	;; check if we are done
	(when (or (< (l2-norm y) (* epsilon (1+ (l2-norm x))))
		  (< (l2-norm fdx) (* delta (1+ fx))))
	  (return-from bfgs-minimize 
	    (values x fx fdx iteration)))))
    (warn "reached maximum number of iterations")
    (return-from bfgs-minimize (values x fx fdx nil))))

(defun negate-multivariate-objective-function (f)
  "Return a function (g x &optional derivative-p), where g(x)=-f(x)
and g'(x)=-f(x).  Useful if you want to use bfgs-minimize for
maximization."
  (lambda (x &optional derivative-p)
    (if derivative-p
	(multiple-value-bind (fx fdx) (funcall f x t)
	  (values (- fx) (array-negate fdx)))
	(- (funcall f x nil)))))
