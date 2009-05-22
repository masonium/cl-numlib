(in-package :cl-numlib)

;; multivariate optimization 

;; Calling convention: (f x &optional derivatives-p).  If
;; derivatives-p, return 2 values: a scalar f(x), and a vector f'(x).
;; If derivatives are not requested, then only the first value may be
;; expected.

;; Line search functions: return an alpha for scaling direction

(declaim (optimize (debug 3) (safety 3)))

(define-condition linesearch-reached-maximum-iterations-error (error)
  ;; !!! might want to give some useful info here
  ())

(defun linesearch (f x0 direction &key
		   (alpha 1d0) (a1 0.1d0) (a2 1.5d0)
		   (m 0.0001) (m1 0.9) (maxiter 50) f0 fd0
		   (alpha-right nil))
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
- alpha-right: the right boundary for the linesearch, nil means infinity

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
	 (hd0 (dot-product fd0 direction)) ; h'(0)
	 alpha-prev			; previous value of alpha
	 (x (make-array n :element-type 'double-float))	; x=x0+alpha*direction,
							; we will
							; reuse it
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
			       ;; will try cubic, if numerically
			       ;; unstable, will revert to quadratic-or-midpoint
			       (handler-case 
				   (let* ((k-prev (- h-alpha-prev f0 
						     (* hd0 alpha-prev)))
					  (alpha-square (square alpha))
					  (alpha-prev-square (square alpha-prev))
					  (denominator (* alpha-square 
							  alpha-prev-square
							  (- alpha alpha-prev)))
					  (a (/ (- (* alpha-prev-square k)
						   (* alpha-square k-prev))
						denominator))
					  (b (/ (- (* alpha-square alpha k-prev)
						   (* alpha-prev-square 
						      alpha-prev k))))
					  (discriminant (- (square b) (* 3 a hd0))))
;;;				 (format t "trying cubic, a=~a~%" a)
				     (if (or (zerop a) (minusp discriminant))
					 (quadratic-or-midpoint)
					 (/ (- (sqrt discriminant) b)
					    (* 3 a))))
				 (arithmetic-error ()
				   (quadratic-or-midpoint)))
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
	(unless h-alpha
	  ;; note: x is not copied, as we are not supposed to recover
	  ;; this error with a restart
	  (error 'function-evaluated-at-infeasible-point :f f :x x
		 :alpha alpha :x0 x0 :direction direction))
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
      (error 'linesearch-reached-maximum-iterations-error))))

(defun linesearch-robust (f x0 direction &key
			  (alpha 1d0) (a1 0.1d0) (a2 1.5d0)
			  (m 0.0001) (m1 0.9) (maxiter 50) f0 fd0
			  (alpha-right nil))
  "Robust linesearch: if plain vanilla linesearch fails, try golden
section minimization.  Arguments and return values are the same as
with linesearch."
  (handler-case (linesearch f x0 direction :alpha alpha :a1 a1 :a2 a2
			    :m m :m1 m1 :maxiter maxiter :f0 f0 :fd0 fd0
			    :alpha-right alpha-right)
    (linesearch-reached-maximum-iterations-error ()
;;      (format t "performing golden section search
;;x0=~a direction=~a alpha=~a~%" x0 direction alpha)
      (flet ((g (alpha)
;;	       (format t "gs alpha=~a~%" alpha)
	       (funcall f (array+ x0 (array-scalar* direction alpha)) nil)))
	(bind ((alpha (golden-section-minimize #'g 0d0 alpha (* alpha 1d-5)))
	       (x (array+ x0 (array-scalar* direction alpha)))
	       ((values fx fdx) (funcall f x t)))
	  (values alpha x fx fdx))))))
	  
    
(defun bfgs-update-H (H dx dg &key (minimum-dxdg 1d-12) (H-scaling 1d0))
  "Update the inverse Hessian matrix used by the BFGS method.  dx is
the change in x, dg is the change in g=f'(x).  If their cross product
is below minimum-dxdg, return an identity-matrix multiplied by
H-scaling, except if H-scaling is nil, then it is left alone."
  (let ((dxdg (dot-product dx dg)))
	(if (<= dxdg minimum-dxdg)
	    ;; when dxdg is very small, H is reset to identity
	    (if H-scaling
		(identity-matrix (length dx) H-scaling)
		H)
	    ;; otherwise, we just proceed with updating H
	    (let* ((dgHdg (quadratic-form dg H))
		   (dxdgH (matrix-vector-multiply (outer-product dx dg) H)))
	      (array-
	       (array+ H (array-scalar* (outer-product dx dx)
					(/ (+ 1d0 (/ dgHdg dxdg)) dxdg)))
	       (array-scalar/ (array+ dxdgH (matrix-transpose dxdgH)) dxdg))))))

(define-condition bfgs-H-not-positive-definite-error (error)
  ((H :initarg :H)))

(define-condition bfgs-stuck-error (error)
  ((H :initarg :H)
   (x :initarg :x)
   (fx :initarg :fx)
   (fdx :initarg :fdx)
   (fx-next :initarg :fx-next)))
	  
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
Springer, 2006.  Chapter 6.  Heuristic for initializing H is from p.
143."
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
      (let ((direction (array-negate (matrix-vector-multiply H fdx))))
	(when (plusp (dot-product direction fdx))
	  (error 'bfgs-H-not-positive-definite-error :H H))
	(bind (((values alpha x-next fx-next fdx-next)
		(linesearch-robust f x direction :a1 linesearch-a1 :a2 linesearch-a2
		 :m linesearch-m :m1 linesearch-m1 :maxiter linesearch-maxiter 
		 :f0 fx :fd0 fdx))
	       (dx (array-scalar* direction alpha)) ; change in x
	       (dg (array- fdx-next fdx)))	    ; change in f'(x)
	  (let ((ssdg (sum-of-squares dg)))
	    (setf H (if (and (zerop iteration) (not H-scaling-p) (not (zerop ssdg)))
			(identity-matrix n (/ (dot-product dx dg) ssdg))
			(bfgs-update-H H dx dg :h-scaling #|(if H-scaling-p
							      H-scaling
							      1d0)|# nil))))
	  ;; check if the iteration is stuck
	  (if (>= fx-next fx)
	      (error 'bfgs-stuck-error :H H :x x :fx fx :fdx fdx :fx-next fx-next))
	  ;; now we update
	  (setf x x-next
		fx fx-next
		fdx fdx-next)
	  ;; check if we are done
	  (when (or (< (l2-norm dx) (* epsilon (1+ (l2-norm x))))
		    (< (l2-norm fdx) (* delta (1+ fx))))
	    (return-from bfgs-minimize 
	      (values x fx fdx iteration))))))
    (warn "bfgs reached maximum number of iterations")
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
