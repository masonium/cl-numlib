(defun golden-section-minimize (f a b tol &key 
				(fa nil fa-given-p)
				(fb nil fb-given-p)
				(m nil m-given-p)
				(fm nil fm-given-p))
  (when (and (not m-given-p) fm-given-p)
    (error "You can't supply fm without m!"))
  (unless m-given-p
    (setf m (golden-section-combination a b)))
  ;; !!! use a macro here
  (unless fa-given-p
    (setf fa (funcall f a)))
  (unless fb-given-p
    (setf fb (funcall f b)))
  (unless fm-given-p
    (setf fm (funcall f m)))
  ;; reorder a and b if necessary
  (when (> a b)
    (rotatef a b)
    (rotatef fa fb))
  ;; check bracketing
  (unless (< fm (min fa fb))
    (error "f(~a)=~a and f(~a)=~a do not bracket a minimum at f(~a)=~a."
	   a fa b fb m fm))
  ;; calculate the other inner point
  (multiple-value-bind (m1 f1 m2 f2) 
      (if (< (- m a) (- b m))
	  (let ((m2 (golden-section-combination m b)))
	    (values m fm m2 (funcall f m2)))
	  (let ((m1 (golden-section-combination m a)))
	    (values m1 (funcall f m1) m fm)))
    (do ()
	((< (abs (- b a)) 
	    (* tol (+ (abs m1) (abs m2))))
	 (if (< f1 f2)
	     (values m1 f1)
	     (values m2 f2)))
;;       (format t "bracket is~%f(~a)=~a~%f(~a)=~a~%f(~a)=~a~%f(~a)=~a~%"
;; 	      a fa m1 f1 m2 f2 b fb)
      (if (< f1 f2)
	  (progn 
	    ;; new bracket is (a,m1,m2)
	    (shiftf b m2 m1 (golden-section-combination m1 a))
	    (shiftf f2 f1 (funcall f m1)))
	  (progn
	    ;; new bracket is (m1,m2,b)
	    (shiftf a m1 m2 (golden-section-combination m2 b))
	    (shiftf f1 f2 (funcall f m2)))))))
(defun linesearch-armijo (f x0 direction &key
			  (alpha 1d0) (rho 0.75d0) (c 0.01)
			  (maxiter 50) f0 fd0)
  "Very simple Armijo-style linesearch, terminates when

f(x0+alpha*direction) >= f(x0) + c alpha < f'(x0), direction >

If this is not satisfied, alpha is set to alpha*rho.  Maxiter is the
maximum possible number of iterations, if reached, we signal an error.
You can supply f0=f(x0) and fd0=f'(x0) for speed, otherwise they are
calculated."
  (assert (and (plusp c) (< c 1d0) (plusp rho) (< rho 1d0) (plusp alpha)))
  (unless (and f0 fd0)
    (setf (values f0 fd0) (funcall f x0 t)))
  (let* ((n (length x0))
	 (tolerance (* c (dot-product direction fd0))))
    (flet ((g (alpha)
	     "Calculate f(x0+alpha*direction)."
	     (let ((x-alpha (make-ffa n :double)))
	       (dotimes (i n)
		 (setf (aref x-alpha i)
		       (+ (aref x0 i) (* alpha (aref direction i)))))
	       (funcall f x-alpha))))
      (dotimes (i maxiter)
	(when (>= (g alpha) (+ f0 (* alpha tolerance)))
	  (return-from linesearch-armijo alpha))
	(setf alpha (* alpha rho)))
      (error "reached maximum number of iterations, alpha=~a is still
not good enough" alpha))))




(defun linesearch (f x0 direction
		   &key 
		   (lower-tolerance 0.0001) (upper-tolerance 0.9999)
		   (lower-bound 0.1) (upper-bound 0.5)
		   (maxiter 50) f0 fd0)
  "return lambda, g(lambda)"
;;  (format t "entering linesearch, x0=~a, direction=~a~%" x0 direction)
  ;; check control parameters
  (assert (and (plusp lower-tolerance) (< lower-tolerance 0.5)
	       (< 0.5 upper-tolerance) (< upper-tolerance 1.0)
	       (plusp lower-bound) (<= lower-bound 0.5)
	       (<= 0.5 upper-bound) (< upper-bound 1.0)))
  ;; if initial values not given, calculate
  (unless (and f0 fd0)
    (setf (values f0 fd0) (funcall f x0 t)))
  (let* ((n (length x0))
	 (g0 f0)
	 (gd0 (dot-product direction fd0))
	 (lower-tol (* lower-tolerance gd0))
	 (upper-tol (* upper-tolerance gd0))
	 lambda
	 lambda-prev 
	 g-lambda
	 g-lambda-prev)
;;    (format t "lower-tol=~a  upper-tol~a~%" lower-tol upper-tol)
    (flet ((g (s)
	     "Calculate f(x0+s*direction), also return x0+s*direction."
	     (let ((xs (make-ffa n :double)))
	       (dotimes (i n)
		 (setf (aref xs i) (+ (aref x0 i) (* s (aref direction i)))))
	       (funcall f xs)))
	   (in-interval-p (x)
	     "Check if x is in the desired tolerance interval."
	     (and (<= lower-tol x) (<= x upper-tol))))
      ;; set lambda-prev to 1d0
      (setf lambda-prev 1d0
	    g-lambda-prev (g 1d0))
      ;; quadratic approximation
      (let ((a (- g-lambda-prev g0 gd0)))
	(setf lambda
	      (cond
		((zerop a)
		 (if (minusp gd0)
		     lower-bound	; b > 0, line slopes down
		     upper-bound))
		((minusp a)
		 (max lower-bound (/ gd0 (twice (- a)))))
		(t			; a > 0, parabola opens upwards!
		 lower-bound)))
	(setf g-lambda (g lambda))
;	(format t "will check if g(~a)=~a is ok~%" lambda g-lambda)
	(if (in-interval-p (/ (- g-lambda g0) lambda))
	  (return-from linesearch (values lambda g-lambda 0)))
;;	(format t "quadratic was not good enough, trying cubic~%")
	;; cubic approximation
	(dotimes (iteration maxiter)
	  (let* ((lambdadiff (- lambda lambda-prev))
		 (h1 (* (square lambda) lambdadiff))
		 (h2 (* (square lambda-prev) lambdadiff))
		 (k1 (/ (- g-lambda (* gd0 lambda) g0) h1))
		 (k2 (/ (- g-lambda-prev (* gd0 lambda-prev) g0) h2))
		 (a (- k1 k2))
		 (b (- (* lambda k2) (* lambda-prev k1))))
	    ;; save previous lambda values
	    (setf lambda-prev lambda
		  g-lambda-prev g-lambda)
	    ;; calculate new lambda and g(lambda)
	    (setf lambda (max 
			  (min 
			   (if (zerop a)
			       (/ gd0 b -2d0)
			       (let ((discriminant (- (square b) (* 3 a gd0))))
				 (when (minusp discriminant)
				   (error "negative discriminant: not
				 supposed to happen."))
				 (/ (+ b (sqrt discriminant)) a -3d0)))
			   (* upper-bound lambda-prev))
			  (* lower-bound lambda-prev)))
	    (setf g-lambda (g lambda))
;;	    (format t "******** linesearch: lambda=~a~%" lambda)
	    ;; check if we are in the desired interval
	    (if (in-interval-p (/ (- g-lambda g0) lambda))
		(return-from linesearch 
		  (values lambda g-lambda (1+ iteration)))))))))
  (error "problem with algorithm: not supposed to get here."))
