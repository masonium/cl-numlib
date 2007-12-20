(in-package :cl-numlib)

(defun root-ridders (f a b &optional
		     (tolerance (* (abs (- b a)) #.(expt double-float-epsilon 0.25)))
		     (epsilon #.(expt double-float-epsilon 0.25)))
  "Find the root of f bracketed between a and b using Ridders' method.
The algorithm stops when either the root is bracketed in an interval
of length 2*tolerance, or root is found such that abs(f(root)) <=
epsilon.

Return three values: the root, the function evaluated at the root, and
a boolean which is true iff abs(f(root)) <= epsilon."
;;   (declare (double-float a b tolerance epsilon)
;; 	   ((function (double-float) double-float) f))
  (flet ((opposite-sign (x y)
	   (minusp (* x y)))
	 (evaluate-and-check-zero (x)
	   (let ((fx (funcall f x)))
	     (if (<= (abs fx) epsilon)
		 (return-from root-ridders (values x fx t))
		 fx))))
    (macrolet ((new-bracket (a b fa fb)
		 `(progn
		    (setf a ,a
			  b ,b
			  fa ,fa
			  fb ,fb)
		    (go top))))
      (let ((fa (evaluate-and-check-zero a))
	    (fb (evaluate-and-check-zero b)))
	(when (< b a)
	  (rotatef a b))
	(unless (opposite-sign fa fb)
	  (error "f(~a)=~a and f(~a)=~a don't bracket a root." a fa b fb))
	(tagbody
	 top
;;	   (format t "~a ~a~%" a b)
	   (let* ((d (/ (- b a) 2))	; half-distance
		  (m (+ a d))		; midpoint
		  (fm (evaluate-and-check-zero m))) ; value at midpoint
	     (when (<= d tolerance)
	       (return-from root-ridders (values m fm nil)))
	     (let* ((w (- (square fm) (* fa fb)))	   ; discriminant
		    (delta (/ (* (signum fa) fm d) (sqrt w))) ; c-m
		    (c (+ m delta))	; interpolated guess
		    (fc (evaluate-and-check-zero c))) ; value at guess
	       (if (minusp delta)
		   ;; c < m
		   (cond
		     ((opposite-sign fm fc) (new-bracket c m fc fm))
		     ((opposite-sign fa fc) (new-bracket a c fa fc))
		     ((opposite-sign fb fc) (new-bracket c b fc fb))
		     (t (error "internal error")))
		   ;; m < c
		   (cond
		     ((opposite-sign fm fc) (new-bracket m c fm fc))
		     ((opposite-sign fb fc) (new-bracket c b fc fb))
		     ((opposite-sign fa fc) (new-bracket a c fa fc))
		     (t (error "internal error")))))))))))
