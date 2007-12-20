(in-package :cl-numlib)

;;;;
;;;;  primitive integration functions, operate on [0,1]
;;;;
    
(defun trapezoidal-closed (f eps &key (minstep 5) (maxstep 10))
  (let ((sum (half (+ (funcall f 0) (funcall f 1))))
	(h 1)
	(n 1)
	(step 0))
    (tagbody
       top
       (when (> step maxstep)
	 (error "exceeded maxstep"))
       (let* ((s (iter
		  (for i :from 0 :below n)
		  (for x :from (half h) :by h)
		  (summing (funcall f x))))
	     (new-sum (half (+ sum (* s h)))))
	 (when (and (>= step minstep)
		    (< (abs (- sum new-sum)) eps))
	   (return-from trapezoidal-closed new-sum))
;;;;	 (format t "sum: ~a s: ~a~%" sum s)
	 (setf n (twice n)
	       h (half h)
	       sum new-sum))
       (incf step)
       (go top))))

;;;;
;;;; transformations to map from [0,1] to domains, incorporating derivative
;;;;

(defun map-a-b (f a b)
  (let* ((multiplier (- b a))
	 (derivative (abs multiplier)))
    (lambda (x)
      (* (funcall f (+ (* multiplier x) a)) derivative))))

(defun map-a-inf (f a)
  "Map [0,1) to [a,inf] using the mapping x/(1-x)+a, also multiply by
the derivative."
  (lambda (x)
    (if (= x 1)
	0
	(let* ((1minusx (- 1 x)))
	  (* (funcall f (+ (/ x 1minusx) a)) (/ (square 1minusx)))))))

(defun map-inf-b (f b)
  "Map [0,1) to [-inf,b] using the mapping 1-1/x+b, also multiply by
the derivative."
  (lambda (x)
    (if (zerop x)
	0
	(* (funcall f (- (+ 1 b) (/ 1 x))) (/ (square x))))))

(defun map-inf-inf (f)
  (declare (ignore f))
  (error "need to write this"))
	
(defun integrate (f a b &key 
		  (method 'trapezoidal-closed)
		  (eps 1e-5)
		  (minstep 5)
		  (maxstep 15))
  (declare (ignore method))
  (cond
    ((and (eq a :minusinfinity) (eq b :plusinfinity))
     (trapezoidal-closed (map-inf-inf f) eps :minstep minstep :maxstep maxstep))
    ((eq a :minusinfinity)
     (trapezoidal-closed (map-inf-b f b) eps :minstep minstep :maxstep maxstep))
    ((eq b :plusinfinity)
     (trapezoidal-closed (map-a-inf f a) eps :minstep minstep :maxstep maxstep))
    (t
     (trapezoidal-closed (map-a-b f a b) eps :minstep minstep :maxstep maxstep))))
