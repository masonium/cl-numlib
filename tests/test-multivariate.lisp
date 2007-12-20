(in-package :cl-numlib)

;;; line search
(defun testfun-quadratic (x &optional derivative-p)
  (declare (ignore derivative-p))
  (let ((x (aref x 0)))
    (format t "~a~%" x)
    (values (square x) (make-ffa 1 :double :initial-contents (list (* 2 x))))))

(defun testfun-sin (x &optional derivative-p)
  (declare (ignore derivative-p))
  (let ((x (aref x 0)))
    (format t "~a~%" x)
    (values (sin x) (make-ffa 1 :double :initial-contents (list (cos x))))))

(defun testfun-kink (x &optional derivative-p)
  (declare (ignore derivative-p))
  (let ((x (aref x 0)))
    (format t "~a~%" x)
    (if (minusp x)
	(values x (make-ffa 1 :double :initial-element -1))
	(values (1- (/ 1 (1+ x))) (make-ffa 1 :double :initial-element
					    (/ -1 (square (1+ x))))))))
(defun testfun-kink2 (x &optional derivative-p)
  (declare (ignore derivative-p))
  (let ((x (aref x 0)))
    (format t "~a~%" x)
    (if (minusp x)
	(values x (make-ffa 1 :double :initial-element -1))
	(values (- (square x)) (make-ffa 1 :double :initial-element (* -2 x))))))
  


(linesearch #'testfun-quadratic #(-1) #(1))
(linesearch #'testfun-sin #(-1) #(-1) :alpha 0.5d0)
(linesearch #'testfun-kink #(-1d0) #(0.7d0) :alpha 1)
(linesearch #'testfun-kink2 #(-1d0) #(0.7d0) :alpha 1)

;;;;

(defun quadratic-0 (x &optional derivative-p)
  "Quadratic around 0."
  (let ((x1 (aref x 0))
	(x2 (aref x 1))
	(fx (l2-norm x)))
  (format t "x=(~a,~a)~%" x1 x2)
  (if derivative-p
      (values fx (vector (* 2d0 x1) (* 2d0 x2)))
      fx)))

(defun banana-function (a &optional derivative-p)
  (declare (ignore derivative-p))
  "\"Banana\" function f(x,y) = -100(y-x^2)^2+(1-x)^2."
  (let ((x (aref a 0))
	(y (aref a 1)))
    (format t "x=(~a,~a)~%" x y)
    (values (- (square (- 1 x))
	       (* -100 (square (- y (square x)))))
	    (vector (+ (* -2 (- 1 x))
		       (* -400 (- y (square x)) x))
		    (* 200 (- y (square x)))))))
  
(bfgs-minimize #'banana-function #(0 0) :H-scaling 4d0)
(bfgs-minimize #'banana-function #(.9 .9))
(bfgs-minimize #'quadratic-0 #(1 9))
