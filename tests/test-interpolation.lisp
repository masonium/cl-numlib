(in-package :cl-numlib)

(bind ((x #(0 1 2))
       (y1 #(3 4 6))
       (y2 #(7 5 4))
       ((y1-int y2-int) (make-interpolating-functions x y1 y2)))
  (defparameter *y1-int* y1-int)
  (defparameter *y2-int* y2-int))

(iter
  (for x :from -1 :to 2.5 :by 0.1)
  (for xx := (coerce (/ (round x 0.1) 10) 'single-float))
  (format t "x=~a  y1=~a  y2=~a~%" xx
	  (funcall *y1-int* xx)
	  (funcall *y2-int* xx)))
 
