(in-package :cl-numlib)

(require :cl-2d)


(defun f (x y)
  (declare (ignore x))
  (vector (* -1.4d0 (aref y 0))))

(bind ((y0 9d0)
       (steps 1000)
       (end 10d0)
       ;; forward
       ((:values x y yp) (rk4 #'f (vector y0) end :steps steps))
       (y (displace-array y steps 0))
       (yp (displace-array yp steps 0))
       (truey (vectorize (x) (* y0 (exp (* -1.4d0 x)))))
       (diff (array- truey y))
       ;; backward
       ((:values xb yb ypb) (rk4 #'f (vector (* y0 (exp (* -1.4d0 end))))
			     0d0 :start 10d0 :steps steps))
       (xb (reverse (flatten-array xb)))
       (yb (reverse (flatten-array yb)))
       (ybp (reverse (flatten-array ypb)))
       (trueyb (vectorize (xb) (* 9d0 (exp (* -1.4d0 xb)))))
       (diffb (array- trueyb y)))
  (defparameter *x* x)
  (defparameter *y* y)
  (defparameter *yp* yp)
  (defparameter *truey* truey)
  (defparameter *diff* diff)
  (defparameter *yb* yb)
  (defparameter *ypb* ypb)
  (defparameter *trueyb* trueyb)
  (defparameter *diffb* diffb))



(defparameter *frame* (cl-2d:context-as-frame
		       (cl-cairo2:create-xlib-image-context 400 400)))
(cl-2d:clear *frame*)
(cl-2d:plot-xy *frame* *x* *y*)
(cl-2d:plot-xy *frame* *x* *truey*)
(cl-2d:plot-xy *frame* *x* *diff*)
