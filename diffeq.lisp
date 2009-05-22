(in-package :cl-numlib)

(defun rk4 (f y0 start end &key (steps 1000))
  "Runge-Kutta 4th order solver.  We are solving dy/dx = f(x,y) on x
\in [start,end], with y(start)=y0, using a given number of steps.  f
accepts x as a double-float and y as a vector of double-floats, and
returns dy/dx as a vector of double-floats.

If start > end, solution will proceed backwards (note that x will be
in descending order, with y and yp ordered accordingly).

Returns the vector x and matrices y and yp (each coordinate in a single ROW),
as (values x y yp).

All matrices and vectors need to be double-float.

The f should have arguments (x y i), where x is a real number, y is a
vector of real numbers, and i, when non-nil, is an integer which gives
the index of the current step (so you can save intermediate
calculations using a closure if desired).

The algorithm used is plain vanilla 4th order Runge-Kutta, without
error checking or adaptive stepsize."
  (assert (typep y0 '(array * (*))))
  (let* ((n (length y0))
	 (start (coerce start 'double-float))
	 (end (coerce end 'double-float))
	 (h (/ (- end start) (1- steps)))
	 (h2 (/ h 2))
	 (h3 (/ h 3))
	 (h6 (/ h 6))
	 (x (make-array steps :element-type 'double-float))
	 (y (make-array (list n steps) :element-type 'double-float))
	 (yp (make-array (list n steps) :element-type 'double-float))
	 (y# (array-copy y0)))
    (iter
      (for i :from 0 :below steps)
      (for x# := (+ start (* h i)))
      (for k1 := (funcall f x# y# i))
      (for k2 := (funcall f (+ x# (/ h 2)) 
			  (vectorize (y# k1) (+ y# (* h2 k1)))
			  nil))
      (for k3 := (funcall f (+ x# (/ h 2)) 
			  (vectorize (y# k2) (+ y# (* h2 k2)))
			  nil))
      (for k4 := (funcall f (+ x# h)
			  (vectorize (y# k3) (+ y# (* h k3)))
			  nil))
      ;; save results
      (setf (aref x i) x#)
      (dotimes (j n)
	(setf (aref y j i) (aref y# j)
	      (aref yp j i) (aref k1 j)))
      ;; update y#
      (vectorize (y# k1 k2 k3 k4)
		 (+ y# (* k1 h6) (* k2 h3) (* k3 h3) (* k4 h6))
		 :into y#))
    ;; return results
    (values x y yp)))
