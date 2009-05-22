
(defun test-binning-function (bins N)
  "Test make-binning-function `N' times with `bins' bins."
  (let ((boundaries (make-ffa (1+ bins) :double)))
    (dotimes (i (length boundaries))
      (setf (aref boundaries i) (random 1000d0)))
    (sort boundaries #'<)
    (format t "boundaries: ~a~%" boundaries)
    (let ((binning-function (make-binning-function boundaries)))
      (dotimes (n N)
	(let* ((x (- (random 1100d0) 50d0))
	       (bin (funcall binning-function x)))
;;	  (format t "x=~a bin=~a~%" x bin)
	  (cond
	    ((zerop bin) (assert (< x (aref boundaries 0))))
	    ((= bin (1+ bins)) (assert (<= (aref boundaries bins) x)))
	    (t (assert (and (<= (aref boundaries (1- bin)) x)
			    (< x (aref boundaries bin)))))))))))

(test-binning-function 40 10000)

(funcall (make-binning-function #(1 2 3)) 0.5) ; => 0
(funcall (make-binning-function #(1 2 3)) 1)   ; => 1
(funcall (make-binning-function #(1 2 3)) 2.5) ; => 2
(funcall (make-binning-function #(1 2 3)) 4.5) ; => 3
