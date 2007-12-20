(in-package :cl-numlib)

(defmacro evaluate-polynomial (x &rest coefficients)
  "Expand the evaluation of coefficients as a polynomial (coefficients
increasing by degree, the constant first) using Horner's method.  The
variable is evaluated once only."
  (once-only (x)
    (labels ((horner-method (coefficients)
	       (cond 
		 ((null coefficients) 0)
		 ((null (rest coefficients)) (first coefficients))
		 (t `(+ ,(first coefficients) 
			(* ,x ,(horner-method (rest coefficients))))))))
      (horner-method coefficients))))

(defun standard-normal-quantile (p)
  "Calculate quantile of a standard normal distribution,
using the approximation of Peter J Acklam,
http://home.online.no/~pjacklam/notes/invnorm/.  Relative error has an
absolute value of less than 1.15e-9 over the entire region."
  (unless (and (plusp p) (< p 1))
    (error "p has to lie in (0,1)"))
  (flet ((a (r)
	   (evaluate-polynomial r
				2.506628277459239d+00 -3.066479806614716d+01
				1.383577518672690d+02 -2.759285104469687d+02
				2.209460984245205d+02 -3.969683028665376d+01))
	 (b (r)
	   (evaluate-polynomial r 1d0
				-1.328068155288572d+01 6.680131188771972d+01
				-1.556989798598866d+02 1.615858368580409d+02
				-5.447609879822406d+01))
	 (c (q)
	   (evaluate-polynomial q
				2.938163982698783d+00 4.374664141464968d+00
				-2.549732539343734d+00 -2.400758277161838d+00
				-3.223964580411365d-01 -7.784894002430293d-03))
	 (d (q)
	   (evaluate-polynomial q 1d0
				3.754408661907416d+00 2.445134137142996d+00
				3.224671290700398d-01 7.784695709041462d-03)))
    (let ((p (coerce p 'double-float)))
      (cond
	((<= p 0.02425)
	 (let ((q (sqrt (* -2d0 (log p)))))
	   (/ (c q) (d q))))
	((>= p #.(- 1 0.02425))
	 (let ((q (sqrt (* -2d0 (log (- 1 p))))))
	   (- (/ (c q) (d q)))))
	(t (let* ((q (- p 0.5d0))
		  (r (square q)))
	     (/ (* (a r) q) (b r))))))))


(defun normal-quantile (p mean standard-deviation)
  "Quantile function for normal distribution."
  (+ (* (standard-normal-quantile p) standard-deviation) mean))
