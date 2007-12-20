;; golden-section-minimize
(golden-section-minimize (lambda (x) (let ((d (- x 2)))
				       (the double-float (- (* d d)))))
			 -5d0 12d0 1d-5) ; should be around 12

(golden-section-minimize (lambda (x) (let ((d (- x 2)))
				       (the double-float (* d d))))
			 -5d0 12d0 1d-5) ; should be around 2

(golden-section-minimize #'sin
			 pi (* 2 pi) 1d-5) ; should be around (* 1.5 pi) \approx 4.71

(golden-section-minimize #'sin 0d0 pi 1d-5) ; should be around pi

;; newton-minimize
(newton-minimize #'(lambda (x) (values (cos x) (- (sin x)))) -1 -10 -1)
