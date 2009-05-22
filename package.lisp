(in-package :cl-numlib-asd)

(defpackage :cl-numlib
  (:use :common-lisp :cl-utilities :iterate :ffa :array-operations :bind
	#|:row-major-indexing|#)
  (:shadowing-import-from :iterate :collecting :collect)
  (:export

   ;; optimization
   
   golden-section-minimize univariate-outside-boundaries
   reached-maximum-iterations univariate-found-local-maximum
   newton-minimize negate-univariate-function negate-univariate-fp

   ;; utilities

   int-sequence num-sequence maxf minf maxf* minf* multf square
   positive-part negative-part twice half mean absolute-difference
   arrays-range matrix-transpose matrix-multiply
   matrix-vector-multiply quadratic-form identity-matrix
   sum-of-squares l1-norm l2-norm sup-norm
   find-point-approaching-boundary quadratic-roots solve-2x2-system
   iter-indexing-vector -+ -* displace-matrix-row slice-matrix-rows
   normalize-vector normalize-probability first-difference
   #|array-finite-difference |# extend-vector pretty-axis-positions
   vector-last make-gamma-transform

   ;; random

   weighted-draw round-probabilities-to-count fill-with-weighted-draw
   replicate make-empirical-inverse-random-sampler exponential-random
   first-exponential-random make-poisson-random

   ;; rootfinding

   boundaries-dont-bracket-root root-ridders root-bisection
   find-satisfactory-fx make-expanding-rule make-contracting-rule
   root-autobracket

   ;; integration

   integrate simpsons-rule-on-index

   ;; functions
   
   standard-normal-quantile normal-quantile

   ;; multivariate-optimization
   
   linesearch-reached-maximum-iterations-error bfgs-stuck-error
   bfgs-H-not-positive-definite-error bfgs-minimize
   negate-multivariate-objective-function
   
   ;; bins
   
   make-binning-function-internal make-binning-function count-in-bins
   histogram-evenly-distributed-breaks histogram-automatic-breaks
   histogram-fixed-breaks histogram breaks counts exponent

   ;; diffeq

   rk4

   ;; interpolation

   make-interpolating-functions make-interpolating-function

   ;; statistics
   
   quantiles quantiles-weighted mean-variance-weighted

   ))
