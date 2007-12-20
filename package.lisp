(defpackage :cl-numlib
  (:use :common-lisp :cl-utilities :iterate :ffa :bind)
  (:shadowing-import-from :iterate :collecting :collect)
  (:export

   ;; optimization
   
   golden-section-minimize univariate-outside-boundaries
   reached-maximum-iterations univariate-found-local-maximum
   newton-minimize negate-univariate-function negate-univariate-fp

   ;; utilities

   int-sequence num-sequence maxf minf multf square positive-part
   negative-part twice half mean arrays-range matrix-multiply
   matrix-vector-multiply identity-matrix sum-of-squares l2-norm
   find-point-approaching-boundary quadratic-roots

   ;; random

   weighted-draw round-probabilities-to-count fill-with-weighted-draw
   replicate

   ;; rootfinding

   root-ridders

   ;; integration

   integrate

   ;; functions
   
   standard-normal-quantile normal-quantile

   ;; multivariate-optimization
   
   bfgs-minimize negate-multivariate-objective-function

   ))
