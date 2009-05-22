(defpackage #:cl-numlib-asd
  (:use :cl :asdf))

(in-package :cl-numlib-asd)

(defsystem :cl-numlib
 :description "Numerical functions (optimization, roofinding, utilities)"
 :version "0.1"
 :author "Tamas K Papp"
 :license "GPL"
 :components ((:file "package")
	      (:file "optimization" :depends-on ("package"))
	      (:file "utilities" :depends-on ("package"))
	      (:file "rootfinding" :depends-on ("utilities"))
	      (:file "integration" :depends-on ("package"))
	      (:file "functions" :depends-on ("utilities"))
	      (:file "multivariate-optimization" :depends-on ("utilities"))
	      (:file "bins" :depends-on ("utilities"))
	      (:file "diffeq" :depends-on ("package"))
	      (:file "interpolation" :depends-on ("bins"))
	      (:file "random" :depends-on ("utilities" "bins"))
	      (:file "statistics" :depends-on ("utilities" "bins")))
 :depends-on (:cl-utilities :iterate :ffa :array-operations :metabang-bind
			    #| :row-major-indexing |#))
