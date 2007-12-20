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
	      (:file "random" :depends-on ("utilities"))
	      (:file "rootfinding" :depends-on ("utilities"))
	      (:file "integration" :depends-on ("package"))
	      (:file "functions" :depends-on ("utilities"))
	      (:file "multivariate-optimization" :depends-on ("utilities")))
 :depends-on (:cl-utilities :iterate :ffa :metabang-bind))