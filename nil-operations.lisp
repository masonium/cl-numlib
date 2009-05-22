(in-package :cl-numlib)

(defmacro define-nil-operation (name function &optional (docstring ""))
  `(defun ,name (&rest arguments)
     ,docstring
     (if (some #'null arguments)
	 nil
	 (apply #',function arguments))))

(define-nil-operation nil- - "- which returns nil if any of the arguments are nil")
