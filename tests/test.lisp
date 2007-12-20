(in-package :cl-numlib)

(root-ridders #'sin -0.3d0 2.9d0)

(root-ridders #'(lambda (x) (if (< x 3) -1 1)) -5 9)
