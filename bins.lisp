(in-package :cl-numlib)

(defun interval-binary-search (boundaries x i j)
  "Return k such that boundaries[k-1] <= x < boundaries[k].  Search
starts assuming that boundaries[i] <= x <= boundaries[j], but this is
not checked.  Also, boundaries should be (weakly) increasing, but this
is not checked, the algorithm terminates anyway.  _If_ these
properties are satisfied, the elements will bracket x.

Note: this is the plainest interval binary search algorithm
imaginable, it is meant to be used by other functions which do tricks,
error checking, etc."
  (assert (and (<= 0 i) (< i j (length boundaries))))
  (tagbody
   top
     (when (= (1+ i) j)
       (return-from interval-binary-search j))
     (let ((m (floor (+ i j) 2)))
       (if (<= (aref boundaries m) x)
	   (setf i m)
	   (setf j m)))
     (go top)))

(defun make-binning-function-internal (boundaries)
  "Internal function for making binning function.  Does not check &
copy its arguments.  Make sure you know what you are doing."
  (let* ((n (length boundaries))
	 (last (1- n))
	 (left (aref boundaries 0))
	 (right (aref boundaries last))
	 (index 0))		; initial guess: below boundary
        (lambda (x)
      (flet ((find-index (i j)
	       (interval-binary-search boundaries x i j)))
;;	(format t "index=~a x=~a~%" index x)
	(cond 
	 ((zerop index)
	   ;; guess: x < left
	   (when (<= left x)
	     ;; no, boundaries[0] <= x
	     (setf index
		   (if (<= right x)
		       ;; right <= x
		       n
		       ;; we are inside
		       (find-index 0 last)))))
	 ((= index n)
	  ;; guess: right <= x
	  (unless (<= right x)
	    ;; no, x < right
	    (setf index
		  (if (<= left x)
		      ;; we are inside
		      (find-index 0 last)
		      ;; x < left
		      0))))
	 ((<= 1 index last)
	  ;; guess: we are inside
	  (unless (and (<= (aref boundaries (1- index)) x)
		       (< x (aref boundaries index)))
	    ;; no, did not find it in previous bracket
	    (setf index
		  (cond
		    ((< x left) 0)
		    ((<= right x) n)
		    (t (find-index 0 last))))))
	  (t (error "internal error, should not happen")))
	index))))

     
(defun make-binning-function (boundaries)
  "Return a function that places its argument in the appropriate bin.
The integer i returned by the function has the following interpretation:

  0 for number < boundaries[0]

  1 <= i < n for boundaries[i-1] <= number < boundaries[i]

  n for boundaries[n-1] <= number 

where n is the length of the vector (or list) `boundaries', which
needs to be strictly increasing.

Examples:

 (funcall (make-binning-function #(1 2 3)) 0.5) ; => 0
 (funcall (make-binning-function #(1 2 3)) 1)   ; => 1
 (funcall (make-binning-function #(1 2 3)) 2.5) ; => 2
 (funcall (make-binning-function #(1 2 3)) 4.5) ; => 3

Notes: the last index is cached in a closure and is always checked
first.  `boundaries' is copied, so it does not share structure with
the original."
  (let* ((boundaries (cond
		       ((vectorp boundaries) (array-copy boundaries))
		       ((listp boundaries) (coerce boundaries 
						   '(simple-array * (*))))
		       (t (error "boundaries should be a vector or a list")))))
    ;; check that boundaries is increasing
    ;; !!!! could make this better using iter's facilities
    (iter
      (for i :from 1 :below (length boundaries))
      (let ((current (aref boundaries i))
	    (previous (aref boundaries (1- i))))
	(unless (< previous current)
	  (error "Nonincreasing boundaries ~a >= ~a at index ~a" 
		 previous current i))))
    ;; finder function
    (make-binning-function-internal boundaries)))

(defun count-in-bins (data &key 
		      (key #'identity) (lowest nil) (highest nil)
		      (weight-key nil)
		      (weight-type (if weight-key
				       (type-of (funcall weight-key (aref data 0)))
				       'fixnum)))
  "Call `key' (determines the bin) and `weight-key' (determines the
weight) on each element of the vector `data', and sum the weights for
each bin.  When `weight-key' is the default, just count elements in
each bin.

Return two values, the value of the smallest integer, and a vector of
counts/weights starting from there.

Example: 

 (count-in-bins #(0 1 2 3 2 -1)) ; => (values -1 #(1 1 1 2 1))

Note: the vector of counts is extended on demand as elements come, but
you can provide a starting value for the boundaries.  Of course, there
is no guarantee that these will not be extended on demand."
  (assert (and (vectorp data) (plusp (length data))))
  (bind (((:values lowest highest)
	  (cond 
	    ((and lowest highest) (assert (and (integerp lowest) 
					       (integerp highest)
					       (<= lowest highest)))
	     (values lowest highest))
	    (lowest (assert (integerp lowest))
		    (values lowest lowest))
	    (highest (assert (integerp highest))
		     (values highest highest))
	    (t (let ((first-bin (funcall key (aref data 0))))
		 (values first-bin first-bin)))))
	 (zero (coerce 0 weight-type))
	 (bin-counts (make-array (1+ (- highest lowest))
		      :element-type weight-type 
		      :initial-element zero)))
    (dotimes (i (length data))
      (let* ((element (aref data i))
	     (bin (funcall key element)))
	(assert (integerp bin))
	(cond 
	  ((< bin lowest)
	   ;; extend to the left
	   (let ((shift (- lowest bin)))
	     (setf bin-counts (extend-vector bin-counts (- shift) zero))
	     (decf lowest shift)))
	  ((> bin highest)
	   (let ((shift (- bin highest)))
	     (setf bin-counts (extend-vector bin-counts shift zero))
	     (incf highest shift))))
	(incf (aref bin-counts (- bin lowest)) (if weight-key 
						   (funcall weight-key element)
						   1))))
    (values lowest bin-counts)))

(defun histogram-pretty-breaks (data key n &optional (extend-right-rel 1.00001))
  "Prettified brakes for data, using n breakpoints (if not given,
computed using Sturges' rule).  Return two values, the breakpoints and
an exponent, see pretty-axis-positions."
  (bind ((n (if n n (1+ (ceiling (log (length data) 2)))))
	 ((lower upper) (array-range data :key key)))
    (pretty-axis-positions lower (* upper extend-right-rel) 
			   n :include-lower-p t :include-upper-p t)))

(defun histogram-automatic-breaks (&optional n)
  "Histogram break function for automatic breaks.  If n is given, at
most this many breaks are generated, putting them at `round'
numbers (see pretty-axis-positions).  If n is not given, Sturges' rule
is used.  Return breaks."
  (lambda (data key)
    (histogram-pretty-breaks data key n)))

(defun histogram-evenly-distributed-breaks (n &optional (extend-right-rel 1.00001))
  "Histogram break function that returns n evenly distributed breaks
over the range of data.  The parameter extend-right-rel will extend
the upper bound so that the maximum still falls in one of the bins."
  (assert (< 1 extend-right-rel))
  (lambda (data key)
    (bind (((lower upper) (array-range data :key key)))
      (num-sequence :from lower :to (* extend-right-rel upper) :length (1+ n)))))

(defun histogram-fixed-breaks (breaks)
  "Histogram break function that just returns the prespecified
breaks (no copy is made of the vector.  Bounds are checked."
  (lambda (data key)
    (bind (((lower upper) (array-range data :key key))
	   ((b-lower b-upper) (array-range breaks)))
      (unless (and (<= b-lower lower) (>= b-upper upper))
	(error "data will not fit inside breaks provided."))
      breaks)))

(defclass histogram ()
  ((breaks :initarg :breaks :documentation "breaks for the histogram"
	   :accessor breaks)
   (counts :initarg :counts :documentation 
	    "counts/weights for each interval, one less than breaks"
	    :accessor counts)
   (exponent :initarg :exponent :documentation
	     "number of digits in breaks, used for formatting, ignore
	     if nil")))

(defmethod print-object ((obj histogram) stream)
  ;;;; !!! todo: variable max column (currently 60), variable width for break edges
  (print-unreadable-object (obj stream :type t)
    (terpri stream)
    (with-slots (breaks counts) obj
      (dotimes (i (length counts))
	(format stream "[~5f,~5f) " (aref breaks i) (aref breaks (1+ i)))
	(let* ((count (aref counts i))
	       (rounded-count (round count)))
	  (if (<= rounded-count 60)
	      (dotimes (j rounded-count)
		(format stream "*"))
	      (format stream "*****..... (~a)" count)))
	(terpri stream)))))

(defun histogram (data
		  &key (breaks-function (histogram-automatic-breaks))
		  (key #'identity) (weight-key nil))
  "Return a histogram of `data'.  The function breaks-function will be
passed `data' and `key', and is expected to return a vector of breaks.
If it is a positive integer, that many prettified bins are created.
For `key' and `weight-key', see count-in-bins."
  (bind (((:values breaks exponent) (funcall breaks-function data key))
	 (binning-function (make-binning-function breaks))
	 (n-breaks (length breaks))
	 ((:values lowest bin-counts)
	  (count-in-bins data :key (lambda (v)
				     (let ((bin (funcall binning-function 
							 (funcall key v))))
				       (unless (< 0 bin n-breaks)
					 (error "value ~a is outside breaks" v))
				       bin))
	   :weight-key weight-key
	   :lowest 0
	   :highest (1- n-breaks))))
    ;; should never fail
    (assert (and (zerop lowest) (= (length bin-counts) n-breaks)))
    ;; return histogram object
    (make-instance 'histogram :breaks breaks 
		   ;; shift
		   :counts (displace-array bin-counts (1- n-breaks) 1)
		   :exponent exponent)))

;; (defparameter *r* (replicate (lambda () (exponential-random 1d0)) 400 'double-float))

;; (defparameter *hist* (histogram *r*))



