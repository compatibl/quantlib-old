; Copyright (C) 2002, 2003 RiskMap srl
;
; This file is part of QuantLib, a free-software/open-source library
; for financial quantitative analysts and developers - http://quantlib.org/
;
; QuantLib is free software developed by the QuantLib Group; you can
; redistribute it and/or modify it under the terms of the QuantLib License;
; either version 1.0, or (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; QuantLib License for more details.
;
; You should have received a copy of the QuantLib License along with this
; program; if not, please email quantlib-dev@lists.sf.net
;
; The QuantLib License is also available at http://quantlib.org/license.html
; The members of the QuantLib Group are listed in the QuantLib License

(load "unittest.scm")
(use-modules (QuantLib))

(greg-assert/display
 "Testing observability of market elements"
 (let ((flag #f))
   (deleting-let ((me (new-SimpleMarketElement 0.0)
                      delete-MarketElement)
                  (obs (new-Observer (lambda () (set! flag #t)))
                       delete-Observer))
     (deleting-let ((temp (MarketElement->Observable me)
                          delete-Observable))
       (Observer-register-with obs temp))
     (SimpleMarketElement-value-set! me 3.14)
     (check flag
            "Observer was not notified of market element change"))))

(greg-assert/display
 "Testing observability of market element handles"
 (let ((flag #f))
   (deleting-let* ((me1 (new-SimpleMarketElement 0.0)
                        delete-MarketElement)
                   (me2 (new-SimpleMarketElement 0.0)
                        delete-MarketElement)
                   (h (new-MarketElementHandle me1)
                      delete-MarketElementHandle)
                   (obs (new-Observer (lambda () (set! flag #t)))
                        delete-Observer))
     (deleting-let ((temp (MarketElementHandle->Observable h)
                          delete-Observable))
       (Observer-register-with obs temp))
     (and
      (begin (SimpleMarketElement-value-set! me1 3.14)
             (check flag
                    "Observer was not notified of market element change"))
      (begin (set! flag #f)
             (MarketElementHandle-link-to! h me2)
             (check flag
                    "Observer was not notified of market element change"))))))

