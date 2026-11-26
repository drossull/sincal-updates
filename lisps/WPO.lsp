;;; WPO: envolvente convexa de la seleccion, sin modificar los originales.
;;; WIPEOUT es poligonal: arcos/elipses se discretizan cada 2 grados como maximo.
(vl-load-com)

(defun WPO:Fail (message)
  (setq WPO:reason message)
  (exit))

(defun WPO:Ceil (value / integer)
  (setq integer (fix value))
  (if (> value integer) (1+ integer) integer))

(defun WPO:Points (en / data kind flags normal bulges item vertex steps index
                          part count start end points point elevation)
  (setq data (entget en) kind (cdr (assoc 0 data))
        normal (cond ((cdr (assoc 210 data))) ('(0.0 0.0 1.0))))
  (if (not (member kind '("LWPOLYLINE" "POLYLINE" "CIRCLE" "ELLIPSE")))
    (WPO:Fail "Seleccione una polilinea 2D cerrada, circulo o elipse; no bloques ni curvas abiertas."))
  (if (= kind "POLYLINE")
    (progn
      (setq flags (cdr (assoc 70 data)))
      (if (/= 0 (logand flags 94))
        (WPO:Fail "No se admiten polilineas 3D, mallas ni polilineas ajustadas con Fit/Spline."))))
  (if (not (vlax-curve-isClosed en))
    (WPO:Fail "El contorno esta abierto. Cierre o una su geometria antes de ejecutar WPO."))
  (cond
    ((= kind "LWPOLYLINE")
      (foreach item data
        (cond ((= (car item) 10) (setq bulges (cons 0.0 bulges)))
              ((= (car item) 42) (setq bulges (cons (cdr item) (cdr bulges))))))
      (setq bulges (reverse bulges)))
    ((= kind "POLYLINE")
      (setq vertex (entnext en))
      (while (= (cdr (assoc 0 (entget vertex))) "VERTEX")
        (setq bulges (cons (cond ((cdr (assoc 42 (entget vertex)))) (0.0)) bulges)
              vertex (entnext vertex)))
      (setq bulges (reverse bulges))))
  (setq points nil count 0 index 0 start (vlax-curve-getStartParam en)
        end (vlax-curve-getEndParam en))
  ;; Preserve every straight polyline corner, including concave corners.
  (if bulges
    (foreach item bulges
      (setq steps (max 1 (WPO:Ceil (/ (abs (* 4.0 (atan item))) (/ pi 90.0))))
            part 0 count (+ count steps))
      (if (> count 2048) (WPO:Fail "Contorno demasiado complejo (mas de 2048 vertices). Simplifiquelo primero."))
      (repeat steps
        (setq points (cons (trans (vlax-curve-getPointAtParam en (+ index (/ (float part) steps))) 0 normal) points)
              part (1+ part)))
      (setq index (1+ index)))
    (progn
      (setq steps 180)
      (repeat steps
        (setq points (cons (trans (vlax-curve-getPointAtParam en (+ start (* (- end start) (/ (float index) steps)))) 0 normal) points)
              index (1+ index)))))
  (setq points (reverse points) elevation (caddr (car points)))
  (list (mapcar '(lambda (p) (list (car p) (cadr p))) points) normal elevation))

(defun WPO:Cross (a b c)
  (- (* (- (car b) (car a)) (- (cadr c) (cadr a)))
     (* (- (cadr b) (cadr a)) (- (car c) (car a)))))

(defun WPO:OnSegment (p a b eps)
  (and (<= (abs (WPO:Cross a b p)) (* eps (distance a b)))
       (<= (- (min (car a) (car b)) eps) (car p) (+ (max (car a) (car b)) eps))
       (<= (- (min (cadr a) (cadr b)) eps) (cadr p) (+ (max (cadr a) (cadr b)) eps))))

(defun WPO:Validate (points checkCross / origin area clean p eps span edges tail a b c d)
  (setq span (max (- (apply 'max (mapcar 'car points)) (apply 'min (mapcar 'car points)))
                  (- (apply 'max (mapcar 'cadr points)) (apply 'min (mapcar 'cadr points))))
        eps (max 1e-9 (* span 1e-10)))
  (foreach p points
    (if (or (null clean) (> (distance p (car clean)) eps)) (setq clean (cons p clean))))
  (if (and clean (equal (car clean) (car (last clean)) eps)) (setq clean (cdr clean)))
  (setq points (reverse clean))
  (if (< (length points) 3) (WPO:Fail "El contorno no encierra un area valida."))
  (setq origin (car points) area 0.0 edges (append points (list origin)) tail edges)
  (while (cadr tail)
    (setq area (+ area (WPO:Cross origin (car tail) (cadr tail))) tail (cdr tail)))
  (if (<= (abs area) (* eps span)) (WPO:Fail "El contorno es degenerado o se cruza sobre si mismo."))
  (if checkCross
    (while (cadr edges)
      (setq a (car edges) b (cadr edges) tail (cddr edges))
      (while (cadr tail)
        (setq c (car tail) d (cadr tail))
        ;; First and last edges share the closing vertex, which is allowed.
        (if (not (equal a d eps))
          (if (or (inters a b c d T) (WPO:OnSegment a c d eps)
                  (WPO:OnSegment b c d eps) (WPO:OnSegment c a b eps)
                  (WPO:OnSegment d a b eps))
            (WPO:Fail "El contorno se cruza o se toca sobre si mismo. Corrijalo antes de crear la mascara.")))
        (setq tail (cdr tail)))
      (setq edges (cdr edges))))
  points)

(defun WPO:Create (en / WPO:reason result data shape points temp mask previous)
  ;; A failure removes only the entities created here, never the source.
  (setq result (vl-catch-all-apply
    '(lambda ()
      (setq data (entget en))
      (if (/= (cdr (assoc 410 data)) (if (= (getvar "CVPORT") 1) (getvar "CTAB") "Model"))
        (WPO:Fail "Seleccione un objeto del espacio actual."))
      (if (/= 0 (logand 4 (cdr (assoc 70 (tblsearch "LAYER" (getvar "CLAYER"))))))
        (WPO:Fail "La capa actual esta bloqueada. Elija una capa desbloqueada para la mascara."))
      (if (/= 0 (logand 4 (cdr (assoc 70 (tblsearch "LAYER" (cdr (assoc 8 data)))))))
        (WPO:Fail "La capa del contorno esta bloqueada; desbloqueela para ordenar la mascara debajo."))
      (setq shape (WPO:Points en)
            points (WPO:Validate (car shape) (member (cdr (assoc 0 data)) '("LWPOLYLINE" "POLYLINE"))))
      (if (not (entmake (append
            (list '(0 . "LWPOLYLINE") '(100 . "AcDbEntity") (cons 8 (getvar "CLAYER"))
                  '(100 . "AcDbPolyline") (cons 90 (length points)) '(70 . 1)
                  '(43 . 0.0) (cons 38 (caddr shape)) (cons 210 (cadr shape)))
            (mapcar '(lambda (p) (cons 10 p)) points))))
        (WPO:Fail "No se pudo construir el contorno temporal."))
      (setq temp (entlast) previous temp)
      (command-s "_.WIPEOUT" "_Polyline" temp "_No")
      (if (and (not (eq previous (entlast))) (= (cdr (assoc 0 (entget (entlast)))) "WIPEOUT"))
        (setq mask (entlast))
        (WPO:Fail "CAD no pudo generar el WIPEOUT."))
      ;; Keep the mask above existing contents, with the outline above it.
      (command-s "_.DRAWORDER" mask "" "_Front")
      (command-s "_.DRAWORDER" en "" "_Above" mask "")
      mask) nil))
  (if (and temp (entget temp)) (entdel temp))
  (if (vl-catch-all-error-p result)
    (progn
      (if (and mask (entget mask)) (entdel mask))
      (list nil (cond (WPO:reason) ((vl-catch-all-error-message result)))))
    (list T result)))

(defun WPO:CurveSpan (en start end steps sweep / i a b p q corner midpoint factor)
  (setq i 0 p (vlax-curve-getPointAtParam en start) WPO:curve (cons p WPO:curve))
  (repeat steps
    (setq a (+ start (* (- end start) (/ (float i) steps)))
          b (+ start (* (- end start) (/ (float (1+ i)) steps)))
          q (vlax-curve-getPointAtParam en b))
    (if sweep
      (progn
        ;; Tangent intersection from the chord midpoint and arc midpoint. Works
        ;; for circles and affine ellipses, and avoids derivative discontinuities
        ;; at polyline vertices. The polygon encloses rather than cuts the arc.
        (setq midpoint (mapcar '(lambda (x y) (/ (+ x y) 2.0)) p q)
              factor (+ 1.0 (/ 1.0 (cos (/ sweep (* 2.0 steps)))))
              corner (mapcar '+ midpoint
                (mapcar '(lambda (x) (* x factor))
                  (mapcar '- (vlax-curve-getPointAtParam en (/ (+ a b) 2.0)) midpoint))))
        (setq WPO:curve (cons corner WPO:curve))))
    (setq WPO:curve (cons q WPO:curve) p q i (1+ i))))

(defun WPO:Curve (en / data kind start end bulges item vertex index bulge steps WPO:curve)
  (setq data (entget en) kind (cdr (assoc 0 data))
        start (vlax-curve-getStartParam en) end (vlax-curve-getEndParam en))
  (cond
    ((= kind "LINE") (WPO:CurveSpan en start end 1 nil))
    ((member kind '("CIRCLE" "ARC" "ELLIPSE"))
      (WPO:CurveSpan en start end (max 1 (WPO:Ceil (/ (- end start) (/ pi 90.0)))) (- end start)))
    (T
      (if (= kind "LWPOLYLINE")
        (foreach item data
          (cond ((= (car item) 10) (setq bulges (cons 0.0 bulges)))
                ((= (car item) 42) (setq bulges (cons (cdr item) (cdr bulges))))))
        (progn
          (if (/= 0 (logand (cdr (assoc 70 data)) 94))
            (WPO:Fail "No se admiten mallas ni polilineas 3D o ajustadas con Fit/Spline."))
          (setq vertex (entnext en))
          (while (= (cdr (assoc 0 (entget vertex))) "VERTEX")
            (setq bulges (cons (cond ((cdr (assoc 42 (entget vertex)))) (0.0)) bulges)
                  vertex (entnext vertex)))))
      (setq bulges (reverse bulges) index 0)
      (while (< index (fix end))
        (setq bulge (nth index bulges)
              steps (max 1 (WPO:Ceil (/ (abs (* 4.0 (atan bulge))) (/ pi 90.0)))))
        (WPO:CurveSpan en index (1+ index) steps (if (/= bulge 0.0) (* 4.0 (atan bulge))))
        (setq index (1+ index)))))
  WPO:curve)

(defun WPO:HalfHull (points / hull p)
  (foreach p points
    (while (and (cadr hull) (<= (WPO:Cross (cadr hull) (car hull) p) 0.0))
      (setq hull (cdr hull)))
    (setq hull (cons p hull)))
  hull)

(defun WPO:Hull (points / sorted lower upper)
  (setq sorted (vl-sort points '(lambda (a b) (if (= (car a) (car b)) (< (cadr a) (cadr b)) (< (car a) (car b)))))
        lower (WPO:HalfHull sorted) upper (WPO:HalfHull (reverse sorted)))
  (append (cdr lower) (cdr upper)))

(defun WPO:Geometry (ss / selected index en)
  (setq selected (ssadd) index 0)
  (repeat (sslength ss)
    (setq en (ssname ss index) index (1+ index))
    (if (member (cdr (assoc 0 (entget en))) '("LINE" "ARC" "CIRCLE" "ELLIPSE" "LWPOLYLINE" "POLYLINE"))
      (ssadd en selected)))
  selected)

(defun WPO:Envelope (ss / WPO:reason WPO:cloud WPO:count result index en data points normal elevation temp mask p originals)
  (setq WPO:count 0)
  (setq result (vl-catch-all-apply
    '(lambda ()
      (setq originals ss ss (WPO:Geometry ss))
      (if (= (sslength ss) 0)
        (WPO:Fail "No hay lineas, polilineas, arcos, circulos ni elipses en la seleccion."))
      (if (/= 0 (logand 4 (cdr (assoc 70 (tblsearch "LAYER" (getvar "CLAYER"))))))
        (WPO:Fail "La capa actual esta bloqueada."))
      (setq index 0)
      (repeat (sslength ss)
        (setq en (ssname ss index) data (entget en) index (1+ index))
        (if (/= (cdr (assoc 410 data)) (if (= (getvar "CVPORT") 1) (getvar "CTAB") "Model"))
          (WPO:Fail "Seleccione geometria del espacio actual."))
        (if (/= 0 (logand 4 (cdr (assoc 70 (tblsearch "LAYER" (cdr (assoc 8 data)))))))
          (WPO:Fail "La seleccion incluye capas bloqueadas. Desbloqueelas para ordenar la mascara debajo."))
        (foreach p (WPO:Curve en)
          (setq WPO:cloud (cons p WPO:cloud) WPO:count (1+ WPO:count)))
        (if (> WPO:count 50000) (WPO:Fail "Seleccione menos elementos (limite de 50000 puntos).")))
      (if (< (length WPO:cloud) 3) (WPO:Fail "Seleccione geometria suficiente para encerrar un area."))
      (setq normal (trans '(0.0 0.0 1.0) 1 0 T)
            points (mapcar '(lambda (p) (trans p 0 normal)) WPO:cloud)
            elevation (caddr (car points)))
      (if (not (vl-every '(lambda (p) (equal (caddr p) elevation 1e-6)) points))
        (WPO:Fail "La seleccion debe ser plana y paralela al XY del UCS actual."))
      (setq points (WPO:Validate (WPO:Hull points) nil))
      (if (> (length points) 2048) (WPO:Fail "Envolvente demasiado compleja."))
      (if (not (entmake (append
        (list '(0 . "LWPOLYLINE") '(100 . "AcDbEntity") (cons 8 (getvar "CLAYER"))
          '(100 . "AcDbPolyline") (cons 90 (length points)) '(70 . 1)
          (cons 38 elevation) (cons 210 normal))
        (mapcar '(lambda (p) (cons 10 (list (car p) (cadr p)))) points))))
        (WPO:Fail "No se pudo crear la envolvente temporal."))
      (setq temp (entlast) result (WPO:Create temp))
      (if (not (car result)) (WPO:Fail (cadr result)))
      (setq mask (cadr result))
      ;; Ignored annotations do not enlarge the envelope, but remain visible.
      (command-s "_.DRAWORDER" originals "" "_Above" mask "")
      mask) nil))
  (if (and temp (entget temp)) (entdel temp))
  (if (vl-catch-all-error-p result)
    (progn
      (if (and mask (entget mask)) (entdel mask))
      (list nil (cond (WPO:reason) ((vl-catch-all-error-message result)))))
    (list T result)))

(defun c:WPO (/ *error* echo opened ss result)
  (defun *error* (msg)
    (if opened (command-s "_.UNDO" "_End"))
    (if echo (setvar "CMDECHO" echo))
    (if (and msg (not (wcmatch (strcase msg) "*CANCEL*,*QUIT*,*EXIT*,*BREAK*")))
      (princ (strcat "\nWPO: " msg)))
    (princ))
  (setq ss (ssget "_I"))
  (if (null ss)
    (progn (sssetfirst nil nil)
           (princ "\nSeleccione la geometria que debe abarcar el WIPEOUT (se ignoran textos, cotas y bloques): ")
           (setq ss (ssget))))
  (if ss
    (progn
      (setq echo (getvar "CMDECHO")) (setvar "CMDECHO" 0)
      (command-s "_.UNDO" "_Begin") (setq opened T)
      (setq result (WPO:Envelope ss))
      (command-s "_.UNDO" "_End") (setq opened nil)
      (setvar "CMDECHO" echo)
      (if (car result)
        (princ "\nWPO: envolvente exterior creada bajo la seleccion. Originales conservados; UNDO permite deshacer.")
        (princ (strcat "\nWPO: " (cadr result))))))
  (princ))

(princ)
