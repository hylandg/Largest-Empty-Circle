function lec_qhull, x, y, plot_tri=plot_tri, plot_vor=plot_vor

; Purpose:
;   Solves the largest empty circle (LEC) problem for a set of 2-d points.
;   This consists of finding the largest circle that contains none of the
;   2-d points and is also centred inside the convex hull of the points.
;   This solution uses the IDL routine QHULL to determine the Delaunay
;   triagulation and the Veronoi diagram of the set of points.
;
; Inputs:
;   x:  An array containing the x locations of the points
;   y:  An array containing the y locations of the points
;
; Keyword Inputs:
;   plot_tri:  Plot the Delaunay triagulation
;   plot_vor:  Plot the Veronoi diagram
;
; Outputs:
;   LEC_QHULL returns a structure containing the tags {radius,x,y} that
;   specifies the radius and centre location of the largest empty circle.
;
; Modification History:
;   Written by: Glenn Hyland, Australian Antarctic Division & Antarctic
;   Climate & Ecosystems CRC, July 2015

qhull, x, y, tr, /delaunay, vdiagram=vdiagram, vvertices=vvert, vnorm=vnorm, $
		bounds=b

; QHULL returns the convex hull points unsorted, but the "inside_poly" function
; needs a sorted polygon. So sort into anti-clockwise order using the angle
; of each point from the centre of the set.

b = b[sort(atan(y[b]-mean(y[b]),x[b]-mean(x[b])))]

; Determine end points for each unbounded Veronoi ridge. This will allow us
; to plot the ridge, and also determine if it intersects with the convex hull.
;
; Note: the solution used is mathematically correct, but due to numerical
; precision issues and the occasional error from QHULL (in vnorm) it will
; give an incorrect end point. A forced check is included later on to try
; and catch these problems.

xmin = min(x[b], max=xmax)
ymin = min(y[b], max=ymax)
vunvert = dblarr(3,n_elements(vnorm[0,*]))

unbnd = where(vdiagram[2,*] lt 0, nunbnd)
for ii = 0, nunbnd-1 do begin

	vdiag = vdiagram[*,unbnd[ii]]

	; For each unbounded ridge QHULL returns the line equation in vnorm,
	; but we need to work out the direction of the line from the single
	; known ridge vertex. Start by locating all points that refer to the
	; known vertex.

	iv = where(vdiagram[2,*] eq vdiag[3] or vdiagram[3,*] eq vdiag[3], nv)
	jv = [vdiagram[0,iv],vdiagram[1,iv]]
	idx = where(jv ne vdiag[0] and jv ne vdiag[1])
	kv = jv[idx[0]]

	; So vdiag[0:1] are the 2 points that the ridge bisects and kv is
	; the next closest point. Now move a set distance from the known 
	; vertex along the ridge and see how this affects the distance back
	; to vdiag[0:1] and kv. If we're heading in the right direction, kv
	; should be the furthest point (by definition of the Voronoi diagram).

	j = -vdiag[2] - 1
	xystart = vvert[*,vdiag[3]]
	if (abs(vnorm[1,j]) gt abs(vnorm[0,j])) then begin
		dx = double(xmax - xmin)
		dy = -vnorm[0,j]*dx/vnorm[1,j]
	endif else begin
		dy = double(ymax - ymin)
		dx = -vnorm[1,j]*dy/vnorm[0,j]
	endelse
	dxkv = x[vdiag[0]] - x[kv]
	dykv = y[vdiag[0]] - y[kv]
	dd = dx*dxkv + dy*dykv

	vunvert[0:1,j] = (dd lt 0.)? xystart-[dx,dy] : xystart+[dx,dy]

	; Save the distance to the next closest point for the forced check on
	; unbound ridge intersections later on.

	vunvert[2,j] = dxkv^2 + dykv^2
endfor

; Optionally plot the Delaunay triagulation (contained in variable "tr")

if (keyword_set(plot_tri)) then begin
	plots, psym=2, x[b], y[b]
	for ii = 0, n_elements(tr[0,*])-1 do $
		plots, x[[tr[*,ii],tr[0,ii]]], y[[tr[*,ii],tr[0,ii]]]
endif

; Optionally plot the Veronoi diagram

if (keyword_set(plot_vor)) then begin
	for ii = 0, n_elements(vdiagram[2,*])-1 do begin
		vdiag = vdiagram[*,ii]
		if (vdiag[2] ge 0) then $
			plots, vvert[*,vdiag[2:3]], linestyle=1, noclip=0 $
		else begin
			j = -vdiag[2] - 1
			plots, [[vvert[*,vdiag[3]]],[vunvert[0:1,j]]], line=2, $
				noclip=0
		endelse
	endfor
endif

; The LEC is always centred at either a Voronoi vertex (within the convex
; hull), or on an intersection between a Voronoi ridge and the convex hull.

s = ''
dmax = 0.

; Firstly, check each vertex. Do this via vdiagram because it already has
; links to vvert (going the other way involves using "where" inefficiently).
; The largest radius is the distance from the vertex to the closest point
; (which is any of the bisecting points associated with the vertex).

k = uniq(vdiagram[3,*], sort(vdiagram[3,*]))
for i = 0, n_elements(k)-1 do begin
	iv = vdiagram[3,k[i]]
	if (inside_poly_crossings(vvert[0,iv], vvert[1,iv], x[b], y[b], /convex)) then begin
		dx = vvert[0,iv] - x[vdiagram[0,k[i]]]
		dy = vvert[1,iv] - y[vdiagram[0,k[i]]]
		dmin = dx^2 + dy^2
		if (dmin gt dmax) then begin
			dmax = dmin
			pmax = reform(vvert[*,iv])
		endif
	endif
endfor
if (keyword_set(plot_tri) or keyword_set(plot_vor)) then begin
	tvcircle, sqrt(dmax), pmax[0], pmax[1], /data, line=2
	plots, pmax, psym=5
endif

; Finally, check the intersections of Voronoi ridges and the convex hull.
; The largest radius is the distance from the intersection to either of
; the bisecting points associated with the ridge.

; Note: We do a special check on unbound ridge intersections because sometimes
; the results from QHULL (in vnorm) can be wrong.

for i = 0, n_elements(vdiagram[2,*])-1 do begin
	vdiag = reform(vdiagram[*,i])
	if (vdiag[2] ge 0) then begin
		unbound = 0
		pv3 = vvert[*,vdiag[2]]
		pv4 = vvert[*,vdiag[3]]
	endif else begin
		unbound = -vdiag[2] - 1
		pv3 = vvert[*,vdiag[3]]
		pv4 = vunvert[0:1,unbound]
;		plots, pv3, psym=6, symsize=2
	endelse

	pb1 = [x[b[-1]],y[b[-1]]]
	for j = 0, n_elements(b)-1 do begin
		pb2 = [x[b[j]],y[b[j]]]
		p = segment_intersect(pb1, pb2, pv3, pv4, unbound=unbound)
		if (n_elements(p) gt 1) then begin

			dmin = (p[0] - x[vdiag[0]])^2 + (p[1] - y[vdiag[0]])^2

			; This is the special check on unbound ridge
			; intersections - if it is valid it should be closer
			; then the next closest point (previously saved in
			; the vunvert array).

			if ((unbound gt 0) and $
				(dmin gt vunvert[2,unbound])) then dmin = 0.

			if (dmin gt dmax) then begin
				dmax = dmin
				pmax = p
			endif
;			if (unbound) then plots, p, psym=4, symsize=2
		endif
		pb1 = pb2
	endfor
endfor

return, {radius:sqrt(dmax), x:pmax[0], y:pmax[1]}
end
