function segment_intersect, p1, p2, p3, p4, unbound=unbound, p3equ=p3equ

; Determine the intersection of the line segments (p1,p2) and (p3,p4). The
; intersection point must be contained within each line segment unless the
; "unbound" keyword is set, in which case the intersection point can extend
; beyond p4 only.
;
; Keyword "p3equ" allows the line equation for (p3,p4) to be specified as
; ux + vy + w = 0 where p3=[u,v,w] (keyword "unbound" and p4 are ignored).
;
; Written by: Glenn Hyland, Australian Antarctic Division & Antarctic
; Climate & Ecosystems CRC, July 2015

if (not keyword_set(p3equ)) then begin

	d12 = p2 - p1
	d34 = p4 - p3

	d13 = p1 - p3

	d1 = d12[0]*d34[1] - d12[1]*d34[0]
	if (d1 eq 0.) then return, -1

	s = (d12[0]*d13[1] - d12[1]*d13[0]) / d1
	t = (d34[0]*d13[1] - d34[1]*d13[0]) / d1

	p = (s ge 0. and (s le 1. or keyword_set(unbound)) and $
		t ge 0. and t le 1.)? p3 + s*d34 : -1

	return, p

endif else begin

	dx = p2[0] - p1[0]
	dy = p2[1] - p1[1]
	if (dx eq 0 and dy eq 0) then return, -1

	dxy = p1[0]*p2[1] - p2[0]*p1[1]

	y = -(p3[0]*dxy + p3[2]*dy) / (p3[0]*dx + p3[1]*dy)

	if (dy ne 0.) then begin
		s = (y - p1[1]) / dy
		if (s lt 0. or s gt 1.) then return, -1
		x = (y*dx + dxy) / dy
	endif else begin
		if (p3[0] eq 0.) then return, -1
		x = -(p3[1]*y + p3[2]) / p3[0]
		s = (x - p1[0]) / dx
		if (s lt 0. or s gt 1.) then return, -1
	endelse

	return, [x,y]

endelse

end
