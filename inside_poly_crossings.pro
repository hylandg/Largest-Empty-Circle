function inside_poly_crossings, xp, yp, vx, vy, convex=convex

; Purpose:
;   Test whether a given point (xp,yp) is inside a polygon (vx,vy) using
;   the "Crossings" algorithm - the point is inside the polygon if, for any
;   ray from this point, there is an odd number of crossings of the ray with
;   the polygon's edges
;
; Inputs:
;   It is assumed that the polygon (vx,vy) has n vertices with n sides and the
;   edges connect the vertices in the order:
;     [(vx1,vy1), (vx2,vy2), ..., (vxn,vyn), (vx1,vy1)].
;   i.e. the last vertex is connected to the first vertex.
;
;   xp:  The x coordinate location of the point
;   yp:  The y coordinate location of the point
;   vx:  A vector of x coordinate locations for the polygon vertices
;   vy:  A vector of y coordinate locations for the polygon vertices
;
; Keyword Inputs:
;   convex:  Indicates the polygon (vx,vy) is convex. This allows a faster
;            algorithm to be used due to the polygon's geometric properties.
;
; Outputs:
;   INSIDE_POLY_CROSSINGS returns TRUE (=1) if the given point is inside the
;   polygon. Otherwise the functions returns FALSE (=0).
;
; Method:
;   Shoot a test ray along +X axis. The strategy is to compare vertex Y values
;   to the testing point's Y and quickly discard edges which are entirely to
;   one side of the test ray. This version is usually somewhat faster than
;   the original published in Graphics Gems IV; by turning the division for
;   testing the X axis crossing into a tricky multiplication test this part
;   of the test became faster on machines where division is expensive. Your
;   mileage may (in fact, will) vary, depending on the machine and the test
;   data, but in general I believe this code is both shorter and faster.
;   This test was inspired by unpublished Graphics Gems submitted by Joseph
;   Samosky and Mark Haigh-Hutchinson.
;
; Modification History:
;   Code taken from ptinpoly.c - point in polygon inside/outside code by
;   Haines, Eric, "Point in Polygon Strategies," Graphics Gems IV, ed. Paul
;   Heckbert, Academic Press, p. 24-46, 1994.
;
;   Converted to IDL by: Glenn Hyland, Australian Antarctic Division &
;   Antarctic Climate & Ecosystems CRC, July 2015

nv = n_elements(vx)

; get initial test bit for above/below X axis
yflag0 = ( vy[nv-1] ge yp )

if (keyword_set(convex)) then line_flag = 0
inside_flag = 0

for j = 0, nv-1 do begin

	yflag1 = ( vy[j] ge yp )

	; Check if endpoints straddle (are on opposite sides) of the X axis
	; (i.e. the Y's differ); if so, the +X ray could intersect this edge.

	if ( yflag0 ne yflag1 ) then begin

		; Check intersection of polygon segment with the +X ray.
		; Note if >= point's X; if so, the ray hits it.
		; A division operation is avoided for the flag test by
		; checking the sign of the first vertex wrt the test point
		; from an idea inspired by Joseph Samosky's and Mark Haigh-
		; Hutchinson's different polygon inclusion tests.

		flag = (vy[j]-yp)*(vx[j-1]-vx[j]) ge (vx[j]-xp)*(vy[j-1]-vy[j])
		if ( flag eq yflag1 ) then inside_flag = 1 - inside_flag

		; For convex polygons, the crossings test can quit as soon as
		; two Y sign difference edges are found, since this is the
		; maximum that a convex polygon can have.

		if (keyword_set(convex)) then begin
	    		if (line_flag) then break
	    		line_flag = 1
		endif
	endif

	yflag0 = yflag1
endfor

return, inside_flag
end
