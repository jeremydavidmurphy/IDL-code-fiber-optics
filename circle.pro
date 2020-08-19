FUNCTION CIRCLE, xcenter, ycenter, radius
points = (2 * !PI / 499.0) * FINDGEN(500)
x = xcenter + radius * COS(points )
y = ycenter + radius * SIN(points )
RETURN, TRANSPOSE([[x],[y]])
END
