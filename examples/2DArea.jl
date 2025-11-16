using MinGal

Algebra(2)

point = (x, y) -> return x*e1 + y*e2

# Calculates the area of a polygon by an ordered set of points in geometric-algebra terms by triangulating 
# the polygon with a fixed reference point and summing the areas of the triangles.
function area(points)
    area = 0
    p1 = points[1]
    for i in 2:length(points)-1
        area += (points[i]-p1)^(points[i+1]-p1)
    end

    return norm(area)/2
end

pts = [point(1,2), point(7.4, 2.3), point(4.17, -4.3), point(1.38, -1.12)]
println(area(pts))