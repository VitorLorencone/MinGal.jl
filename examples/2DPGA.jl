using MinGal, CairoMakie

Algebra(2,0,1)

# -------------------------- Constructions --------------------------

mutable struct Point
    mv::GAType
end

mutable struct Line
    mv::GAType
    start_point::Point
    end_point::Point
end

# P = (x, y)
point = (x, y) -> return Point(-x*e0e2 + y*e0e1 + e1e2)

# r: ax + by + c = 0
line = (a, b, c) -> Line(a*e1 + b*e2 + c*e0, nothing, nothing)

# -------------------------- Functions --------------------------

function join_line(a::Point, b::Point)
    res = a.mv & b.mv
    return Line(res, a, b)
end

# Translate 'a' by (tx,ty)
function translate(a, tx, ty)
    t = tx*e1 + ty*e2
    T = 1 + t*e0/2
    return T * a * invert(T)
end

# Rotate 'a' by an angle 'n' rad from the origin (0,0)
function rotate(a, n)
    rotor = exp_ga(-e1e2*n/2)
    return rotor * a * invert(rotor)
end

# -------------------------- Scene --------------------------

fig = Figure()
ax = Axis(fig[1,1])

# Construct a house with segments
objects = Line[]
push!(objects, join_line(point(1, 2), point(1, 8)))
push!(objects, join_line(point(1, 8), point(-3, 8)))
push!(objects, join_line(point(-3, 8), point(-3, 2)))
push!(objects, join_line(point(-3, 2), point(1, 2)))
push!(objects, join_line(point(-3, 2), point(-6, 5)))
push!(objects, join_line(point(-6, 5), point(-3, 8)))
push!(objects, join_line(point(1, 6), point(-2, 6)))
push!(objects, join_line(point(1, 4), point(-2, 4)))
push!(objects, join_line(point(-2, 6), point(-2, 4)))

# Translate everything (-1, -2) -> Origin
for obj in objects
    obj.mv = translate(obj.mv, -1, -2)
    obj.start_point = Point(translate(obj.start_point.mv, -1, -2))
    obj.end_point = Point(translate(obj.end_point.mv, -1, -2))
end

# Rotate everything -pi/2 rad -> Right direction
for obj in objects
    obj.mv = rotate(obj.mv, -pi/2)
    obj.start_point = Point(rotate(obj.start_point.mv, -pi/2))
    obj.end_point = Point(rotate(obj.end_point.mv, -pi/2))
end

# Render everything with CairoMakie
for obj in objects
    lx = [-obj.start_point.mv[e0e2], -obj.end_point.mv[e0e2]]
    ly = [obj.start_point.mv[e0e1], obj.end_point.mv[e0e1]]
    lines!(ax, lx, ly, color = :red, linewidth = 2)
end

# The house after the transformations
save("examples/2DPGA.png", fig)