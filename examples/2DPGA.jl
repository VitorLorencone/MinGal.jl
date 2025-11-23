using MinGal, CairoMakie

Algebra(2,0,1)

# -------------------------- Constructions --------------------------

mutable struct Point
    mv::GAType
end

mutable struct Segment
    mv::GAType
    start_point::Point
    end_point::Point
end

# P = (x, y)
point = (x::Number, y::Number) -> return Point(-x*e0e2 + y*e0e1 + e1e2)

# S: Line from point `a` to point `b`
segment = (a::Point, b::Point) -> return Segment(a.mv & b.mv, a, b)

# -------------------------- Functions --------------------------

# Translate Point 'a' by (tx,ty)
function translate(a::Point, tx, ty)
    T = 1 + (tx*e1 + ty*e2)*e0/2
    return Point(T * a.mv * T^-1)
end

# Translate Segment 'a' by (tx,ty)
function translate(a::Segment, tx, ty)
    T = 1 + (tx*e1 + ty*e2)*e0/2
    return Segment(T*a.mv*T^-1, Point(T*a.start_point.mv*T^-1), Point(T*a.end_point.mv*T^-1))
end

# Rotate Point 'a' by an angle 'n' rad from the origin (0,0)
function rotate(a::Point, n)
    R = exp_ga(-e1e2*n/2)
    return Point(R * a.mv * R^-1)
end

# Rotate Segment 'a' by an angle 'n' rad from the origin (0,0)
function rotate(a::Segment, n)
    R = exp_ga(-e1e2*n/2)
    return Segment(R*a.mv*R^-1, Point(R*a.start_point.mv*R^-1), Point(R*a.end_point.mv*R^-1))
end

# -------------------------- Scene --------------------------

# Starts the scene
fig = Figure()
ax = Axis(fig[1,1])

# Construct a house with segments
objects = Segment[]
push!(objects, segment(point(1, 2), point(1, 8)))
push!(objects, segment(point(1, 8), point(-3, 8)))
push!(objects, segment(point(-3, 8), point(-3, 2)))
push!(objects, segment(point(-3, 2), point(1, 2)))
push!(objects, segment(point(-3, 2), point(-6, 5)))
push!(objects, segment(point(-6, 5), point(-3, 8)))
push!(objects, segment(point(1, 6), point(-2, 6)))
push!(objects, segment(point(1, 4), point(-2, 4)))
push!(objects, segment(point(-2, 6), point(-2, 4)))

# Translate everything (-1, -2) -> Origin
objects = map(a -> translate(a, -1, -2), objects)

# Rotate everything -pi/2 rad -> Right direction
objects = map(a -> rotate(a, -pi/2), objects)

# Render everything with CairoMakie
for obj in objects
    lx = [-obj.start_point.mv[e0e2], -obj.end_point.mv[e0e2]]
    ly = [obj.start_point.mv[e0e1], obj.end_point.mv[e0e1]]
    lines!(ax, lx, ly, color = :red, linewidth = 2)
end

# The house after the transformations
save("examples/2DPGA.png", fig)