abstract type Galaxy end

struct Spiral <: Galaxy
	name::Symbol
	diskHeight::Number
	radius::Number
	units::Quantity
end