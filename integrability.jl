using Pkg
Pkg.add("SymPy")

using SymPy

# Define symbolic variables
@syms x y

# Define polynomials
p1 = x^2 + 2 * x + 1
p2 = x - 1

# Define a rational polynomial
rational_poly = p1 / p2

# Display the rational polynomial
println("Rational Polynomial: ", rational_poly)

# Perform symbolic operations
expanded_poly = expand(rational_poly)
simplified_poly = simplify(rational_poly)

# Display the results
println("Expanded Polynomial: ", expanded_poly)
println("Simplified Polynomial: ", simplified_poly)

# Define Markov's equation
@syms a b c
markovs_equation = a^2 + b^2 + c^2 - 3 * a * b * c

# Display Markov's equation
println("Markov's Equation: ", markovs_equation)

# Define a function to check if a polynomial morphism keeps Markov's equation invariant
function check_invariance(f, g, h)
    substituted_eq = markovs_equation(a => f, b => g, c => h)
    return simplify(substituted_eq) == markovs_equation
end

# Define low degree polynomial morphisms
@syms u v w
morphisms = [
    (u, v, w),
    (u + 1, v, w),
    (u, v + 1, w),
    (u, v, w + 1),
    (u^2, v, w),
    (u, v^2, w),
    (u, v, w^2)
]

# Check each morphism for invariance
for (f, g, h) in morphisms
    is_invariant = check_invariance(f, g, h)
    println("Morphism (f, g, h) = ($f, $g, $h) keeps Markov's equation invariant: ", is_invariant)
end
