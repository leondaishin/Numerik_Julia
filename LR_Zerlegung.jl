using LinearAlgebra

# LR-Zerlegung: A = L * R
function LRdecomposition(A)
    n = size(A, 1)

    L = Matrix{Float64}(I, n, n)
    R = Float64.(A)

    for j in 1:n-1
        abs(R[j, j]) > eps(Float64) ||
            error("Null-Pivot in Spalte $j. Die LR-Zerlegung ohne Pivotisierung ist nicht möglich.")

        for i in j+1:n
            faktor = R[i, j] / R[j, j]
            L[i, j] = faktor

            for k in j:n
                # Gauß-Elimination
                R[i, k] = R[i, k] - faktor * R[j, k]
            end
        end
    end

    return L, R
end

# Vorwärtssubstitution: Löse L*y = b
function vorSub(L, b)
    n = length(b)
    y = zeros(Float64, n)

    for i in 1:n
        s = 0.0
        for j in 1:i-1
            s += L[i, j] * y[j]
        end
        y[i] = (b[i] - s) / L[i, i]
    end

    return y
end

# Rückwärtssubstitution: Löse R*x = c
function rueckSub(R, c)
    n = length(c)
    x = zeros(Float64, n)

    for i in n:-1:1
        s = 0.0
        for j in i+1:n
            s += R[i, j] * x[j]
        end
        x[i] = (c[i] - s) / R[i, i]
    end

    return x
end

# Eingabe der Matrix A und des Vektors b
function eingabe()
    println("Wie groß soll n sein für A ∈ K^(n×n)?")
    n = parse(Int, readline())

    A = zeros(Float64, n, n)

    println("\nGib jetzt die Einträge der Matrix A ein:")
    for i in 1:n
        for j in 1:n
            print("a", i, j, " = ")
            A[i, j] = parse(Float64, readline())
        end
    end

    b = zeros(Float64, n)

    println("\nGib jetzt die Einträge des Vektors b ein:")
    for i in 1:n
        print("b", i, " = ")
        b[i] = parse(Float64, readline())
    end

    return A, b
end

# Hauptprogramm
A, b = eingabe()

L, R = LRdecomposition(A)

y = vorSub(L, b)
x = rueckSub(R, y)

println("\nMatrix A:")
display(A)

println("\nMatrix L:")
display(L)

println("\nMatrix R:")
display(R)

println("\nLösung y von L*y = b:")
display(y)

println("\nLösung x von R*x = y:")
display(x)