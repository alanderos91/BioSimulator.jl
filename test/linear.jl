# From Cao, Li, & Petzold 2004
function linear(M, x0)
    m = Network("Linear Chain System")

    m <= Species(:S1, x0, istracked=true)

    for i = 2:(M+1)
        m <= Species(symbol("S$(i)"), 0, istracked=false)
    end

    for i = 1:M
        m <= Reaction(symbol("R$(i)"), symbol("k$(i)"),
            r=(symbol("S$(i)") => 1),
            p=(symbol("S$(i+1)") => 1),
            istracked=false
        )
        m <= parameter(symbol("k$(i)"), 1.0)
    end

    return m
end
