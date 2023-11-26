function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end