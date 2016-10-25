check_strata <-
function(x, mat)
{
    Str <- as.integer(strata(x))
    BBstr <- apply(data.matrix(mat), 2, function(z) Str[z])
    nstr <- apply(BBstr, 2, function(z) length(unique(z)))
    structure(nstr == length(unique(Str)),
        nx = length(unique(Str)),
        nmat = nstr)
}
