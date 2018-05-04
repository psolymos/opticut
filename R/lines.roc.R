lines.roc
function(x, ...)
{
    lines(x$fpr, x$tpr, ...)
    invisible(x)
}
