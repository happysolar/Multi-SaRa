llce <- function(y,h) {
    yy = c(rep(0,h-1),y,rep(0,h))            # add some zeros on the head and tail
    n = length(y)
    z = rep(0,n)
    ans<-.C("diagnosticValue", yy=as.double(yy), h=as.integer(h), 
            n=as.integer(n), z=as.double(z))

    return (ans$z)                                                                 #
}
