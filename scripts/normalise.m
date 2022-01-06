function Y = normalise(X,lower,upper)

thismin = min(X(:));
thismax = max(X(:));

Y = (X - thismin) / ( thismax - thismin );
Y = Y*(upper-lower) + lower;

end