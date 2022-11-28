function out = returngabor(A,B,t,f,w,X)
    out = exp(-(X-t).^2/(2*(w)^2)).*(A*cos(2*pi*(X-t)*f) + B*sin(2*pi*(X-t)*f));%windowed over one cycle
end