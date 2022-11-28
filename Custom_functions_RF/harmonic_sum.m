function f = harmonic_sum(amp, len)

    if isrow(amp)

        amp = amp';

    end

    x = linspace(0,pi, len);
    f = sum(cell2mat(arrayfun(@(n,a) a*sin(n*x), (1:length(amp))', amp, 'UniformOutput', false)))';
    %x = linspace(-1,1, len);
    %f = sum(cell2mat(arrayfun(@(n,a) a*chebyshevT(n, x), (1:length(amp))', amp, 'UniformOutput', false)))';

end