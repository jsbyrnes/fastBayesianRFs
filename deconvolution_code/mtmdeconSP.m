function [result, fft_result] = mtmdeconSP(source_timeseries,converted_phases_timeseries, n, number_of_tapers, nmtw, E, V)

segments = length(0:number_of_tapers/4:n-number_of_tapers);
tapers   = zeros(segments,n,nmtw);

[m, ~] = size(source_timeseries);

for i = 1:length(V)
    
    for j = 1:segments
        
        tmp = 1:number_of_tapers;
        
        tapers(j, tmp + number_of_tapers*(j-1)/4,i) = E(:,i);
    
    end
    
end

first_tapered = zeros(n,nmtw);
second_tapered = zeros(n,nmtw);

for i = 1:m
    
    [converted_gridded, ~] = meshgrid(converted_phases_timeseries(i,:),1:segments);
    
    [source_gridded, ~] = meshgrid(source_timeseries(i,:),1:segments);
    
    for j = 1:length(V)
        
        first_tapered(:,j) = sum(fft(tapers(:,:,j).*source_gridded*V(j),[],2),1);
        
        second_tapered(:,j) = sum(fft(tapers(:,:,j).*converted_gridded*V(j),[],2),1);
                
    end

    num=sum(second_tapered.*conj(first_tapered),2);
    denom=sum(first_tapered.*conj(first_tapered),2);
    
    fft_result = num./(denom);
        
    result(:, i) = real(ifft(fft_result));
    
end

result = sum(result, 2);