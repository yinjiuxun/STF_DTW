function normalized_series = series_normalization(series_M,normalization_methods)
% normalize the series (in matrix form) amplitude
% normalized_series = series_normalization(T,series_M)
% series_M is a Matrix of series with size (N_series,N_pts)
% normalization_methods can be: "max_amp", "integral" (area)

switch normalization_methods
    case "max_amp"
        normalized_series = series_M./max(series_M,[],2);
        
    case "integral"
        normalized_series = series_M./trapz(series_M,2);
        
end
        


