function out = gaussianComposition(height, mu, sig, x)

out = zeros(length(x),1)';

for n = 1:length(height)
    
    out = out + height(n)*exp(-((x - mu(n))/sig(n)).^2);
    
end
    
    