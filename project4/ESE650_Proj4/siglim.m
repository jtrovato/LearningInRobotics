function y = siglim(x)
y = 1./(1+exp(-x));
if y == 1
    y = 0.99999;
elseif y == 0
    y = 0.000001;
end

end

