function colors = getColors(z)
    colors = 0.5*ones(length(z), 3);
    colors(z > 5,:) = [0 0 0];
    colors(z<-5,:) = [1 1 1];
end

