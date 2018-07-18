function [ color ] = coloring( on,req)
for i = 1:length(on)
    if on(i) && req(i)
        color(i) = 10;
    elseif (~req(i)) && on(i)
        color(i) = 0;
    else
        color(i) = 5;
    end
end
end

