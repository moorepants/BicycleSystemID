load fits

for j = 1:6
    figure()
    hold on
    for i = 1:length(speeds)
        if fits(i) > 50
            if strcmp(riders{i}, 'Jason')
                color = 'b.';
            elseif strcmp(riders{i}, 'Charlie')
                color = 'g.';
            elseif strcmp(riders{i}, 'Luke')
                color = 'r.';
            else
                color = 'k.';
            end
            plot(speeds(i), parameters(i, j), color)
        end
    end
    hold off
end
