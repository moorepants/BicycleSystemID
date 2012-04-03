load fits

outputs = {'phiDot'};

for i = 1:length(speeds)
    if fits(i) > 30
        [z, speed, rider] = build_id_data(matFiles{i}, outputs, 'pavilion/mat');
        guess = parameters(i, :);
        m = bicycle_grey(['Rigid' rider], speed, outputs, guess(1:5), guess(6));
        compare(z, m);
        pause
        close all
    end
end
