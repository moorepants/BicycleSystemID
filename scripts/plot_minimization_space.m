outputs = {'phiDot'};
[z, speed, rider] = build_id_data('00700.mat', outputs, 'pavilion/mat');

% this was the best solution I found so far
guess = [9.7676, -1.5008, 3.0039, 0.7115, 0.1670, 43.3128];
m = bicycle_grey(['Rigid' rider], speed, outputs, guess(1:5), guess(6));

figure()
compare(z, m)

par1 = linspace(-1.8, -1, 30);
par2 = linspace(0.05, 0.3, 30);
s = zeros(length(par1), length(par2));

for i = 1:length(par1)
    for j = 1:length(par2)
        display(sprintf('Simulating with %1.5f and %1.5f.', par1(i), par2(j)))
        guess = [9.7676, par1(i), 3.0039, 0.7115, par2(j), 43.3128];
        m.ParameterVector = guess;
        y = sim(m, z);
        s(i, j) = sum((y.OutputData - z.OutputData).^2);
    end
end

figure()
surf(par1, par2, s)
