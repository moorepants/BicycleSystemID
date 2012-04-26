outputs = {'delta'};
inputs = {'fB'};

dir = '../../scripts/statespaceid/exports';

data = build_id_data('00700.mat', outputs, inputs, dir, true);

[newData, percent] = remove_process_noise(data);
percent

plot(data, newData)
