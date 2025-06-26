function runCircuitAnalysis(inputFile)
    % Symbolic variables for Laplace transform
    s = sym('s');
    t = sym('t');


    % Read and parse input file
    try
        fid = fopen(inputFile, 'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        lines = lines{1};
        fclose(fid);
    catch
        error('Error reading input file');
    end

    % Parse circuit parameters
    numNodes = str2double(lines{1});
    if isnan(numNodes) || numNodes < 1
        error('Invalid number of nodes');
    end

    % Parse time parameters
    timeParamsLine = find(contains(lines, 'TIME_PARAMS:'), 1);
    if isempty(timeParamsLine)
        error('Time parameters not found in input file');
    end
    timeParams = sscanf(lines{timeParamsLine + 1}, '%f %f %f');
    tStart = timeParams(1);
    tEnd = timeParams(2);
    timeStep = timeParams(3);

    % Parse components
    componentsLine = find(contains(lines, 'COMPONENTS:'), 1);
    if isempty(componentsLine)
        error('Components section not found in input file');
    end

    components = {};
    i = componentsLine + 1;
    while i <= length(lines) && ~contains(lines{i}, 'VOLTAGE_SOURCE:')
        tokens = strsplit(strtrim(lines{i}));
        if length(tokens) < 4
            error('Invalid component format at line %d.', i);
        end
        components{end + 1} = {str2double(tokens{1}), str2double(tokens{2}), tokens{3}, str2double(tokens{4})};
        i = i + 1;
    end

    % Parse multiple voltage sources
    voltageSourceLine = find(contains(lines, 'VOLTAGE_SOURCE:'), 1);
    voltageSources = [];
    i = voltageSourceLine + 1;
    while i <= length(lines) && ~contains(lines{i}, 'MEASURE_VOLTAGE:')
        voltageData = sscanf(lines{i}, '%d %d %f');
        voltageSources = [voltageSources; voltageData'];
        i = i + 1;
    end

    % Initialize symbolic matrices
    G = sym(zeros(numNodes-1));  % Conductance matrix
    I = sym(zeros(numNodes-1, 1));  % Current vector

    % Process components and build system matrices
    for idx = 1:length(components)
        comp = components{idx};
        node1 = comp{1};
        node2 = comp{2};
        type = comp{3};
        value = comp{4};

        switch type
            case 'R'  % Resistor
                admittance = 1 / value;
            case 'L'  % Inductor
                admittance = 1 / (s * value);
            case 'C'  % Capacitor
                admittance = s * value;
            otherwise
                error('Unknown component type: %s', type);
        end
        addComponentToMatrix(admittance, node1, node2);
    end

    % Handle voltage sources in MNA
    numG = size(G, 1);
    G_ext = [G, zeros(numG, size(voltageSources, 1)); zeros(size(voltageSources, 1), numG + size(voltageSources, 1))];
    I_ext = [I; sym(zeros(size(voltageSources, 1), 1))];

    for idx = 1:size(voltageSources, 1)
        node1 = voltageSources(idx, 1);
        node2 = voltageSources(idx, 2);
        value = voltageSources(idx, 3) / s;  % Step input in Laplace

        if node1 > 0
            G_ext(node1, numG + idx) = 1;
            G_ext(numG + idx, node1) = 1;
        end
        if node2 > 0
            G_ext(node2, numG + idx) = -1;
            G_ext(numG + idx, node2) = -1;
        end
        I_ext(numG + idx) = value;
    end

    % Solve system
    nodeVoltagesLaplace = G_ext \ I_ext;

    % Generate time vector
    t_vector = tStart:timeStep:tEnd;

    % Convert node voltages to time domain
    nodeVoltagesTime = zeros(numNodes-1, length(t_vector));
    for i = 1:numNodes-1
        voltageSymbolic = ilaplace(nodeVoltagesLaplace(i), s, t);
        for j = 1:length(t_vector)
            nodeVoltagesTime(i, j) = double(subs(voltageSymbolic, t, t_vector(j)));
        end
    end

    % Calculate component currents
    componentCurrentsTime = zeros(length(components), length(t_vector));
    for idx = 1:length(components)
        comp = components{idx};
        node1 = comp{1};
        node2 = comp{2};
        type = comp{3};
        value = comp{4};

        if node1 == 0  % If node1 is ground
    voltageDiff = -nodeVoltagesTime(node2, :);
elseif node2 == 0  % If node2 is ground
    voltageDiff = nodeVoltagesTime(node1, :);
else  % If both nodes are non-ground
    voltageDiff = nodeVoltagesTime(node1, :) - nodeVoltagesTime(node2, :);
end

        switch type
            case 'R'
                componentCurrentsTime(idx, :) = voltageDiff / value;
            case 'L'
    % Handle inductor current calculation using Euler's method (discrete integration)
     current = zeros(1, length(voltageDiff));
     for k = 2:length(voltageDiff)
         % Apply the discrete Euler method for integrating the voltage across the inductor
        current(k) = current(k-1) + (voltageDiff(k-1) / value) * timeStep;
     end
     componentCurrentsTime(idx, :) = current;

            case 'C'
              componentCurrentsTime(idx, :) = value * [0 diff(voltageDiff) / timeStep];
        end
    end

    % Plot node voltages
    plotResults(t_vector, nodeVoltagesTime, 'Node Voltages', 'Voltage (V)');

    % Prompt user for component current to plot
    fprintf('Available components to plot current:\n');
    for idx = 1:length(components)
        fprintf('%d: %s between nodes %d and %d\n', idx, components{idx}{3}, components{idx}{1}, components{idx}{2});
    end
    choice = input('Enter the component number to plot its current: ');
    
    % Validate and plot chosen component current
    if choice >= 1 && choice <= length(components)
        plotResults(t_vector, componentCurrentsTime(choice, :), ['Current through Component ', num2str(choice)], 'Current (A)');
    else
        error('Invalid component choice.');
    end

    % Nested function to add component to conductance matrix
    function addComponentToMatrix(admittance, node1, node2)
        if node1 > 0
            G(node1, node1) = G(node1, node1) + admittance;
        end
        if node2 > 0
            G(node2, node2) = G(node2, node2) + admittance;
        end
        if node1 > 0 && node2 > 0
            G(node1, node2) = G(node1, node2) - admittance;
            G(node2, node1) = G(node2, node1) - admittance;
        end
    end
end

function plotResults(timeVector, data, titleText, yLabelText)
    figure('Name', titleText);
    plot(timeVector, data, 'LineWidth', 2);
    title(titleText);
    xlabel('Time (s)');
    ylabel(yLabelText);
    grid on;
end
