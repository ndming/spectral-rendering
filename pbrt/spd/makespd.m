function makespd(file)
    arguments
        file {mustBeFile}
    end

    I = readmatrix(file);
    [m, ~] = size(I);

    spd = zeros(1, 2 * m);
    spd(1:2:end) = I(:,1);
    spd(2:2:end) = I(:,2);

    match = wildcardPattern + ".";
    out = extract(file, match) + "txt";
    rep = extract(file, match) + "spd";
    writematrix(spd, out, Delimiter=" ");
    movefile(out, rep);
end