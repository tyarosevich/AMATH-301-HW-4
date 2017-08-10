function dy = chaos_ode(t, y)
a = 16; b = 5; c = 10; d = 6; e = 18;f = .05;

dy = [
    -a * y(1,:) + f * y(2,:) .* y(3,:);
    c * y(2,:) - d * y(1,:) .* y(3,:);
    -b * y(3,:) + e * y(2,:).^2];

