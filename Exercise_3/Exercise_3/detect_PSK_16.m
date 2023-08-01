function [est_X, est_bit_seq] = detect_PSK_16(Y)
            
    % Define the constellation points of 16-PSK
    mi = cos((0:15) * 2 * pi / 16);
    mq = sin((0:15) * 2 * pi / 16);

    % Initialize variables
    L = length(Y);
    est_X = zeros(1, L);
    est_bit_seq = zeros(4 * L, 1);

    for i = 1:length(Y)
        Yi = real(Y(i));
        Yq = imag(Y(i));

        % Compute the distances to each constellation point
        distances = sqrt((Yi - mi).^2 + (Yq - mq).^2);

        % Find the index of the closest constellation point
        [~, index] = min(distances);

        % Estimate the symbol and corresponding bits
        est_X(i) = mi(index) + 1i * mq(index);
        est_bit_seq(4*i-3:4*i) = de2bi(index-1, 4, 'left-msb')';
    end
end








