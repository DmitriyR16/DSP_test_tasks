function symbols = Modulation(bits, mod_type)
% This function transforms bit sequence into modulation symbols
%
% Input:
% 
% bits - a vector of input bit sequence 
% mod_type - a string with the name of modulation type (QPSK, 16QAM and 64QAM are implemented)
%
% Output:
% 
% symbols - an output sequence of constellation points 

    switch mod_type
        case 'QPSK'
            symbols = (-1).^bits(1:2:end) + 1i*(-1).^bits(2:2:end);
        case '16QAM'
            symbols = 3 * (-1).^(bits(1:4:end)+1) .* (1/3).^bits(2:4:end) + 3i * (-1).^(bits(3:4:end)+1) .* (1/3).^bits(4:4:end);
        case '64QAM'
            bit_in_sym = 6;
            symbols = zeros(1, length(bits)/bit_in_sym); 
            for i = 1:length(symbols)
                bit_word = bits(1+(i-1)*bit_in_sym:i*bit_in_sym);
                inphase_bits = bit_word(5:6);
                inphase_bits_str = num2str(inphase_bits);
                inphase_dec = bin2dec(inphase_bits_str);
                inphase_factor = 1;
                switch inphase_dec
                    case 0
                        inphase_factor = 7;
                    case 1
                        inphase_factor = 5;
                    case 3
                        inphase_factor = 3;
                end
                quadr_bits = bit_word(2:3);
                quadr_bits_str = num2str(quadr_bits);
                quadr_dec = bin2dec(quadr_bits_str);                
                quadr_factor = 1;
                switch quadr_dec
                    case 0
                        quadr_factor = 7;
                    case 1
                        quadr_factor = 5;
                    case 3
                        quadr_factor = 3; 
                end
                symbols(i) = (-1)^(bit_word(1)+1) * quadr_factor + 1j * (-1)^(bit_word(4)+1) * inphase_factor;
           end
    end
end
