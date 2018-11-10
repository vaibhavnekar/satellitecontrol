function output_double = custom_str2double(input_str)
i = 0.0;
output_double = 0;
npast = 0;
num_sign = 1;
input_str = input_str(input_str ~= ' ');
for i = 1:length(input_str)
    if input_str(i) == '-'
        num_sign = -1;
    elseif input_str(i) == '.' && npast == 0
        npast = 1;
    elseif npast >= 1
        output_double = output_double + (input_str(i) - '0')*10^-npast;
        npast = npast+1;
    else
        % Value is greater than 1
        output_double = output_double*10 + (input_str(i) - '0');
    end
end
output_double = num_sign*output_double;
end