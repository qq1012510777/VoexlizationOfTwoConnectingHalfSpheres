function f = Num2Str_Set_Width(n, width_x)
    kl = num2str(n, '%32.16e');

    kl_1 = kl(size(kl, 2) - 2:size(kl, 2));

    exponent = str2num(kl_1);

    if (exponent == 0)
        f = num2str(n);
        return
    end

    kl_2 = kl(1:size(kl, 2) - 4);

    real_part = str2num(kl_2);

    if (abs(exponent) >= 100)
        error('The exponent is too large!');
    end

    sign__exponent = [];

    if (exponent > 0)
        sign__exponent = '+';
    else
        sign__exponent = '-';
    end

    CharArray_exponent = [sign__exponent, num2str(abs(exponent))];

    remainning_width = width_x - size(CharArray_exponent, 2);

    if (real_part < 0)
        remainning_width = remainning_width - 2;
        % real_part = round(real_part, remainning_width - 1);

        % f = ['-', num2str(abs(real_part)), CharArray_exponent];

        if (remainning_width > size(kl_2, 2))
            f = [kl_2([1:end]), CharArray_exponent];
        else
            f = [kl_2([1:remainning_width]), CharArray_exponent];
        end

    else
        remainning_width = remainning_width - 1;
        % real_part = round(real_part, remainning_width - 1);

        % f = [num2str(abs(real_part)), CharArray_exponent];
        if (remainning_width > size(kl_2, 2))
            f = [kl_2([1:end]), CharArray_exponent];
        else
            f = [kl_2([1:remainning_width]), CharArray_exponent];
        end

    end

end
