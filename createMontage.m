function montage = createMontage(img, columns, rotate, str, clrmp, mask)

montage = struct;

[Ni, Nj, Nk] = size(img);

if rotate
    for p = 1:Nk
        img(:,:,p) = rot90(img(:,:,p));
    end
end


filler = zeros(Ni, Nj);
rows = floor(Nk/columns);
fill = mod(Nk, columns);
if fill == 0
    N_fill = 0;
else
    N_fill = columns - mod(Nk, columns);
end

montage.rows = rows;
montage.columns = columns;
montage.N_fill = N_fill;

assignin('base', 'montagexxx', montage)

parts = {};

for i = 1:rows
    for j = 1:columns
        if j ==1
            parts{i} = img(:,:,(columns*(i-1)+j));
        else
            parts{i} = cat(2, parts{i}, img(:,:,(columns*(i-1)+j)));
        end
    end
    if i ==1
        whole = parts{i};
    else
        whole = cat(1, whole, parts{i});
    end
end

if N_fill ~= 0
    % last row
    last_parts = img(:,:,(rows*columns+1));
    for k = (rows*columns+2):Nk
        last_parts = cat(2, last_parts, img(:,:,k));
    end
    for m = 1:N_fill
        last_parts = cat(2, last_parts, filler);
    end
    montage.image = cat(1, whole, last_parts);
else
    montage.image = whole;
end



f = figure;imagesc(montage.image); colormap(clrmp); colorbar;
title(str);
montage.f = f;





