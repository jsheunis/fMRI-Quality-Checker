function f = displayMaskContour(epi_img, mask, rotate, Nskip)

[Ni, Nj, Nk] = size(epi_img);

if rotate
    for p = 1:Nk
        epi_img(:,:,p) = rot90(epi_img(:,:,p));
    end
end

bound = cell(Nk,1);
for k = 1:Nk
    bound{k,1} = bwboundaries(squeeze(mask(:, :, k)));
end


f = figure;
im = cell(9,1);
for r = 1:3
    for c = 1:3
        el = 3*(r-1)+c;
        subplot(3,3,el);
        im{el,1} = imagesc(epi_img(:, :, Nskip*el));
        colormap gray;
        hold on;
        bmax = numel(bound{Nskip*el, 1});
        for b = 1:bmax
            ax = plot(bound{Nskip*el,1}{b,1}(:,2), bound{Nskip*el,1}{b,1}(:,1), 'r', 'LineWidth', 2.5);
            if rotate
                direction = [0 0 1];
                rotate(ax,direction,90);
            end
        end
        hold off;
    end
end

subplot(3,3,2); title('Brain mask contours overlayed on mean EPI')
