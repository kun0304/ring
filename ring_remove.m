%   ==========
%   The Code is created based on the method described in the following paper 
%   "Computed Tomography Ring Artifact Correction Method with Super-Pixel and Adaptive Relative Total Variation", Na Li.
%   The paper was under review in the journal
%   The code and the algorithm are for non-comercial use only.
%  
%   Na Li (na.li3@gdmu.edu.cn)
%   Date  : 17 Feb 2024
%   Version : 1.0 
%   Copyright 2024, Guangdong Medical University.

function [correct2_img] = ring_remove(ring_img)
    % Determine data type
    dataType = class(ring_img);
    dt = str2func(dataType);
    
    % Convert input image to double
    ring_img_o = double(ring_img);
    totalSlice = size(ring_img,3);

    % Normalize image intensity to [0, 255]
    img_temp = ring_img_o;
    img_temp1 = 255 * (img_temp - min(img_temp(:))) / (max(img_temp(:)) - min(img_temp(:)));

    % Superpixel segmentation
    nnn = 250;
    if size(img_temp, 3) == 1
        [super_ring, Nn] = superpixels(img_temp1, nnn);
    else
        [super_ring, Nn] = superpixels3(img_temp1, nnn);
    end

    % Replace superpixels with mean intensity
    for ii = 1:Nn
        img_temp(super_ring == ii) = mean(img_temp(super_ring == ii));
    end

    % Obtain ring component
    ring_t = ring_img_o - img_temp;
    ring = zeros(size(ring_img));

    % Process each slice
    for i = 1:totalSlice
        ring_i = ring_t(:, :, i);
        [polar_r_I, pol_s] = im_cart2pol(ring_i, size(ring_i), [], 'cubic');
        polar_r_I = repmat(mean(polar_r_I, 2), [1 size(polar_r_I, 2)]);
        ring(:, :, i) = im_pol2cart(polar_r_I, pol_s, size(ring_i));
    end

    % Correct image by subtracting ring component
    correct1_img = ring_img_o - ring;

    % Select slices for further processing
    select_slice_num = 3;
    if totalSlice < select_slice_num
        select_slice = 1:totalSlice;
    else
        select_slice = round(linspace(1, totalSlice, select_slice_num + 2));
        select_slice = select_slice(2:end-1);
    end

    % Initialize lambda array
    lambda = zeros(1, select_slice_num);
    
    % Process selected slices
    for k = 1:numel(select_slice)
        i = select_slice(k);
        r_I = correct1_img(:, :, i);
        r_I_t = imresize(r_I, [144 144]);
        t_max = max(r_I_t(:)); t_min = min(r_I_t(:));
        r_I_t = (r_I_t - t_min) / (t_max - t_min);

        % Superpixel segmentation and smoothing
        [super_ring, Nn] = superpixels(r_I_t, nnn);
        for ii = 1:Nn
            super_ring(super_ring == ii) = mean(r_I_t(super_ring == ii));
        end

        % Optimize lambda
        objectiveFcn = @(lambda)rtv_lambda(lambda, r_I_t, super_ring);
        lb = 0.0035; ub = 0.015; x0 = 0.0065;
        optsm = optimoptions("patternsearch", 'Algorithm', 'nups-mads', 'MaxIterations', 50);
        lambda(k) = patternsearch(objectiveFcn, x0, [], [], [], [], lb, ub, [], optsm);
    end

    % Interpolate lambda values
    if totalSlice > 2
        lambda_n = interp1([1, select_slice, totalSlice], [lambda(1), lambda, lambda(end)], 1:totalSlice);
    else
        lambda_n = lambda;
    end

    % Obtain second ring component
    ring2 = zeros(size(ring_img));
    for i = 1:totalSlice
        r_I = correct1_img(:, :, i);
        t_max = max(r_I(:)); t_min = min(r_I(:));
        rtv_I = tsmooth((r_I - t_min) / (t_max - t_min), lambda_n(i));
        rtv_I = rtv_I * (t_max - t_min) + t_min;

        [residual_polar, pol_s] = im_cart2pol(r_I - rtv_I, size(r_I) * 2, [], 'cubic');
        hsize = [1 round(size(ring_img, 1) / 7)];
        h = fspecial('average', hsize);
        residual = [residual_polar, residual_polar, residual_polar];
        residual = imfilter(residual, h);
        residual = residual(:, size(residual, 2) / 3 + 1:2 * size(residual, 2) / 3);
        ring2(:, :, i) = im_pol2cart(residual, pol_s, size(r_I));
    end

    % Corrected image
    correct2_img = correct1_img - ring2;
    correct2_img = correct2_img - (mean(single(correct2_img(:))) - mean(single(ring_img(:))));
    correct2_img = dt(correct2_img);
end

function smoothed_image = tsmooth(input_image, lambda)
    % Total variation smoothing algorithm.
    % Inputs:
    %   - input_image: the input grayscale image to be smoothed
    %   - lambda: parameter controlling the amount of smoothing
    % Outputs:
    %   - smoothed_image: the smoothed output image
    
    % Initialize parameters
    initial_sigma = 3.0; % Initial value for sigma
    sharpness = 0.02; % Sharpness parameter
    max_iterations = 3; % Maximum number of iterations
    current_image = input_image;
    current_sigma = initial_sigma;
    lambda = lambda / 2.0; % Adjust lambda
    
    % Iterative smoothing process
    for iter = 1:max_iterations
        % Compute texture weights
        [wx, wy] = computeTextureWeights(current_image, current_sigma, sharpness);
        
        % Solve linear equation to update image
        current_image = solveLinearEquation(input_image, wx, wy, lambda);
        
        % Update sigma
        current_sigma = current_sigma / 2.0;
        
        % Ensure sigma doesn't go below a minimum value
        if current_sigma < 0.5
            current_sigma = 0.5;
        end
    end
    
    % Output the smoothed image
    smoothed_image = current_image;
end

function [retx, rety] = computeTextureWeights(fin, sigma, sharpness)
    % computeTextureWeights: Computes texture weights for an input image.
    % Arguments:
    %   - fin: Input image.
    %   - sigma: Standard deviation for Gaussian filter.
    %   - sharpness: Sharpness parameter.
    % Returns:
    %   - retx: Texture weights in x-direction.
    %   - rety: Texture weights in y-direction.

    % Set sharpness threshold
    vareps_s = sharpness;
    % Small constant for numerical stability
    vareps = 0.001;

    % Compute gradients in x and y directions
    fx = diff(fin, 1, 2);
    fx = padarray(fx, [0, 1, 0], 'post');
    fy = diff(fin, 1, 1);
    fy = padarray(fy, [1, 0, 0], 'post');

    % Compute weights based on gradient magnitudes
    wto = max(sum(sqrt(fx.^2 + fy.^2), 3) / size(fin, 3), vareps_s).^(-1);

    % Apply Gaussian filter to input image
    fbin = lpfilter(fin, sigma);
    % Compute gradients of filtered image
    gfx = diff(fbin, 1, 2);
    gfx = padarray(gfx, [0, 1], 'post');
    gfy = diff(fbin, 1, 1);
    gfy = padarray(gfy, [1, 0], 'post');

    % Compute weights based on filtered gradient magnitudes
    wtbx = max(sum(abs(gfx), 3) / size(fin, 3), vareps).^(-1);
    wtby = max(sum(abs(gfy), 3) / size(fin, 3), vareps).^(-1);

    % Combine weights in x and y directions
    retx = wtbx .* wto;
    rety = wtby .* wto;

    % Ensure boundary weights are set to zero
    retx(:, end) = 0;
    rety(end, :) = 0;
end

function result = conv2_sep(image, sigma)
    % Determine the kernel size based on the provided sigma
    kernel_size = round(5 * sigma);
    kernel_size = max(1, bitshift(kernel_size, -1)) * 2 + 1; % Ensure kernel size is odd
    
    % Generate the Gaussian kernel for horizontal convolution
    horizontal_kernel = fspecial('gaussian', [1, kernel_size], sigma);
    
    % Convolve the image with the horizontal Gaussian kernel
    convolved_image = conv2(image, horizontal_kernel, 'same');
    
    % Generate the Gaussian kernel for vertical convolution
    vertical_kernel = fspecial('gaussian', [kernel_size, 1], sigma);
    
    % Convolve the horizontally convolved image with the vertical Gaussian kernel
    result = conv2(convolved_image, vertical_kernel, 'same');
end

function filtered_image = lpfilter(input_image, sigma)
    % Initialize the filtered image with the input image
    filtered_image = input_image;
    
    % Apply separable Gaussian filtering to each channel of the input image
    for channel = 1:size(filtered_image, 3)
        filtered_image(:, :, channel) = conv2_sep(input_image(:, :, channel), sigma);
    end
end

function OUT = solveLinearEquation(IN, wx, wy, lambda)
% Solve a linear equation system for each channel of the input image
% using the conjugate gradient method with preconditioning.

% Extract dimensions of the input image
[r, c, ch] = size(IN);
k = r * c;

% Compute the weighted differences in x and y directions
dx = -lambda * wx(:);
dy = -lambda * wy(:);

% Construct the sparse matrix A for the linear system
B = zeros(size(dx, 1), 2);
B(:, 1) = dx;
B(:, 2) = dy;
d = [-r, -1];
A = spdiags(B, d, k, k);

% Compute the east, west, south, and north terms for the diagonal of A
e = dx;
w = padarray(dx, r, 'pre');
w = w(1:end - r);
s = dy;
n = padarray(dy, 1, 'pre');
n = n(1:end - 1);
D = 1 - (e + w + s + n);

% Incorporate the diagonal terms into the sparse matrix A
A = A + A' + spdiags(D, 0, k, k);

% Preconditioning using incomplete Cholesky factorization
L = ichol(A, struct('michol', 'on'));

% Initialize the output matrix with the same size as the input image
OUT = IN;

% Solve the linear equation system for each channel of the input image
for ii = 1:ch
    tin = IN(:, :, ii);
    [tout, ~] = pcg(A, tin(:), 0.1, 100, L, L');
    OUT(:, :, ii) = reshape(tout, r, c);
end

end

function y = inpaintn(x, n, y0, m)
    % Convert input to double
    x = double(x);
    
    % Set default value for 'n' if not provided
    if nargin == 1 || isempty(n)
        n = 100;
    end

    % Get size and dimensions of input
    sizx = size(x);
    d = ndims(x);
    
    % Compute Lambda matrix
    Lambda = zeros(sizx);
    for i = 1:d
        siz0 = ones(1, d);
        siz0(i) = sizx(i);
        Lambda = bsxfun(@plus, Lambda, cos(pi * (reshape(1:sizx(i), siz0) - 1) / sizx(i)));
    end
    Lambda = 2 * (d - Lambda);

    % Identify finite elements in x
    W = isfinite(x);
    
    % Set initial guess based on y0 or calculate it
    if nargin == 3 && ~isempty(y0)
        y = y0;
        s0 = 3;
    else
        if any(~W(:))
            [y, s0] = InitialGuess(x, isfinite(x));
        else
            y = x;
            return
        end
    end
    x(~W) = 0;

    % Set default value for 'n' if not provided
    if isempty(n) || n <= 0
        n = 100;
    end
    
    % Generate logarithmically spaced values for iteration
    s = logspace(s0, -6, n);

    % Set RF value
    RF = 2;

    % Set default value for 'm' if not provided
    if nargin < 4 || isempty(m)
        m = 2;
    end
    Lambda = Lambda .^ m;

    % Perform iterative inpainting
    for i = 1:n
        Gamma = 1 ./ (1 + s(i) * Lambda);
        y = RF * idctn(Gamma .* dctn(W .* (x - y) + y)) + (1 - RF) * y;
    end

    % Replace inpainted elements in y with original values from x
    y(W) = x(W);
end

function [z, s0] = InitialGuess(y, I)
    % Check if Image Processing Toolbox is available
    if license('test', 'image_toolbox')
        [~, L] = bwdist(I);
        z = y;
        z(~I) = y(L(~I));
        s0 = 3;
    else
        % Warn if Image Processing Toolbox is not available
        warning('MATLAB:inpaintn:InitialGuess', ...
            ['BWDIST (Image Processing Toolbox) does not exist. ', ...
            'The initial guess may not be optimal; additional', ...
            ' iterations can thus be required to ensure complete', ...
            ' convergence. Increase N value if necessary.']);
        z = y;
        z(~I) = mean(y(I));
        s0 = 6;
    end
end

% Define forward Discrete Cosine Transform
function y = dctn(y)
    % Perform DCT along each dimension
    y = double(y);
    sizy = size(y);
    y = squeeze(y);
    dimy = ndims(y);
    
    % Check if y is a vector
    if isvector(y)
        dimy = 1;
        if size(y, 1) == 1
            y = y.';
        end
    end

    % Generate exponential weights for DCT
    w = cell(1, dimy);
    for dim = 1:dimy
        n = (dimy == 1) * numel(y) + (dimy > 1) * sizy(dim);
        w{dim} = exp(1i * (0:n-1).' * pi / 2 / n);
    end

    % Compute DCT along each dimension
    if ~isreal(y)
        y = complex(dctn(real(y)), dctn(imag(y)));
    else
        for dim = 1:dimy
            siz = size(y);
            n = siz(1);
            y = y([1:2:n, 2 * floor(n / 2):-2:2], :);
            y = reshape(y, n, []);
            y = y * sqrt(2 * n);
            y = ifft(y, [], 1);
            y = bsxfun(@times, y, w{dim});
            y = real(y);
            y(1, :) = y(1, :) / sqrt(2);
            y = reshape(y, siz);
            y = shiftdim(y, 1);
        end
    end

    y = reshape(y, sizy);
end

% Define inverse Discrete Cosine Transform
function y = idctn(y)
    % Perform IDCT along each dimension
    y = double(y);
    sizy = size(y);
    y = squeeze(y);
    dimy = ndims(y);
    
    % Check if y is a vector
    if isvector(y)
        dimy = 1;
        if size(y, 1) == 1
            y = y.';
        end
    end

    % Generate exponential weights for IDCT
    w = cell(1, dimy);
    for dim = 1:dimy
        n = (dimy == 1) * numel(y) + (dimy > 1) * sizy(dim);
        w{dim} = exp(1i * (0:n-1).' * pi / 2 / n);
    end

    % Compute IDCT along each dimension
    if ~isreal(y)
        y = complex(idctn(real(y)), idctn(imag(y)));
    else
        for dim = 1:dimy
            siz = size(y);
            n = siz(1);
            y = reshape(y, n, []);
            y = bsxfun(@times, y, w{dim});
            y(1, :) = y(1, :) / sqrt(2);
            y = ifft(y, [], 1);
            y = real(y * sqrt(2 * n));
            I = (1:n) * 0.5 + 0.5;
            I(2:2:end) = n - I(1:2:end-1) + 1;
            y = y(I, :);
            y = reshape(y, siz);
            y = shiftdim(y, 1);
        end
    end

    y = reshape(y, sizy);
end

function im_cart = im_pol2cart(im_pol, sz_pol, sz_cart)
    % Compute the origin of the Cartesian image
    origin = sz_cart ./ 2;

    % Compute the origin in Cartesian coordinates
    origin_c = origin;

    % Compute the maximum radius and theta in Cartesian coordinates
    [r_max_c, t_max_c] = compute_max_c(origin_c, sz_cart);

    % Generate the x and y basis for Cartesian coordinates
    x_basis_i = 1 : sz_cart(2);
    x_basis_c = (x_basis_i - origin_c(2)) ./ r_max_c;
    y_basis_i = 1 : sz_cart(1);
    y_basis_c = (y_basis_i - origin_c(1)) ./ r_max_c;
    [x_c, y_c] = meshgrid(x_basis_c, y_basis_c);

    % Generate the theta and radius basis for polar coordinates
    t_basis_c = 1 : sz_pol(2);
    r_basis_c = 1 : sz_pol(1);

    % Generate the meshgrid for polar coordinates and resize it
    [r_c, t_c] = meshgrid(r_basis_c, t_basis_c);
    r_c = imresize(r_c,[size(im_pol,2) size(im_pol,1)]);
    t_c = imresize(t_c,[size(im_pol,2) size(im_pol,1)]);

    % Compute conversion factors from image indices to Cartesian coordinates
    [r_i_to_c, t_i_to_c] = compute_i_to_c(r_max_c, t_max_c);

    % Compute the query points in Cartesian coordinates
    r_query_i = sqrt(x_c.^2 + y_c.^2);
    r_query_c = r_query_i ./ r_i_to_c;
    t_query_i = atan2(y_c, x_c) + pi;
    t_query_c = t_query_i ./ t_i_to_c;

    % Interpolate the polar image to Cartesian coordinates
    im_cart = griddata(r_c, t_c, double(im_pol)', r_query_c, t_query_c);

    % Fix the output type
    im_cart = fix_output_type(im_pol, im_cart, 0.5);

    % Rotate the Cartesian image by 180 degrees
    im_cart = imrotate(im_cart, 180);

    % Inpaint missing values in a specific region of the Cartesian image
    im_cart(round(size(im_cart,1)/2)-5:round(size(im_cart,1)/2)+5,round(size(im_cart,1)/2)-5:end) =...
        inpaintn(im_cart(round(size(im_cart,1)/2)-5:round(size(im_cart,1)/2)+5,round(size(im_cart,1)/2)-5:end), 10);

    % Set NaN and Inf values to zero
    im_cart(isnan(im_cart)) = 0;
    im_cart(isinf(im_cart)) = 0;

    % Translate the Cartesian image by [-1, -1]
    im_cart = imtranslate(im_cart, [-1, -1]);
end

function [im_pol, sz_pol] = im_cart2pol(im_cart, polar_size, origin, interp_method)
    % Compute the maximum and minimum values of the Cartesian image
    im_max = max(single(im_cart(:)));
    im_min = min(single(im_cart(:)));

    % Set default values if not provided
    if nargin < 3
        origin = [];
    end
    if nargin < 4
        interp_method = "makima";
    end

    % Get the size of the Cartesian image
    sz_cart = size(im_cart);

    % Compute the origin if not provided
    if isempty(origin)
        origin = sz_cart ./ 2;
    end

    % Convert interpolation method to string
    interp_method = string(interp_method);

    % Compute the origin in Cartesian coordinates
    origin_c = origin;

    % Generate the x and y basis for Cartesian coordinates
    x_basis_c = 1 : sz_cart(2);
    y_basis_c = 1 : sz_cart(1);
    [x_c, y_c] = meshgrid(x_basis_c, y_basis_c);

    % Compute the maximum radius and theta in Cartesian coordinates
    [r_max_c, t_max_c] = compute_max_c(origin_c, sz_cart);

    % Compute the size of the polar image
    sz_pol = ceil([r_max_c, t_max_c]);

    % Compute conversion factors from image indices to Cartesian coordinates
    [r_i_to_c, t_i_to_c] = compute_i_to_c(r_max_c, t_max_c);

    % Generate the theta and radius basis for polar coordinates
    t_basis_i = 1 : sz_pol(2);
    t_basis_c = t_basis_i .* t_i_to_c;
    r_basis_i = 1 : sz_pol(1);
    r_basis_c = r_basis_i .* r_i_to_c;

    % Generate the meshgrid for polar coordinates
    [t_c, r_c] = meshgrid(t_basis_c, r_basis_c);

    % Compute the query points in Cartesian coordinates
    x_query_c = r_c .* cos(t_c) .* r_max_c + origin_c(2);
    y_query_c = r_c .* sin(t_c) .* r_max_c + origin_c(1);
    x_query_c = imresize(x_query_c, polar_size);
    y_query_c = imresize(y_query_c, polar_size);

    % Interpolate the Cartesian image to polar coordinates
    im_pol = interp2(x_c, y_c, double(im_cart), x_query_c, y_query_c, interp_method);

    % Fix the output type
    im_pol = fix_output_type(im_cart, im_pol, 0.5);

    % Clip the polar image to the maximum and minimum values
    im_pol(im_pol > im_max) = im_max;
    im_pol(im_pol < im_min) = im_min;
end

function im_out = fix_output_type(im_in, im_out, logical_threshold)
    % Convert the output type based on the input type
    t = string(class(im_in));
    if t == "logical"
        im_out = im_out > logical_threshold;
    else
        im_out = cast(im_out, t);
    end
end

function [r_max_c, t_max_c] = compute_max_c(origin_c, sz_cart)
    % Compute the maximum radius and theta in Cartesian coordinates
    corners_c = [...
        1, 1; ...
        sz_cart(1), 1; ...
        1, sz_cart(2); ...
        sz_cart ...
    ];
    sums = sum((origin_c - corners_c) .^ 2, 2);
    r_max_c = sqrt(max(sums));
    t_max_c = 2 * pi * r_max_c;
end

function [r_i_to_c, t_i_to_c] = compute_i_to_c(r_max_c, t_max_c)
    % Compute conversion factors from image indices to Cartesian coordinates
    r_i_to_c = 1 ./ r_max_c;
    t_i_to_c = 2 .* pi ./ t_max_c;
end

function y = rtv_lambda(lambda, r_I, super_ring)
rtv_I = tsmooth(r_I, lambda);
y = abs(rtv_I - super_ring);
y = sum(y(:));
end
