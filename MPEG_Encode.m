function [readyForTransmittion,dict,quant_level,L_Stream, frameRate, L, W] = MPEG_Encode(videoSrc,quant_level, numOfFrames)
switch nargin
  case 0
    videoSrc = "xylophone.mp4";
    quant_level = 1;
    numOfFrames = -1;
  case 1
    quant_level = 1;
    numOfFrames = -1;
  case 2
     numOfFrames = -1;
  case 3
  otherwise
    error('3 inputs only are accepted.')
end
% Read
vidobj = VideoReader(videoSrc);
if(numOfFrames == -1)
    sz = vidobj.NumFrames; % number of frames
else
    sz = numOfFrames;
end
frameRate = vidobj.FrameRate;
% Encoder
outputHuff = {};
outputMotion = {};
finalOutput = {}; %combine both cells {Huff1, MV1, Huff2, MV2, ...}
for i = 1 : sz
    %encoder
    current = read(vidobj, i); % Get current frame
    current = fixLW(double(rgb2gray(current)), 8); % Change it into Gray
    [L, W] = size(current);
    if (rem(i - 1, 10) == 0) % I frame is regenerated each 10 frames
        reference = current;
        % First we need to JPEG encode it to get the quantized error
        % version
        rlc_array = JPEG_Encoder_No_Huff(reference, quant_level);
        output = DCT_Quant(reference, 1);
        reference = IDCT_Quant(output, 1); % set a new reference for P frames
        % Send null in motion vectors
        outputMotion{end + 1} = [0 0];
    else
        motion_vectors = motion_estimation(reference, current, 16, 4);
        prediction = motion_compensation(reference, motion_vectors, 16);
        diff = current - prediction;
        rlc_array = JPEG_Encoder_No_Huff(diff, quant_level); %%%%%%%%%%% <======= shelt diff 7atet curr
        % Send motion vectors as huffman
         relMVs = cell2mat(motion_vectors);
%         cod = huffmanenco(relMVs, dict);
        outputMotion{end + 1} = relMVs;
    end
    % Push the frame to the set to be transmitted later
    outputHuff{end + 1} = rlc_array;

    finalOutput{end+1} = rlc_array;
    finalOutput{end+1} = outputMotion{end};
end

% Apply huffman on all frames
AllBlocks = cell2mat(finalOutput);
% Huffman Encoder
% 1. get the proba
[counts,elements] = groupcounts(AllBlocks');
proba = counts ./ sum(counts);
% 2. construct dict
dict = myHuffDict(elements,proba);
% Apply on each cell in outputHuff (The reason we don't apply on all the
% frames is to be able to reconstruct it in the decoder)
calcHuff = @(x) huffmanenco(x, dict);
finalOutput = cellfun(calcHuff, finalOutput, 'UniformOutput',false);
% Outputs
readyForTransmittion = cell2mat(finalOutput);
L_Stream = cellfun('length',finalOutput);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Other Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tmp = fixLW(Image, n)
% Fix 8's block and complete rest by 0's
[L, W] = size(Image);
originalR = L;
originalC = W;
% Get the first number divisible by 8 in rows and cols
while(true)
    if(mod(L,n) ~= 0)
        L = L+1;
    else
        break;
    end
end
while(true)
    if(mod(W,n) ~= 0)
        W = W+1;
    else
        break;
    end
end
% pad the image with zeros to be divisible by 8
tmp = zeros(L,W);
tmp(1:originalR, 1:originalC) = Image;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function array_rlc = JPEG_Encoder_No_Huff(Image, quant_level)

% Fix 8's block and complete rest by 0's
[L, W] = size(Image);
originalR = L;
originalC = W;
% Get the first number divisible by 8 in rows and cols
while(true)
    if(mod(L,8) ~= 0)
        L = L+1;
    else
        break;
    end
end
while(true)
    if(mod(W,8) ~= 0)
        W = W+1;
    else
        break;
    end
end
% pad the image with zeros to be divisible by 8
tmp = zeros(L,W);
tmp(1:originalR, 1:originalC) = Image;
Image = tmp;

%%%%%%%%%%%%%%%%  Start Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

block_no = 1; %the number of block I'm working on right now
array = zeros(1, L * W); %1-D array of the whole image after DCT and quantization

for r = 1 : 8 : L
for c = 1 : 8 : W
    Block = double(Image(r : r + 7, c : c + 7)); %capturing an 8x8 block
    Block_DCT = Block8_DCT(Block); %applying DCT
    Block_quant = round(Block_DCT./QuantizationTable(quant_level)); %quantizing with the chosen level
    array((block_no - 1) * 64 + 1 : block_no * 64) = serpentine(Block_quant); %applying the zigzag method
    block_no = block_no + 1;
end
end

array_rlc = run_length_encoding(array); %apply run length encoding on the whole 1-D image array
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = DCT_Quant(Image, quant_level)
% Fix 8's block and complete rest by 0's
[L, W] = size(Image);
originalR = L;
originalC = W;
% Get the first number divisible by 8 in rows and cols
while(true)
    if(mod(L,8) ~= 0)
        L = L+1;
    else
        break;
    end
end
while(true)
    if(mod(W,8) ~= 0)
        W = W+1;
    else
        break;
    end
end
% pad the image with zeros to be divisible by 8
tmp = zeros(L,W);
tmp(1:originalR, 1:originalC) = Image;
Image = tmp;

%%%%%%%%%%%%%%%%  Start Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1 : 8 : L
for c = 1 : 8 : W
    Block = double(Image(r : r + 7, c : c + 7)); %capturing an 8x8 block
    Block_DCT = Block8_DCT(Block); %applying DCT
    Block_quant = round(Block_DCT./QuantizationTable(quant_level)); %quantizing with the chosen level
    output(r : r + 7, c : c + 7) = Block_quant;
end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DCT_result = Block8_DCT(block)
basis = zeros(8); %basis matrix
DCT_result = zeros(8); %final DCT result matrix
for u = 0 : 7
    for v = 0 : 7
        for x = 0 : 7
            for y = 0 : 7
                basis(x + 1, y + 1) = cos((2 * x + 1) * u * pi / 16) * cos((2 * y + 1) * v * pi / 16); %constructing the basis matrix(u, v)
            end
        end
        DCT_result(u + 1, v + 1) = sum(sum(basis.*block)); %summing the multiplication of the basis with the block into one entry in the DCT result matrix
    end
end

%scaling the DCT result matrix
DCT_result(1, :) = DCT_result(1, :) / 2;
DCT_result(:, 1) = DCT_result(:, 1) / 2;
DCT_result(:, :) = DCT_result(:, :) / 16;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Table = QuantizationTable(x)
if x==1  % for Low compression
    Table = [1 1 1 1 1 2 2 4
        1 1 1 1 1 2 2 4
        1 1 1 1 2 2 2 4
        1 1 1 1 2 2 4 8
        1 1 2 2 2 2 4 8
        2 2 2 2 2 4 8 8
        2 2 2 4 4 8 8 16
        4 4 4 4 8 8 16 16];
elseif x==2  % for high compression
    Table=[1 2 4 8 16 32 64 128;
        2 4 4 8 16 32 64 128;
        4 4 8 16 32 64 128 128;
        8 8 16 32 64 128 128 256;
        16 16 32 64 128 128 256 256;
        32 32 64 128 128 256 256 256;
        64 64 128 128 256 256 256 256;
        128 128 128 256 256 256 256 256];
else
    disp('Error in choose X');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = IDCT_Quant(Image, quant_level)
[L, W] = size(Image);

for r = 1 : 8 : L
for c = 1 : 8 : W
    block_deQuant = round(Image(r : r + 7, c : c + 7).*QuantizationTable(quant_level)); %dequantizing the 8x8 block
    block = round(Block8_IDCT(block_deQuant)); %applying inverse DCT
    output(r : r + 7, c : c + 7) = block; %adding the retrieved block back into the image
end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IDCT_result = Block8_IDCT(block)
basis = zeros(8); %basis matrix
IDCT_result = zeros(8); %final IDCT result matrix
for u = 0 : 7
    for v = 0 : 7
        for x = 0 : 7
            for y = 0 : 7
                basis(x + 1, y + 1) = cos((2 * x + 1) * u * pi / 16) * cos((2 * y + 1) * v * pi / 16); %constructing the basis matrix(u, v)
            end
        end
        IDCT_result = IDCT_result + basis*block(u + 1, v + 1); %appending the output of scaled basis to the final result
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [motion_vectors] = motion_estimation(reference, current, block_size, search_size)

%block_size is 16*16, resolution is one pixel (no subpixel estimation)
[L, W] = size(current);
motion_vectors = {};

%i and j are the top left corners of the block
for i = 1 : block_size : L
for j = 1 : block_size : W
    if (i + block_size - 1 > L || j + block_size - 1 > W)
        continue;
    end
    block = current(i : i + block_size - 1, j : j + block_size - 1); %check lw b3ady boundries fe akher itteration
    margin = search_size/2; %search size = 4
    max_corr = 0;
    %r and c are motion vector for single block

    for m1 = -margin : margin
    for m2 = -margin : margin
        x1 = i + m1;
        y1 = j + m2;
        x2 = x1 + block_size - 1;
        y2 = y1 + block_size - 1;
        if x1 < 1 || x2 > L || y1 < 1 || y2 > W
            continue;
        end
        corr = correlation(block, reference(x1 : x2, y1 : y2));
        if corr > max_corr
            r = m1;
            c = m2;
            max_corr = corr;
        end
    end
    end
    motion_vectors{end + 1} =[r c];
end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corr = correlation(block1, block2)

[L, W] = size(block1);
corr = 0;

for i = 1 : L
for j = 1 : W
    corr = corr + (block1(i, j) * block2(i, j));
end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prediction = motion_compensation(reference, motion_vectors, block_size)

[L, W] = size(reference);
prediction = zeros([L, W]);

idx = 1;
for i = 1 : block_size : L
for j = 1 : block_size : W
    x1 = i;
    x2 = i + block_size - 1;
    y1 = j;
    y2 = j + block_size - 1;
    if(idx > length(motion_vectors))
        continue;
    end
    motion_vector = cell2mat(motion_vectors(idx));
    r = motion_vector(1);
    c = motion_vector(2);
    if (x2 > L || y2 > W || x2 + r > L || y2 + c > W)
        continue;
    end
    prediction(x1 : x2, y1 : y2) = reference(x1 + r : x2 + r, y1 + c : y2 + c); %check lw b3ady boundries fe akher itteration
    idx = idx + 1;
end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function array = serpentine(block)
[L, W] = size(block);
if L ~= W %check if the block is square
    error("Hal anta 3abet?!");
end
array = zeros(1, L * W); %the output 1-D array
loops = L + W - 1; %number of loops (diagonal moves) for the block

x = 1; %index to hold the row
y = 1; %index to hold the column
idx = 1; %index of the output array

%I have constructed this algorithm based on the tracing I did in "Serpentine algo tracing.txt"
for i = 1 : loops
    if ceil(i / 2) == i / 2 %if loop is even, we start with x as minimum and y as maximum
            x = 0;
            y = i + 1;
            even = 1; % a flag to tell whether the loop is even or odd
    else                    %if loop is odd, we start with x as maximum and y as minimum
            x = i + 1;
            y = 0;
            even = 0;
    end
    for j = 1 : i %within each loop, we iterate by its number
        if even == 1    %if loop is even, we increment the x and decrement y
            x = x + 1;
            y = y - 1;
        else            %if loop is even, we decrement the x and increment y
            x = x - 1;
            y = y + 1;
        end
        if (x > L || y > W) %we did not go out of bound
            continue;
        end
        array(idx) = block(x, y); %assign the corresponding pixel in its place in the 1-D array
        idx = idx + 1;  %increment the index of the 1-D array
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = run_length_encoding(input)
N = length(input);

c = 0; % counts consecutive zeros
idx = 1;
for i = 1 : N
    if input(i) ~= 0 %if the value is not zero
        if c ~= 0 %replace all counted zeros with 0 and the number of zeros
            output(idx) = 0;
            idx = idx + 1;
            output(idx) = c;
            idx = idx + 1;
            c = 0; %reset number of zeros to zero
        end
        output(idx) = input(i); %in case not zero, copy the value as it is
        idx = idx + 1;
    else %if the value is zero, count consecutive zeros
        c = c + 1;
    end
end

%in case, we have been counting zeros but the loop has ended
if c ~= 0 %replace all zeros with 0 and the number of zeros
            output(idx) = 0;
            idx = idx + 1;
            output(idx) = c;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = myHuffDict(symbolsVector,probaVector)
% The documentation for this function is attached in
% "A1_myHuff_Function_Documentation.pdf"

% First we crate an index for each symbol to track it
symbolIndex = 1:length(symbolsVector);
dict = table( symbolIndex', symbolsVector, probaVector, strings([length(symbolsVector), 1]));
output = table(symbolIndex', symbolsVector, strings([length(symbolsVector), 1]), ...
    'VariableNames',["index", "Symbol","Code"]);
for i = length(symbolsVector) :-1: 2
    % 1. sort
    dict = sortrows(dict, "probaVector", "descend");
    % 2. Detect the index of the last 2 elements
    index1 = dict{i,1};
    index2 = dict{i-1,1};
    % 3. Calculate the new probability
    newProba = dict{i,3} + dict{i-1, 3};
    % 4. Add coding bits:
    % 5.Add bit '1' for the index i 
    if(dict{i, 1} == -1)
        trackArray = str2num(dict{i,4});
        for j = 1 : length(trackArray)
            opIndex1 = find(ismember(output.index, trackArray(j), "rows"), 1, 'first');
            output(opIndex1,3) = num2cell(strcat('1',output{opIndex1,3}));
        end
    else
        opIndex1 = find(ismember(output.index, index1, "rows"), 1, 'first');
        output(opIndex1,3) = num2cell(strcat('1',output{opIndex1,3}));
    end
    % 6. Add bit '1' for the index i-1
    if(dict{i-1, 1} == -1)
        trackArray = str2num(dict{i-1,4});
        for j = 1 : length(trackArray)
            opIndex2 = find(ismember(output.index, trackArray(j), "rows"), 1, 'first');
            output(opIndex2,3) = num2cell(strcat('0',output{opIndex2,3}));
        end
    else
        opIndex2 = find(ismember(output.index, index2, "rows"), 1, 'first');
        output(opIndex2,3) = num2cell(strcat('0',output{opIndex2,3}));
    end
    % 7. Update the trackArray (queue)
    if(dict{i-1, 4} == "" && index1 ~= -1 && index2 ~= -1)
        dict{i-1, 4} = {sprintf("%d %d", index1, index2)};
    else
        if(index1 == -1 && index2 ~= -1)
            dict{i-1, 4} = {strcat(dict{i, 4}, " ", int2str(index2))};
        elseif(index2 == -1 && index1 ~= -1)
            dict{i-1, 4} = {strcat(dict{i-1, 4}, " ", int2str(index1))};
        else
            dict{i-1, 4} = {strcat(dict{i-1, 4}, " ", dict{i, 4})};
        end
    end
    % 8. Update the dict table
    dict(i-1, 3) = num2cell(newProba);
    dict(i-1, 1) = num2cell(-1);
    dict(i,:) = [];
end

out = output(:,2:3);
out = table2cell(out);
for i = 1 : length(out)
    mystr = convertStringsToChars(out{i,2});
    x = mat2cell(mystr, 1, ones(1,numel(mystr)));
    outy = cell2mat(cellfun(@str2num, x, 'uniform', 0));
    out{i,2} = outy;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = myHuffEnc(input,dict)
encodedString = [];
for i = 1 : length(input)
    index = find(ismember(dict.Symbol, input(i)),1,'first');
    bits = convertStringsToChars(dict{index, 2});
    encodedString = [encodedString bits];
end
x = mat2cell(encodedString, 1, ones(1,numel(encodedString)));
output = cell2mat(cellfun(@str2num, x, 'uniform', 0));
end





