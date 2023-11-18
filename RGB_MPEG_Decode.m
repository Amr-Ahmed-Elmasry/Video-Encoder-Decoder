function  RGB_MPEG_Decode(code,dict,quant_level,L_Stream, fps, L, W, outputName)
readyForTransmittion = code;
outputVideo = VideoWriter(fullfile(outputName), "MPEG-4");
outputVideo.FrameRate = fps;
open(outputVideo);

outputRGB = zeros(L,W,3);
outputRGBIndex = 1; % 1, 2, 3 => R, G, B

skipParam = 2*L_Stream(2); % The first frame is guareented to be I then the next is null
counter = 1;
rgbcounter = 1;

for i = 1 : 2 : length(L_Stream)
    writenowFlag = (rem(rgbcounter, 3) == 0);
    n = L_Stream(i); % Length of this stream
    m = L_Stream(i+1);

    % Get the stream and remove the bits
    stream = readyForTransmittion(1:n);
    stream = huffmandeco(stream, dict);
    readyForTransmittion = [readyForTransmittion(n+1:end)];

    % Get the MV and remove its bits
    mV = readyForTransmittion(1:m);
    readyForTransmittion = [readyForTransmittion(m+1:end)];
    % decode MV if != 1
    if (m > skipParam)
        decodd = huffmandeco(mV, dict);
        I = ones(1,length(decodd)/2)*2;
        motion_vectors = mat2cell(decodd,1,I);
    end

    if (rem(counter - 1, 10) == 0)
        reference_back = JPEG_Decoder_No_Huff(stream, quant_level, L, W);
        if(writenowFlag)
            outputRGB(:,:, outputRGBIndex) = uint8(reference_back);
            outputRGBIndex = 1;

            writeVideo(outputVideo,uint8(outputRGB));

            outputRGB = zeros(L,W,3);
        else
            outputRGB(:,:, outputRGBIndex) = uint8(reference_back);
            outputRGBIndex = outputRGBIndex + 1;
        end
    else
        diff_back = JPEG_Decoder_No_Huff(stream, quant_level, L, W);
        %using motion vectors get the prediction
        prediction_back = motion_compensation(reference_back, motion_vectors, 16);
        frame_back = prediction_back + diff_back;
        if(writenowFlag)
            outputRGB(:,:, outputRGBIndex) = uint8(frame_back);
            outputRGBIndex = 1;

            writeVideo(outputVideo,uint8(outputRGB));

            outputRGB = zeros(L,W,3);
        else
            outputRGB(:,:, outputRGBIndex) = uint8(frame_back);
            outputRGBIndex = outputRGBIndex + 1;
        end
    end
    if(writenowFlag)
        counter = counter + 1;
    end
    rgbcounter = rgbcounter + 1;
end

if(outputName ~= "")
    close(outputVideo);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Other Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Image_back = JPEG_Decoder_No_Huff(huff_deco, quant_level, L, W)

array_rld = run_length_decoding(huff_deco); %Run Length Decoding

block_no = 1; %the number of block I'm working on right now
%Two loops to hold the indices of the 8x8 block we are recovering right now
for r = 1 : 8 : L
    for c = 1 : 8 : W
        array_1D = array_rld((block_no - 1) * 64 + 1 : block_no * 64); %capturing each 1x64 block after Run Length Decoding (rld)
        block2d_rev = serpentine_rev(array_1D); %reversing the effect of the zigzag (from 1-D to 2-D)
        block2d_deQuant = round(block2d_rev.*QuantizationTable(quant_level)); %dequantizing the 8x8 block
        block2d = round(Block8_IDCT(block2d_deQuant)); %applying inverse DCT
        Image_back(r : r + 7, c : c + 7) = block2d; %adding the retrieved block back into the image
        block_no = block_no + 1;
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = run_length_decoding(input)

N = length(input);

i = 1; %index of the input
idx = 1; %index of the output

%iterate over the input array
while i <= N
    if input(i) ~= 0 %if the value is not zero, copy it as it is to output
        output(idx) = input(i);
        idx = idx + 1;
    else %otherwise, insert all the amount of zeros
        i = i + 1;
        zero_repeated = input(i); %amount of zeros (beside the zero)

        %iterate and insert zeros
        for j = 1 : zero_repeated
            output(idx) = 0;
            idx = idx + 1;
        end
    end
    i = i + 1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function block = serpentine_rev(array)

%same as serpentine.m but only differs in assigning array values into block

N = length(array); %dimension of the block

block = zeros(sqrt(N)); %the output MxM block (M = sqrt(N))
loops = 2 * sqrt(N) - 1; %number of loops (diagonal moves) required to construct the block

x = 1; %index to hold the row
y = 1; %index to hold the column
idx = 1; %index of array
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
    for j = 1 : i
        if even == 1    %if loop is even, we increment the x and decrement y
            x = x + 1;
            y = y - 1;
        else            %if loop is even, we decrement the x and increment y
            x = x - 1;
            y = y + 1;
        end
        if (x > sqrt(N) || y > sqrt(N)) %we did not go out of bound
            continue;
        end
        block(x, y) = array(idx); %assign the corresponding pixel in its place in the 2-D block
        idx = idx + 1; %increment the index of the 1-D array
    end
end

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
function output = myHuffDec(input,dict)
if(istable(dict))
    dict = table2cell(dict);
end
codeWords = dict(:,2)';
output = [];
c = 1;
while(true)
    if(isempty(input))
        break;
    end
    initialWord = input(1:c);
    initialWordStr = string(sprintf('%d', initialWord));

    if(class(codeWords{1}) == 'string')
        index = find(cellfun(@(x) x == initialWordStr, codeWords, 'UniformOutput', 1),1,'first');
    else
        index = find(cellfun(@(x) string(sprintf('%d', x)) == initialWordStr, codeWords, 'UniformOutput', 1),1,'first');
    end

    if(index)
        output = [output dict(index,1)];
        input = [input(c+1:end)];
        c = 1;
    else
        c = c + 1;

    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%













