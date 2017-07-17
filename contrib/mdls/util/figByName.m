function varargout = figByName(string)
% function figByName(string)
%
% My kind of figures, it add the possibility of call a figure with a string
% instead of a number.
%
% Author: federico lecumberry <fefo at fing dot edu dot uy>
%
% Version: $Id: Figure.m 202 2008-05-16 22:22:40Z fefo $
%
if (nargout >= 1)
  returnHandle = true;
else
  returnHandle = false;
end
if (nargin == 0)
  if returnHandle
    hf = figure;
  else
    figure
  end
elseif isnumeric(string)
  if returnHandle
    hf = figure(string);
  else
    figure(string)
  end
else
  nmbr = hashstr(string);
  hf = figure(nmbr);
  set(hf,'NumberTitle','off')
  titleText = sprintf('%s',string);
  set(hf,'Name', titleText)
end

if returnHandle
  varargout{1} = hf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hh = hashstr(X)
% y = MessageDigest(X) Computes the 160-bit Hex values of the
% secure hash algorithm for the ASCII string X into
y = MessagePad(X);
M = size(y,1);
H0 = hex2dec('67452301');
H1 = hex2dec('EFCDAB89');
H2 = hex2dec('98BADCFE');
H3 = hex2dec('10325476');
H4 = hex2dec('C3D2E1F0');
for i = 1:M
  W(1:16) = y(i,1:16);
  for t = 17:80
    W(t) = S(bitxor(bitxor(bitxor(W(t-3),W(t-8)),W(t-14)),W(t-16)),1);
  end
  A = H0;B = H1;C = H2;D = H3;E = H4;
  for t = 1:80
    TEMP = safe_add(safe_add(safe_add(safe_add(S(A,5),func(B,C,D,t-1)),E),W(t)),kay(t-1));
    E = D;D = C;C = S(B,30);B = A;A = TEMP;
  end
  H0 = safe_add(H0,A);
  H1 = safe_add(H1,B);
  H2 = safe_add(H2,C);
  H3 = safe_add(H3,D);
  H4 = safe_add(H4,E);
end
%Z = [dec2hex(H0,8);dec2hex(H1,8);dec2hex(H2,8);dec2hex(H3,8);dec2hex(H4,8)];
Z = H0;
hh = rem(Z,999999999);

function z = safe_add(x,y)
% internal function of secure hash algorithm
z = x+y; if z > 2^32, z = z-2^32; end

function y = MessagePad(X)
% y = MessagePad(X) Converts the ASCII string X into a matrix y
% Each row contains 16 32-bit non negative integers and each row
% represts 512 bit of the ASCII message. The message is then padded
% with 1 followed by m 0's followed by L which is a 64-bit integer
% representing the length of the message
M = length(X);L = M*8;m = mod(512-mod(L+65,512),512);%m is number of zero bits
y = [];W = dec2bin(L,64);
i = 1;j = 1;tmp = M;
while tmp>= 0
  if tmp>= 4
    y(i,j) = double(X(i+4*j-4))*16777216+double(X(i+4*j-3))*65536+double(X(i+4*j-2))*256+double(X(i+4*j-1));
  else
    switch tmp
      case 3
        y(i,j) = double(X(i+4*j-4))*16777216+double(X(i+4*j-3))*65536+double(X(i+4*j-2))*256+128;
        m = m-7;
      case 2
        y(i,j) = double(X(i+4*j-4))*16777216+double(X(i+4*j-3))*65536+32768;
        m = m-15;
      case 1
        y(i,j) = double(X(i+4*j-4))*16777216+8388608;
        m = m-23;
      case 0
        y(i,j) = 2^31;
        m = m-31;
    end
  end
  j = j+1;
  if j>16
    i = i+1;j = 1;
  end
  tmp = tmp-4;
end
tmp = m;
while tmp>= 0
  if tmp == 0
    y(i,j) = bin2dec(W(1:32));
    j = j+1;
    if j>16
      i = i+1;j = 1;
    end
    y(i,j) = bin2dec(W(33:64));
  end
  j = j+1;
  if j>16
    i = i+1;j = 1;
  end
  tmp = tmp-32;
end

function y = S(X,n)
% y = S(X,n) Circular left bit shift operator of 32-bit nonnegative integer X
% by n bits If n is negative the shift will be a right shift
m = mod(n,32);
y = bitor(bitshift(X,m,32),bitshift(X,m-32,32));

function f = func(B,C,D,t)
% internal function of secure hash algorithm
f = 0;
if (t >= 0) && (t <= 19)
    f = bitor(bitand(B,C),bitand(bitcmp(B,32),D));
elseif (t >= 20) && (t <= 39)
    f = bitxor(B,bitxor(C,D));
elseif (t >= 40) && (t <= 59)
    f = bitor(bitor(bitand(B,C),bitand(B,D)),bitand(C,D));
elseif (t >= 60) && (t <= 79)
    f = bitxor(B,bitxor(C,D));
end

function k = kay(t)
% internal function of secure hash algorithm
k = 0;
if (t >= 0) && (t <= 19)
    k = hex2dec('5A827999');
elseif (t >= 20) && (t <= 39)
    k = hex2dec('6ED9EBA1');
elseif (t >= 40) && (t <= 59)
    k = hex2dec('8F1BBCDC');
elseif (t >= 60) && (t <= 79)
    k = hex2dec('CA62C1D6');
end
