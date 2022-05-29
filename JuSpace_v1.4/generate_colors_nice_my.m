function colors = generate_colors_nice_my(N)
% function colors = generate_colors_nice_my(N)
% from pa_LCH2RGB
N = N+1;
L = repmat(65,N,1);
C = repmat(75,N,1);
H = linspace(25,385,N)';


a = cos(H*pi/180).*C;
b = sin(H*pi/180).*C;

%% Different L*ab -> XYZ from http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html

f_y=(L+16)/116;
f_x	= f_y+a/500;
f_z	= f_y-b/200;

epsilon=0.008856;
kappa=903.3;

X=f_x.^3;
sel=f_x.^3<=epsilon;
X(sel)=(116*f_x(sel)-16)/kappa;

Y=((L+16)/116).^3;
sel=L<=(kappa*epsilon);
Y(sel)=L(sel)/kappa;

Z=f_z.^3;
sel=f_z<=epsilon;
Z(sel)=(116*f_z(sel)-16)/kappa;

XYZ = [X Y Z];

ref	=  [95.05, 100.000, 108.90]/100;
XYZ =  bsxfun(@times,XYZ,ref);


    T=[3.2406, -1.5372, -0.4986;...
        -0.9689, 1.8758, 0.0415;...
        0.0557, -0.2040, 1.0570];
    
    RGB		= XYZ*T';
    
    % Gamma correction to convert RGB to sRGB
    sel			= RGB>0.0031308; % colorspace: 0.0031306684425005883
    RGB(sel)	= 1.055*(RGB(sel).^(1/2.4))-0.055;
    RGB(~sel)	= 12.92*RGB(~sel);
    
   
RGB(RGB>1)	= 1;
RGB(RGB<0)	= 0;
colors = RGB(1:end-1,:);
