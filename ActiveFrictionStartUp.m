%Setup file for isotropic active matter model with anisotropic friction patterning

%Note on bulk screening length/viscosity: The solver explicitly solves
%Stokes equation with a viscous stress Lhs^2 ( grad(v) + grad(v)^T ) +
%Lhb^2 div(v) I. In order to match the viscous stress used in the text
%(which is more common in the literature) we set, for most simulations
%Lhb = Lhs*(1i) so that Lhb^2 = -Lhs^2 so that \ell_b^2 as defined in the
%text is \ell_b^2 = 0. This can be changed as needed for different values
%of \ell_b.

sd = 4148; %rng seed 
rng(sd);

%Model Parameters
D = 1; %Diffusion, sets time scale
Pe = 50; %Peclet Number
Lhs = 0.075; %Shear screening length
Lhb = Lhs*(1i); %Bulk screening length: See above Note
delt = 3; %Friction anisotropy
dt = 0.01; %Time Step
NumIts = 500; %Total number of iterations
Theta0 = pi/2;%linspace(0,pi/2,101);

%Create Domain Omega
L = 2; %System size (note from zero will be L/2)
[TRI,VTX] = bcc_triangle_mesh(150,150);
VTX(:,1) = 0.5*L*(2*VTX(:,1) - 1);
VTX(:,2) = 0.5*L*(2*VTX(:,2) - 1);
Mesh = MeshTriangle(TRI,VTX,'Omega');
clear TRI VTX

%Get Boundaries of Mesh
Bdy_Edges = Mesh.freeBoundary();
Mesh = Mesh.Append_Subdomain('1D','Bdy',Bdy_Edges);
Bdy_Verts = unique(Bdy_Edges(:));

%Get Left,Right,Top,Bottom Edges (For PBC)
Bdy_Edge_Centers = 0.5*(Mesh.Points(Bdy_Edges(:,1),:) + Mesh.Points(Bdy_Edges(:,2),:));
Left_Mask = (Bdy_Edge_Centers(:,1) < -0.5*L + 1e-5);
Right_Mask = (Bdy_Edge_Centers(:,1) > 0.5*L - 1e-5);
Top_Mask = (Bdy_Edge_Centers(:,2) > 0.5*L - 1e-5);
Bottom_Mask = (Bdy_Edge_Centers(:,2) < -0.5*L + 1e-5);
Left_Edges = Bdy_Edges(Left_Mask,:);
Right_Edges = Bdy_Edges(Right_Mask,:);
Top_Edges = Bdy_Edges(Top_Mask,:);
Bottom_Edges = Bdy_Edges(Bottom_Mask,:);
Left_Verts = unique(Left_Edges(:));
Right_Verts = unique(Right_Edges(:));
Top_Verts = unique(Top_Edges(:));
Bot_Verts = unique(Bottom_Edges(:));
%Get Left->Right map by sorting y components
yL = Mesh.Points(Left_Verts,2);
yR = Mesh.Points(Right_Verts,2);
%Get Bot->Top map by sorting x components
xB = Mesh.Points(Bot_Verts,1);
xT = Mesh.Points(Top_Verts,1);
if (length(yL) ~= length(yR) | length(xB) ~= length(xT))
    error('Edges not commensurate for PBC!');
end
[~,IL] = sort(yL);
[~,IR] = sort(yR);
[~,IB] = sort(xB);
[~,IT] = sort(xT);
Left_Verts = Left_Verts(IL);
Right_Verts = Right_Verts(IR);
Bot_Verts = Bot_Verts(IB);
Top_Verts = Top_Verts(IT);
CornerVerts = [Bot_Verts(1);Bot_Verts(end);Top_Verts(1);Top_Verts(end)];
%CornerVerts = [];
Left_Verts = Left_Verts(2:end-1);
Right_Verts = Right_Verts(2:end-1);
Top_Verts = Top_Verts(2:end-1);
Bot_Verts = Bot_Verts(2:end-1);
L2RMap = [Left_Verts,Right_Verts];
%L2RMap = [];
B2TMap = [Bot_Verts,Top_Verts];
%FixedDoF = [Left_Verts,Right_Verts];
FixedDoF = [];


%Set up anisotropic friction (Make sure Delta*g and Delta*h satisfy PBCs)
%+1 Defect

r = sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2);
g = (1./(r + eps)).*(Mesh.Points(:,1)*cos(Theta0(kk)) - Mesh.Points(:,2)*sin(Theta0(kk)));
h = (1./(r + eps)).*(Mesh.Points(:,1)*sin(Theta0(kk)) + Mesh.Points(:,2)*cos(Theta0(kk)));
Delta = 0*Mesh.Points(:,1);
DMask = r < 0.95;
Delta(DMask) = delt(gg);

%Horizontal Lines
%{
g = ones(size(Mesh.Points(:,1),1),1);
h = zeros(size(Mesh.Points(:,1),1),1);
Delta = delt*ones(size(Mesh.Points(:,1),1),1);
%}

%Chevron Pattern
%g = (cos(pi/4))*ones(size(Mesh.Points(:,1),1),1);
%h = (sin(pi/4))*ones(size(Mesh.Points(:,1),1),1);
%mxMask = Mesh.Points(:,1) < 0;
%h(mxMask) = -h(mxMask);
%Delta = delt*ones(size(Mesh.Points(:,1),1),1);

%ZigZag Chevron Pattern
%{
q1mask = (Mesh.Points(:,1) >= -Mesh.Points(:,2) + 0.5) & (Mesh.Points(:,2) >= 0);
q2mask = (Mesh.Points(:,1) < -Mesh.Points(:,2) + 0.5) & (Mesh.Points(:,2) >= 0);
q3mask = (Mesh.Points(:,1) >= Mesh.Points(:,2) + 0.5) & (Mesh.Points(:,2) < 0);
q4mask = (Mesh.Points(:,1) < Mesh.Points(:,2) + 0.5) & (Mesh.Points(:,2) < 0);
%g(mxMask) = -g(mxMask);
g(q1mask) = 0;
h(q1mask) = 1;
g(q2mask) = 1;
h(q2mask) = 0;
g(q3mask) = 1;
h(q3mask) = 0;
g(q4mask) = 0;
h(q4mask) = 1;
dmask = abs(Mesh.Points(:,1) - 0.5 + abs(Mesh.Points(:,2))) > 0.2;
%Delta = delt*ones(size(Mesh.Points(:,1),1),1);
Delta(dmask) = 0;

%}

%Set initial conditions on c (note that the total c will be conserved)
CC = 0.25;
r = sqrt(Mesh.Points(:,1).^2 + Mesh.Points(:,2).^2);
c0 = exp(-((Mesh.Points(:,1)).^2 + Mesh.Points(:,2).^2)./0.1) + CC;

%Offcenter Gaussians
%c0 = exp(-((Mesh.Points(:,1) - 0.25).^2 + Mesh.Points(:,2).^2)./0.01) + exp(-((Mesh.Points(:,1) - 1.25).^2 + Mesh.Points(:,2).^2)./0.01) + exp(-((Mesh.Points(:,1) + 0.75).^2 + Mesh.Points(:,2).^2)./0.01) + exp(-((Mesh.Points(:,1) - 0.25).^2 + (Mesh.Points(:,2) - 1).^2)./0.01) + exp(-((Mesh.Points(:,1) - 0.25).^2 + (Mesh.Points(:,2) + 1).^2)./0.01) + exp(-((Mesh.Points(:,1) - 1.25).^2 + (Mesh.Points(:,2) - 1).^2)./0.01) + exp(-((Mesh.Points(:,1) - 1.25).^2 + (Mesh.Points(:,2) + 1).^2)./0.01) + exp(-((Mesh.Points(:,1) + 0.75).^2 + (Mesh.Points(:,2) - 1).^2)./0.01) + exp(-((Mesh.Points(:,1) + 0.75).^2 + (Mesh.Points(:,2) + 1).^2)./0.01) + CC;

%Respect PBC
c0(Left_Verts) = c0(Right_Verts);
%c0(Left_Verts) = 0;
%c0(Right_Verts) = 0;
c0(Top_Verts) = c0(Bot_Verts);
c0(CornerVerts(2:4)) = c0(CornerVerts(1));


[c,v] = BJGSolver2D(Mesh,FixedDoF,c0,D,Pe(ii),Lhs(jj),Lhb,Delta,dt,NumIts,g,h,L2RMap,B2TMap,CornerVerts);

%Save the results
dir = ''; %Directory to save results
filename = [dir,'','.mat']; %Filename of results
save(filename,'Mesh','c','v'); %Saves Mesh, density field, and velocity field data

clear
