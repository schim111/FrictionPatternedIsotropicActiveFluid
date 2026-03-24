%BJGSolver2D solves from an initial 2D c0 on a Mesh, the time evolution of
%the system according to the BJG active fluid, currently has f(c) = c/(1+c)
%hardcoded. Outputs c(t) and v(t) (scalar and vector field). Inputs are the
%Mesh, the Fixed Degrees of Freedom if Dirichlet conditions are used,
%c(0)=c0, Diffusion coefficient D, Peclet Number Pe, Shear screening length
%Lhs, Bulk screening length Lhb, Anisotropy strength Delta (scalar field),
%timestep dt, number of iterations NumIts, anisotropy components g,h (must
%have g^2 + h^2 = 1 or will not work right), PBC maps L2RMap, B2TMap,
%CornerVerts (setting these to [] will use free conditions on
%the boundary)

function [c_F,v_F] = BJGSolver2D(Mesh,FixedDoF,c0,D,Pe,Lhs,Lhb,Delta,dt,NumIts,g,h,L2RMap,B2TMap,CornerVerts)
    
    %Get size of mesh
    M = size(Mesh.Points,1);

    %Get DoFMap
    DoFMap = uint32(Mesh.ConnectivityList);

    %Initialize FEM Mats
    FEMV = [];
    FEMC = [];

    %Correct PBC Maps for Velocity DoF
    VL2RMap = [L2RMap;L2RMap + M];
    VB2TMap = [B2TMap;B2TMap + M];
    VCornerVerts = [CornerVerts;CornerVerts + M];

    c_F = zeros(M,NumIts + 1);
    v_F = zeros(M,2,NumIts);
    c_F(:,1) = c0;

    %Get subdomain embedding data for free/periodic BC
    DoI_Names = {'Omega','Bdy'};
    Subdomain_Embed = Mesh.Generate_Subdomain_Embedding_Data(DoI_Names);

    for its = 1:NumIts

        %Display progress to keep track
        if mod(its,10) == 0
            disp(num2str(its));
        end

        %Active stress
        f = c_F(:,its)./(1 + c_F(:,its));

        %Get Velocity solution matrix
        VFreeDoF = setdiff((1:1:2*M)',[FixedDoF;FixedDoF + M]);

        FEMV = VEquationAniso2D(FEMV,Mesh.Points,DoFMap,[],Subdomain_Embed,DoFMap,DoFMap,D,Delta,Lhb,Lhs,Pe,f,g,h);

        LHS = FEMV(3).MAT + FEMV(5).MAT;
        RHS = FEMV(4).MAT;
        
        %Apply PBC
        if ~isempty(VL2RMap) || ~isempty(VB2TMap)
            LHS = LHS';
            if ~isempty(VL2RMap)
                LHS(:,VL2RMap(:,1)) = LHS(:,VL2RMap(:,1)) + LHS(:,VL2RMap(:,2));
                RHS(VL2RMap(:,1)) = RHS(VL2RMap(:,1)) + RHS(VL2RMap(:,2));
                RHS(VL2RMap(:,2)) = 0;
                LHS(:,VL2RMap(:,2)) = 0;
            end
            if ~isempty(VB2TMap)
                LHS(:,VB2TMap(:,1)) = LHS(:,VB2TMap(:,1)) + LHS(:,VB2TMap(:,2));
                RHS(VB2TMap(:,1)) = RHS(VB2TMap(:,1)) + RHS(VB2TMap(:,2));
                RHS(VB2TMap(:,2)) = 0;
                LHS(:,VB2TMap(:,2)) = 0;
            end
            if ~isempty(VCornerVerts)
                LHS(:,VCornerVerts(1)) = LHS(:,VCornerVerts(1)) + LHS(:,VCornerVerts(2)) + LHS(:,VCornerVerts(3)) + LHS(:,VCornerVerts(4));
                RHS(VCornerVerts(1)) = RHS(VCornerVerts(1)) + RHS(VCornerVerts(2)) + RHS(VCornerVerts(3)) + RHS(VCornerVerts(4));
                RHS(VCornerVerts(2:end)) = 0;
                LHS(:,VCornerVerts(2:end)) = 0;
            end
            for row = 1:length(VL2RMap)
                i = VL2RMap(row,1);
                j = VL2RMap(row,2);
                LHS(j,j) = 1;
                LHS(i,j) = -1;
            end
            for row = 1:length(VB2TMap)
                i = VB2TMap(row,1);
                j = VB2TMap(row,2);
                LHS(j,j) = 1;
                LHS(i,j) = -1;
            end
            for row = 2:length(VCornerVerts)
                i = VCornerVerts(1);
                j = VCornerVerts(row);
                LHS(j,j) = 1;
                LHS(i,j) = -1;
            end
            LHS = LHS';
        end

        %Solve for velocity
        vv = zeros(2*M,1);
        vv(VFreeDoF,1) = LHS(VFreeDoF,VFreeDoF) \ RHS(VFreeDoF);
        v_F(:,:,its) = [vv(1:M,1),vv(M+1:end,1)];

        %Get updated density matrix
        FreeDoF = setdiff((1:1:M)',FixedDoF);

        FEMC = CEquationAniso2D(FEMC,Mesh.Points,DoFMap,[],Subdomain_Embed,DoFMap,DoFMap,D,c_F(:,its),dt,v_F(:,:,its));

        LHS = FEMC(2).MAT + FEMC(4).MAT;
        RHS = FEMC(3).MAT;


        %Apply PBC
        if ~isempty(L2RMap) || ~isempty(B2TMap)
            LHS = LHS';
            if ~isempty(L2RMap)
                LHS(:,L2RMap(:,1)) = LHS(:,L2RMap(:,1)) + LHS(:,L2RMap(:,2));
                RHS(L2RMap(:,1)) = RHS(L2RMap(:,1)) + RHS(L2RMap(:,2));
                RHS(L2RMap(:,2)) = 0;
                LHS(:,L2RMap(:,2)) = 0;
            end
            if ~isempty(B2TMap)
                LHS(:,B2TMap(:,1)) = LHS(:,B2TMap(:,1)) + LHS(:,B2TMap(:,2));
                RHS(B2TMap(:,1)) = RHS(B2TMap(:,1)) + RHS(B2TMap(:,2));
                RHS(B2TMap(:,2)) = 0;
                LHS(:,B2TMap(:,2)) = 0;
            end
            if ~isempty(CornerVerts)
                LHS(:,CornerVerts(1)) = LHS(:,CornerVerts(1)) + LHS(:,CornerVerts(2)) + LHS(:,CornerVerts(3)) + LHS(:,CornerVerts(4));
                RHS(CornerVerts(1)) = RHS(CornerVerts(1)) + RHS(CornerVerts(2)) + RHS(CornerVerts(3)) + RHS(CornerVerts(4));
                RHS(CornerVerts(2:end)) = 0;
                LHS(:,CornerVerts(2:end)) = 0;
            end
            for row = 1:length(L2RMap)
                i = L2RMap(row,1);
                j = L2RMap(row,2);
                LHS(j,j) = 1;
                LHS(i,j) = -1;
            end
            for row = 1:length(B2TMap)
                i = B2TMap(row,1);
                j = B2TMap(row,2);
                LHS(j,j) = 1;
                LHS(i,j) = -1;
            end
            for row = 2:length(CornerVerts)
                i = CornerVerts(1);
                j = CornerVerts(row);
                LHS(j,j) = 1;
                LHS(i,j) = -1;
            end
            LHS = LHS';
        end

        %Solve for new density
        c_F(FreeDoF,its+1) = LHS(FreeDoF,FreeDoF) \ RHS(FreeDoF,1);

        %Check for convergence
        if max(abs(c_F(:,its+1) - c_F(:,its))) < 1e-5
            disp(['Converged at timestep ',num2str(its)])
            c_F(:,its+2:end) = [];
            v_F(:,:,its+1:end) =[];
            break;
        end
    end

end