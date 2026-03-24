function MATS = MatAssemVAniso2D()

    %Domain
    Omega = Domain('triangle');
    Bdy = Domain('interval') < Omega;

    %Define FE Spaces
    Scalar_P1 = Element(Omega,lagrange_deg1_dim2,1);
    Vector_P1 = Element(Omega,lagrange_deg1_dim2,2);

    %test and trial functions
    vv = Trial(Vector_P1);
    ww = Test(Vector_P1);

    %Coef funcs
    f = Coef(Scalar_P1);
    g = Coef(Scalar_P1);
    h = Coef(Scalar_P1);
    Delta = Coef(Scalar_P1);
    Lhs = Constant(Omega);
    Lhb = Constant(Omega);
    D = Constant(Omega);
    Pe = Constant(Omega);
    
    %Geofunc on Sigma
    gSigma = GeoFunc(Bdy);

    %Define FE Matrices
    Mass_Matrix = Bilinear(Vector_P1,Vector_P1);
    I1 = Integral(Omega,vv.val' * ww.val);
    I2 = Integral(Omega,Delta.val * (h.val*h.val*vv.val(1)*ww.val(1) + g.val*g.val*vv.val(2)*ww.val(2) - g.val*h.val*(vv.val(1)*ww.val(2) + vv.val(2)*ww.val(1))));
    Mass_Matrix = Mass_Matrix + I1 + I2;

    Stiff_Matrix = Bilinear(Vector_P1,Vector_P1);
    I3 = Integral(Omega,Lhs.val * 2 * (vv.grad(1,1)*ww.grad(1,1) + vv.grad(2,2)*ww.grad(2,2)));
    I4 = Integral(Omega,Lhs.val * (vv.grad(1,2) + vv.grad(2,1)) * (ww.grad(1,2) + ww.grad(2,1)));
    I5 = Integral(Omega,Lhb.val * (vv.grad(1,1) + vv.grad(2,2)) * (ww.grad(1,1) + ww.grad(2,2)));
    Stiff_Matrix = Stiff_Matrix + I3 + I4 + I5;

    Bdy_Matrix = Bilinear(Vector_P1,Vector_P1);
    I7 = Integral(Bdy,-Lhs.val * (2 * gSigma.N(1) * vv.grad(1,1) * ww.val(1) + 2 * gSigma.N(2) * vv.grad(2,2) * ww.val(2) + (vv.grad(1,2) + vv.grad(2,1)) * (gSigma.N(1) * ww.val(2) + gSigma.N(2) * ww.val(1))));
    I8 = Integral(Bdy,-Lhb.val * (vv.grad(1,1) + vv.grad(2,2)) * gSigma.N' * ww.val);
    Bdy_Matrix = Bdy_Matrix + I7 + I8;

    RHS = Linear(Vector_P1);
    I6 = Integral(Omega,-D.val * Pe.val * f.val * (ww.grad(1,1) + ww.grad(2,2)));
    RHS = RHS + I6;

    Bdy_RHS = Linear(Vector_P1);
    I9 = Integral(Bdy,D.val * Pe.val * f.val * gSigma.N' * ww.val);
    Bdy_RHS = Bdy_RHS + I9;

    Quad_Order = 4;
    
    G1 = GeoElement(Omega);

    MATS = Matrices(Quad_Order,G1);

    %Collect Matrices
    MATS = MATS.Append_Matrix(Mass_Matrix);
    MATS = MATS.Append_Matrix(Stiff_Matrix);
    MATS = MATS.Append_Matrix(RHS);
    MATS = MATS.Append_Matrix(Bdy_Matrix);
    MATS = MATS.Append_Matrix(Bdy_RHS);

end