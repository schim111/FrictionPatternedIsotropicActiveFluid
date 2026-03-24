function MATS = MatAssemCAniso2D()

    %define domains
    Omega = Domain('triangle');
    Bdy = Domain('interval') < Omega;

    %Define FE Spaces
    Scalar_P1 = Element(Omega,lagrange_deg1_dim2,1);
    Vector_P1 = Element(Omega,lagrange_deg1_dim2,2);

    %test and trial functions
    cc = Trial(Scalar_P1);
    dd = Test(Scalar_P1);

    %Coef funcs
    c_k = Coef(Scalar_P1);
    v = Coef(Vector_P1);
    dt = Constant(Omega);
    D = Constant(Omega);

    %Geofunc on Sigma
    gSigma = GeoFunc(Bdy);

    %Define FE Matrices
    Mass_Matrix = Bilinear(Scalar_P1,Scalar_P1);
    I1 = Integral(Omega,cc.val*dd.val);
    Mass_Matrix = Mass_Matrix + I1;

    Stiff_Matrix = Bilinear(Scalar_P1,Scalar_P1);
    I2 = Integral(Omega,-dt.val * (cc.val * v.val' - D.val * cc.grad') * dd.grad);
    Stiff_Matrix = Stiff_Matrix + I2;

    Bdy_Matrix = Bilinear(Scalar_P1,Scalar_P1);
    I4 = Integral(Bdy,dt.val * gSigma.N' * (cc.val * v.val - D.val * cc.grad) * dd.val);
    Bdy_Matrix = Bdy_Matrix + I4;

    RHS = Linear(Scalar_P1);
    I3 = Integral(Omega, c_k.val * dd.val);
    RHS = RHS + I3;

    Quad_Order = 4;

    G1 = GeoElement(Omega);

    MATS = Matrices(Quad_Order,G1);

    %Collect Matrices
    MATS = MATS.Append_Matrix(Mass_Matrix);
    MATS = MATS.Append_Matrix(Stiff_Matrix);
    MATS = MATS.Append_Matrix(RHS);
    MATS = MATS.Append_Matrix(Bdy_Matrix);


end