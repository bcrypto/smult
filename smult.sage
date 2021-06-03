#field ops
def MontInv(u):
    keys = tuple(u.keys())    
    v = [None] * len(keys)
    v[0] = u[keys[0]]
    for i in range(1, len(keys)):
        v[i]  = v[i-1]*u[keys[i]]
      
    t = 1 / v[len(keys) - 1]
    ui = {}
    for i in range(len(keys) - 1, 0, -1):
        ui[keys[i]] = v[i-1] * t
        t = t * u[keys[i]]
    ui[keys[0]] = t
    return ui

def mul2(a, b, a2, b2):
    return ((a + b)^2 - a2 - b2) / 2

#ec points ops
def ecX(P):
    return P[0]

def ecY(P):
    return P[1]

def ecZ(P):
    return P[2]

def ecJToA(P):
    ZI = 1 / (ecZ(P))
    ZIZI = ZI^2
    return (ecX(P)*ZIZI, ecY(P)*ZIZI*ZI, 1)

def ecHToA(P):
    ZI = 1 / (ecZ(P))
    return (ecX(P)*ZI, ecY(P)*ZI, 1)

def ecJToH(P):
    return (ecX(P)*ecZ(P), ecY(P), ecZ(P)^3)

def ecEqPointAA(P1, P2):
    return tuple(P1)==tuple(P2)

def ecEqPointAJ(P1, P2):
    return ecEqPointAA(P1, ecJToA(P2))
    
def allEq(points1, points2, cmp):
    if (len(points1) != len(points2)):
        return False
    for n in points1.keys():
        if not cmp(points1[n], points2[n]):
            return False
    return True

#https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#addition-add-2007-bl
def ecAddJJ(E, P1, P2, check_same=False):
    X1, Y1, Z1 = P1
    X2, Y2, Z2 = P2
    Z1Z1 = Z1^2
    Z2Z2 = Z2^2
    U1 = X1*Z2Z2
    U2 = X2*Z1Z1
    S1 = Y1*Z2*Z2Z2
    S2 = Y2*Z1*Z1Z1
    if check_same and U1 == U2 and S1 == S2:
        print(P1, P2)
        raise ValueError('Adding equivalent points')
    H = U2-U1
    I = (2*H)^2
    J = H*I
    r = 2*(S2-S1)
    V = U1*I
    X3 = r^2-J-2*V
    Y3 = r*(V-X3)-2*S1*J
    Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H
    return (X3, Y3, Z3)

#https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl
def ecDblJ(E, P):
    a = E.a4()
    X1, Y1, Z1 = P
    XX = X1^2
    YY = Y1^2
    YYYY = YY^2
    ZZ = Z1^2
    S = 2*((X1+YY)^2-XX-YYYY)
    M = 3*XX+a*ZZ^2
    T = M^2-2*S
    X3 = T
    Y3 = M*(S-T)-8*YYYY
    Z3 = (Y1+Z1)^2-YY-ZZ
    return (X3, Y3, Z3)

def ecDblJAddA(E, P1J, P2):
    t1, t2, t3 = P1J
    tx, ty, _ = P2
    t4 = t3^2
    t5 = tx * t4
    t5 = t5 - t1
    t6 = t3 + t5
    t6 = t6^2
    t6 = t6 - t4
    t4 = t3 * t4
    t4 = ty * t4
    t4 = t4 - t2
    t3 = t5^2
    t6 = t6 - t3
    t1 = t1 * t3
    t1 = 4 * t1
    t3 = t3 * t5
    t2 = t2 * t3
    t2 = 8*t2
    t5 = t4^2
    t3 = t5 - t3
    t3 = 4 * t3
    t3 = t3 - t1
    t3 = t3 - t1
    t3 = t3 - t1
    t4 = t3 + t4
    t4 = t4^2
    t4 = t5 - t4
    t4 = t4 - t2
    t4 = t4 - t2
    t5 = t3^2
    t4 = t4 + t5
    t1 = t1 * t5
    t5 = t3 * t5
    t3 = t3 * t6
    t2 = t2 * t5
    #36
    t6 = t4^2
    t6 = t6 - t5
    #36
    t5 = 3*t1
    t5 = t5 - t6
    t4 = t4*t5
    t2 = t4 - t2
    t1 = t1 - t5
    return (t1, t2, t3)

def ecDblAdd(E, P1J, P2, smult_alg):
    if smult_alg == SmallMultA:
        return ecDblJAddA(E, P1J, P2)
    if smult_alg == SmallMultJ:
        return ecAddJJ(E, ecDblJ(E, P1J), P2, check_same=True)
    return None

#https://eprint.iacr.org/2015/1060.pdf
#Algorithm 1
def ecAddJJ_complete(E, P1, P2):
    X1, Y1, Z1 = ecJToH(P1)
    X2, Y2, Z2 = ecJToH(P2)
    a = E.a4()
    b3 = 3 * E.a6()
    t0 = X1 * X2 
    t1 = Y1 * Y2 
    t2 = Z1 * Z2
    t3 = X1 + Y1 
    t4 = X2 + Y2 
    t3 = t3 * t4
    t4 = t0 + t1 
    t3 = t3 - t4 
    t4 = X1 + Z1
    t5 = X2 + Z2 
    t4 = t4 * t5 
    t5 = t0 + t2
    t4 = t4 - t5 
    t5 = Y1 + Z1 
    X3 = Y2 + Z2
    t5 = t5 * X3 
    X3 = t1 + t2 
    t5 = t5 - X3
    Z3 = a * t4 
    X3 = b3 * t2 
    Z3 = X3 + Z3
    X3 = t1 - Z3 
    Z3 = t1 + Z3 
    Y3 = X3 * Z3
    t1 = t0 + t0 
    t1 = t1 + t0 
    t2 = a * t2
    t4 = b3 * t4 
    t1 = t1 + t2 
    t2 = t0 - t2
    t2 = a * t2 
    t4 = t4 + t2 
    t0 = t1 * t4
    Y3 = Y3 + t0 
    t0 = t5 * t4 
    X3 = t3 * X3
    X3 = X3 - t0 
    t0 = t3 * t1 
    Z3 = t5 * Z3
    Z3 = Z3 + t0
    return (X3, Y3, Z3)

#https://eprint.iacr.org/2015/1060.pdf
#Algorithm 2
def ecAddJA_complete(E, P1, P2):
    X1, Y1, Z1 = ecJToH(P1)
    X2, Y2, Z2 = P2
    a = E.a4()
    b3 = 3 * E.a6()
    t0 = X1 * X2 
    t1 = Y1 * Y2 
    t3 = X2 + Y2
    t4 = X1 + Y1 
    t3 = t3 * t4 
    t4 = t0 + t1
    t3 = t3 - t4 
    t4 = X2 * Z1 
    t4 = t4 + X1
    t5 = Y2 * Z1 
    t5 = t5 + Y1 
    Z3 = a * t4
    X3 = b3 * Z1 
    Z3 = X3 + Z3 
    X3 = t1 - Z3
    Z3 = t1 + Z3 
    Y3 = X3 * Z3 
    t1 = t0 + t0
    t1 = t1 + t0 
    t2 = a * Z1 
    t4 = b3 * t4
    t1 = t1 + t2 
    t2 = t0 - t2 
    t2 = a * t2
    t4 = t4 + t2 
    t0 = t1 * t4 
    Y3 = Y3 + t0
    t0 = t5 * t4 
    X3 = t3 * X3 
    X3 = X3 - t0
    t0 = t3 * t1
    Z3 = t5 * Z3
    Z3 = Z3 + t0
    return (X3, Y3, Z3)

def ecAdd_complete(E, P1, P2, smult_alg):
    if smult_alg == SmallMultA:
        return ecAddJA_complete(E, P1, P2)
    if smult_alg == SmallMultJ:
        return  ecAddJJ_complete(E, P1, P2)
    return None

def ecNeg(P):
    X, Y, Z = P
    return (X, -Y, Z)

#Small mults
def SmallMultA(E, P, w, check_divp=False):
    a = E.a4()
    b = E.a6()
    F = E.base_field()
    y2 = ecY(P)^2
    dy2 = 4 * y2
    a2 = a^2
    x2 = ecX(P)^2
    bx = b*ecX(P)
    t = (x2 + a)^2
    W = {1: F(1), 2: F(1), 3: 3*t - 4*(a2-3*bx)}
    ax = a * ecX(P)
    x3 = y2 - ax - b
    x6 = x3^2
    W[4] = 2 * (x6 + 4*bx*(5*x2-a)+5*ax*(x3-ax)-(a2*a+8*b*b))
    W2 = {2: F(1), 3: W[3]^2, 4: W[4]^2}
    WW = {1: W[1]*W[3], 2: W[2]*W[4]}
    WWy2 = {2: dy2 * WW[2]}
    WWy4 = {2: dy2 * WWy2[2]}
    W[5] = WWy4[2] - WW[1]*W2[3]
    W2[5] = W[5]^2
    for n in range(3, 2^(w-1) + 1):
        WW[n] = mul2(W[n], W[n+2], W2[n], W2[n+2])
        W[2*n] = WW[n]*W2[n-1] - WW[n-2]*W2[n+1]
        W2[2*n] = W[2*n]^2
        if is_odd(n):
            W[2*n + 1] = WW[n]*W2[n] - WWy4[n-1]*W2[n+1]
        else:
            WWy2[n] = dy2 * WW[n]
            WWy4[n] = dy2 * WWy2[n]
            W[2*n + 1] =  WWy4[n]*W2[n] - WW[n-1]*W2[n+1]
        if n != 2^(w-1):
            W2[2*n + 1] = W[2*n + 1]^2
    
    if (check_divp):
        expectedPoly = dict((k, E.division_polynomial(k, ecX(P), two_torsion_multiplicity=0)) for k in range(1,2^w + 2))
        if W != expectedPoly:
            print('calculated polynomials :', W)
            print('expected polynomials   :', expectedPoly)
            raise ValueError('Invalid division polynomials')

    W2i = MontInv(dict((k, W2[k]) for k in W2.keys() if is_odd(k)))
    X = {}
    for n in range(3, 2^(w-1) + 1 + 1, 2):
        X[n] = ecX(P) - WWy2[n-1]*W2i[n]
    for n in range(2^(w-1) + 3, 2^w - 1 + 1, 2):
        t = mul2(W[n-1], W[n+1], W2[n-1], W2[n+1])
        X[n] = ecX(P) - dy2*t*W2i[n]

    Y = {}
    for n in range(3, 2^(w-1) - 1 + 1, 2):
        Y[n] = ecY(P) * W[2*n] * W2i[n]^2
    for n in range(2^(w-1) + 1, 2^w - 3 + 1, 2):
        WW[n] = mul2(W[n], W[n+2], W2[n], W2[n+2])
        Y[n] = ecY(P) * (WW[n]*W2[n-1] - WW[n-2]*W2[n+1]) * W2i[n]^2
    Y[2^w - 1] = ecY(P) * (W[2^w - 1]*W[2^w + 1]*W2[2^w - 2] - WW[2^w - 3]*W2[2^w]) * W2i[2^w - 1]^2

    return dict((n, (X[n], Y[n], 1)) for n in range(3, 2^w, 2))

def SmallMultJ(E, P, w, check_divp=False):
    a = E.a4()
    b = E.a6()
    F = E.base_field()
    y2 = ecY(P)^2
    dy2 = 4 * y2
    a2 = a^2
    x2 = ecX(P)^2
    bx = b*ecX(P)
    t = (x2 + a)^2
    W = {1: F(1), 2: F(1), 3: 3*t - 4*(a2-3*bx)}
    ax = a * ecX(P)
    x3 = y2 - ax - b
    x6 = x3^2
    W[4] = 2 * (x6 + 4*bx*(5*x2-a)+5*ax*(x3-ax)-(a2*a+8*b*b))
    W2 = {2: F(1), 3: W[3]^2, 4: W[4]^2}
    WW = {1: W[1]*W[3], 2: W[2]*W[4]}
    WWy2 = {2: dy2 * WW[2]}
    WWy4 = {2: dy2 * WWy2[2]}
    W[5] = WWy4[2] - WW[1]*W2[3]
    W2[5] = W[5]^2
    WWW = {}
    for n in range(3, 2^(w-1) + 1):
        WWW[n] = W[n+2]*W2[n-1] - W[n-2]*W2[n+1]
        W[2*n] = W[n]*WWW[n]
        W2[2*n] = W[2*n]^2
        WW[n] = mul2(W[n], W[n+2], W2[n], W2[n+2])
        if is_odd(n):
            W[2*n + 1] = WW[n]*W2[n] - WWy4[n-1]*W2[n+1]
        else:
            WWy2[n] = dy2 * WW[n]
            WWy4[n] = dy2 * WWy2[n]
            W[2*n + 1] =  WWy4[n]*W2[n] - WW[n-1]*W2[n+1]
        if n != 2^(w-1):
            W2[2*n + 1] = W[2*n + 1]^2
    
    if (check_divp):
        expectedPoly = dict((k, E.division_polynomial(k, ecX(P), two_torsion_multiplicity=0)) for k in range(1,2^w + 2))
        if W != expectedPoly:
            print('calculated polynomials :', W)
            print('expected polynomials   :', expectedPoly)
            raise ValueError('Invalid division polynomials')

    X = {}
    for n in range(3, 2^(w-1) + 1 + 1, 2):
        X[n] = ecX(P)*W2[n] - WWy2[n-1]
    for n in range(2^(w-1) + 3, 2^w - 1 + 1, 2):
        t = mul2(W[n-1], W[n+1], W2[n-1], W2[n+1])
        X[n] = ecX(P)*W2[n] - dy2*t
    Y = {}
    for n in range(3, 2^(w-1) - 1 + 1, 2):
        Y[n] = ecY(P) * WWW[n]
    for n in range(2^(w-1) + 1, 2^w - 1 + 1, 2):
        Y[n] = ecY(P) * (W[n+2]*W2[n-1] - W[n-2]*W2[n+1])

    return dict((n, (X[n], Y[n], W[n])) for n in range(3, 2^w, 2))

#ScalarMult
def smultWidth(l, smult_alg):
    if smult_alg == SmallMultA:
        if l <= 256:
            return min(l-2, 4)    
        return 5
    if smult_alg == SmallMultJ:
        if l <= 256:
            return min(l-2, 5)    
        return 6
    return None

def digitize(n, base, k):
    while k:
        n, d = divmod(n, base)
        k -= 1
        yield d

def ScalarMult(E, P, d, smult_alg):
    d_even = d % 2
    d = (1 - d_even) * E.order() + (2*d_even - 1)*d  
    l = int(E.order()).bit_length()
    w = smultWidth(l, smult_alg)
    pre = smult_alg(E, P, w)
    pre[1] = tuple(P)
    
    for n in range(1, 2^w, 2):
        pre[-n] = ecNeg(pre[n])
    
    k = ceil(l / w)
    digits = list(digitize(d, 2^w, k))
    digits[k-1], digits[k-2] = digits[k-1] + is_even(digits[k-1]), digits[k-2] - is_even(digits[k-1])*2^w
    Q = pre[digits[k-1]]
    for n in range(k-2, 0, -1):
        digits[n], digits[n-1] = digits[n] + is_even(digits[n]), digits[n-1] - is_even(digits[n])*2^w
        for _ in range(w-1):
            Q = ecDblJ(E, Q)
        Q = ecDblAdd(E, Q, pre[digits[n]], smult_alg)
    
    for _ in range(w-1):
        Q = ecDblJ(E, Q)
    
    if is_odd(E.order() >> w):
        Q = ecDblJ(E, Q)
        Q = ecAdd_complete(E, Q, pre[digits[0]], smult_alg)
        Q = ecHToA(Q)
    else:
        Q = ecDblAdd(E, Q, pre[digits[0]], smult_alg)
        Q = ecJToA(Q)

    Q = ecChangeSign(Q, bool(d_even))
    return Q

def ecChangeSign(P, change_sign):
    X, Y, Z = P
    return (X, (-1)^(1-change_sign) * Y, Z)

def testSmallMult(E, P):
    w = 6
    smultA = SmallMultA(E, P, w, check_divp=True)
    #print(smultA)
    smultJ = SmallMultJ(E, P, w, check_divp=True)
    #print(smultJ)
    expected = dict((k, k*P) for k in range(3, 2^w, 2))
    #print(expected)
    print('SmallMultA: ', allEq(expected, smultA, ecEqPointAA))
    print('SmallMultJ: ', allEq(expected, smultJ, ecEqPointAJ))

def testScalarMult(E, P):    
    t = 60
    for j in range(1, t):
        print('ScalarMult[{d}, {smult_alg}]'.format(d=j, smult_alg='SmallMultA'), ecEqPointAA(j*P, ScalarMult(E, P, j, SmallMultA)))
        print('ScalarMult[{d}, {smult_alg}]'.format(d=j, smult_alg='SmallMultJ'), ecEqPointAA(j*P, ScalarMult(E, P, j, SmallMultJ)))
    for j in range(E.order() - t, E.order()):
        print('ScalarMult[{d}, {smult_alg}]'.format(d=j, smult_alg='SmallMultA'), ecEqPointAA(j*P, ScalarMult(E, P, j, SmallMultA)))
        print('ScalarMult[{d}, {smult_alg}]'.format(d=j, smult_alg='SmallMultJ'), ecEqPointAA(j*P, ScalarMult(E, P, j, SmallMultJ)))

def testCurve(E):    
    P = E.random_element()
    while P.is_zero():
        P = E.random_element() 
    print(P)
    testSmallMult(E, P)
    testScalarMult(E, P)

    
#STB  l = 128 curve
p_256 = 2^256 - 189

a_256 = 2^256 - 192

b_256 = int.from_bytes([
	0xF1, 0x03, 0x9C, 0xD6, 0x6B, 0x7D, 0x2E, 0xB2,
	0x53, 0x92, 0x8B, 0x97, 0x69, 0x50, 0xF5, 0x4C,
	0xBE, 0xFB, 0xD8, 0xE4, 0xAB, 0x3A, 0xC1, 0xD2,
	0xED, 0xA8, 0xF3, 0x15, 0x15, 0x6C, 0xCE, 0x77
], byteorder='little', signed=False)

q_256 = 2^256 - 51_359303463_308904523_350978545_619999225
    
#STB  l = 192 curve
p_384 = 2^384 - 317

a_384 = 2^384 - 320

b_384 = int.from_bytes([
    0x64, 0xBF, 0x73, 0x68, 0x23, 0xFC, 0xA7, 0xBC,
	0x7C, 0xBD, 0xCE, 0xF3, 0xF0, 0xE2, 0xBD, 0x14,
	0x3A, 0x2E, 0x71, 0xE9, 0xF9, 0x6A, 0x21, 0xA6,
	0x96, 0xB1, 0xFB, 0x0F, 0xBB, 0x48, 0x27, 0x71,
	0xD2, 0x34, 0x5D, 0x65, 0xAB, 0x5A, 0x07, 0x33,
	0x20, 0xEF, 0x9C, 0x95, 0xE1, 0xDF, 0x75, 0x3C
], byteorder='little', signed=False)

q_384 = 2^384 - 9886_438520659_958522437_788006980_660965037_549058207_958390857

#STB  l = 256 curve
p_512 = 2^512 - 569

a_512 = 2^512 - 572

b_512 = int.from_bytes([
	0x90, 0x9C, 0x13, 0xD6, 0x98, 0x69, 0x34, 0x09,
	0x7A, 0xA2, 0x49, 0x3A, 0x27, 0x22, 0x86, 0xEA,
	0x43, 0xA2, 0xAC, 0x87, 0x8C, 0x00, 0x33, 0x29,
	0x95, 0x5E, 0x24, 0xC4, 0xB5, 0xDC, 0x11, 0x27,
	0x88, 0xB0, 0xAD, 0xDA, 0xE3, 0x13, 0xCE, 0x17,
	0x51, 0x25, 0x5D, 0xDD, 0xEE, 0xA9, 0xC6, 0x5B,
	0x89, 0x58, 0xFD, 0x60, 0x6A, 0x5D, 0x8C, 0xD8,
	0x43, 0x8C, 0x3B, 0x93, 0x44, 0x59, 0xB4, 0x6C
], byteorder='little', signed=False)

q_512 = int.from_bytes([
	0xF1, 0x8E, 0x06, 0x0D, 0x49, 0xAD, 0xFF, 0xDC,
	0x32, 0xDF, 0x56, 0x95, 0xE5, 0xCA, 0x1B, 0x36,
	0xF4, 0x13, 0x21, 0x2E, 0xB0, 0xEB, 0x6B, 0xF2,
	0x4E, 0x00, 0x98, 0x01, 0x2C, 0x09, 0xC0, 0xB2,
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
	0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
], byteorder='little', signed=False)


curve_256 = EllipticCurve(Zmod(p_256), [a_256, b_256])
curve_256.set_order(q_256)

curve_384 = EllipticCurve(Zmod(p_384), [a_384, b_384])
curve_384.set_order(q_384)

curve_512 = EllipticCurve(Zmod(p_512), [a_512, b_512])
curve_512.set_order(q_512)

testCurve(curve_256)
testCurve(curve_384)
testCurve(curve_512)
